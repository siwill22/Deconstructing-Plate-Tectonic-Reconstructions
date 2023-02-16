import pygplates
import pygmt
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from gprm import SubductionConvergence, PointDistributionOnSphere
from gprm.utils import pmag
from gprm.datasets import Reconstructions
from gprm import MotionPathFeature
from gprm.utils.create_gpml import gpml2gdf, gdf2gpml
import helper_functions as hf
import copy
import itertools


DEFAULT_REGION = 'd'
DEFAULT_PROJECTION = 'W6i'
DEFAULT_FRAME = 'lrtb'


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def time_stamp(fig, reconstruction_time, x=0.25, y=0.06, font='20p'):
    fig.text(x=x,y=y,text='{:0.0f} Ma'.format(reconstruction_time),
                 region='0/1/0/1', projection='x5c', font=font, no_clip=True)


def paleogeography(fig, raster_series, reconstruction_time,
                   region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):
    
    reconstruction_time = find_nearest(list(raster_series.keys()), reconstruction_time)
    
    with pygmt.config(COLOR_NAN='white'):
        pygmt.makecpt(cmap='geo', series=[-7000,5000,200], background='o')
        fig.grdimage(raster_series[reconstruction_time], 
                    region=region, projection=projection, cmap=True)
        
    fig.basemap(region=region, projection=projection, frame=frame)
    
    
    
def continental_drift(fig, reconstruction_model, reconstruction_time,
                      region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):
    
    if reconstruction_time==0:
        reconstruction_time=0.001

    if isinstance(reconstruction_model.continent_polygons, pygplates.FeatureCollection):
        cloned_polygons = reconstruction_model.continent_polygons.clone()
    else:
        cloned_polygons = pygplates.FeatureCollection(
            list(itertools.chain([list(p) for p in reconstruction_model.continent_polygons]))[0]
        ).clone()
    saved_polygons = reconstruction_model.continent_polygons
    reconstruction_model.continent_polygons = cloned_polygons

    reconstructed_continents = reconstruction_model.polygon_snapshot('continents', reconstruction_time)
    #reconstructed_continents = reconstruction_model.reconstruct(cloned_polygons, reconstruction_time)
    reconstructed_plates = reconstruction_model.plate_snapshot(reconstruction_time)

    # Note this is a trick to allow the colors of the polygons
    # to be set in a way that is consistent with the topological polygons,
    # even though the topologies are not visible
    plate_partitioner = pygplates.PlatePartitioner(reconstructed_plates.resolved_topologies, reconstruction_model.rotation_model)
    for p in reconstructed_continents.reconstructed_polygons:
        result = plate_partitioner.partition_point(p.get_reconstructed_geometry().get_interior_centroid())
        if result is None:
            p.get_feature().set_reconstruction_plate_id(0)
        else:
            p.get_feature().set_reconstruction_plate_id(result.get_resolved_feature().get_reconstruction_plate_id()) #, verify_information_model=pygplates.VerifyInformationModel.no)
    
    fig.basemap(region=region, projection=projection, frame=frame)


    reconstructed_continents.plot(fig, color='+z', cmap='./plate_id_regular.cpt',
                                  aspatial='Z=PLATEID1', pen='0.2p,grey50',
                                  transparency=10, region=region, projection=projection)
    
    reconstruction_model.continent_polygons = saved_polygons

    #fig.text(x=0.0,y=1.05,text='(a) Continental Drift', justify='LM',
    #             region='0/1/0/1', projection='x5c', font='14p', no_clip=True)
    fig.basemap(region=region, projection=projection, frame=frame)


def topological_model(fig, reconstruction_model, reconstruction_time,
                      region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME, 
                      legend=False):
    
    if reconstruction_time==0:
        reconstruction_time=0.001
    
    reconstructed_continents = reconstruction_model.polygon_snapshot('continents', reconstruction_time)
    reconstructed_plates = reconstruction_model.plate_snapshot(reconstruction_time)

    velocity_domain = PointDistributionOnSphere(distribution_type='healpix', N=8)
    velocity_field = reconstructed_plates.velocity_field(velocity_domain_features=[velocity_domain.meshnode_feature])

    fig.basemap(region='d', projection=projection, frame=frame)
    
    reconstructed_continents.plot(fig, color='gray60', region=region, projection=projection)

    reconstructed_plates.plot_polygons(fig, color='+z', cmap='./plate_id_regular.cpt',#cmap=True, reduce_plate_ids=True,
                                       aspatial='Z=PLATEID1', close=True,
                                       transparency=55, region=region, projection=projection)
    
    reconstructed_plates.plot_subduction_zones(fig, pen='0.75p,black', gap=8, size=4, region=region, projection=projection, label='Subduction Zones')
    reconstructed_plates.plot_mid_ocean_ridges(fig, pen='0.75p,red', region=region, projection=projection, label='Mid-Ocean Ridges')
    reconstructed_plates.plot_other_boundaries(fig, pen='0.75p,gray', region=region, projection=projection)
    velocity_field.plot(fig, style='V0.08c+e+a45+ggray40+n', scaling=300., pen="0.5p,black", color="gray", 
                        transparency=30, region=region, projection=projection)

    fig.basemap(region=region, projection=projection, frame=frame)
    
    #fig.legend(transparency=20, position='JBR+jBR+o-1.5c/-0.4c', box='+gwhite+p0.5p')
    if legend:
        fig.legend(transparency=20, position='JBR+jBR+o0.0c/-0.4c', box='+gwhite+p0.4p')

    #fig.text(x=0.0,y=1.05,text='(b) Plate Tectonics', justify='LM',
    #             region='0/1/0/1', projection='x5c', font='14p', no_clip=True)
    


def deforming_model(fig, reconstruction_model, reconstruction_time, time_spans,
                    region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME, 
                    legend=False):
    
    if reconstruction_time==0:
        reconstruction_time=1
    
    crustal_thickness_delta = hf.get_deltas(time_spans, reconstruction_time, delta_time=-1)

    reconstructed_continents = reconstruction_model.polygon_snapshot('continents', reconstruction_time)
    reconstructed_plates = reconstruction_model.plate_snapshot(reconstruction_time)

    velocity_domain = PointDistributionOnSphere(distribution_type='healpix', N=8)
    velocity_field = reconstructed_plates.velocity_field(velocity_domain_features=[velocity_domain.meshnode_feature])


    fig.basemap(region=region, projection=projection, frame=frame)

    with pygmt.config(COLOR_NAN='white'):
        pygmt.makecpt(cmap='polar', series=[-1000,1000], background='o')
        fig.grdimage(crustal_thickness_delta, 
                    region=region, projection=projection, cmap=True)

    reconstructed_continents.plot(fig, color='gray70', pen='0.2p,gray50', projection=projection,
                                  region=region, transparency=65)
    
    velocity_field.plot(fig, scaling=300., pen="0.5p,grey", color="grey", 
                        style='V0.08c+e+a45+ggray40+n', region=region, projection=projection,
                        transparency=20)
    reconstructed_plates.plot_deformation_zones(fig, color='-', pen='0.2p,grey', region=region, projection=projection)
    reconstructed_plates.plot_subduction_zones(fig, region=region, projection=projection)
    reconstructed_plates.plot_mid_ocean_ridges(fig, region=region, projection=projection)
    reconstructed_plates.plot_other_boundaries(fig, region=region, projection=projection)
    
    if legend:
        with pygmt.config(FONT_ANNOT_PRIMARY='1p,white', FONT_LABEL='10p', MAP_TICK_LENGTH_PRIMARY='2p', MAP_FRAME_PEN='1p,black'):
            fig.colorbar(position='JBR+jBR+o0.2c/-0.2c+w4.0c/0.33c+h', frame=['+n','xa1000+lContraction             Extension '])

    #fig.text(x=0.0,y=1.05,text='(d) Deforming Plates', justify='LM',
    #             region='0/1/0/1', projection='x5c', font='14p', no_clip=True)
    
    
def optimization(fig, reconstruction_model, reconstruction_time,
                 region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME, 
                 legend=False):
    
    if reconstruction_time==0:
        reconstruction_time=0.001
    
    subduction_kinematics = SubductionConvergence(reconstruction_model, reconstruction_time, threshold_sampling_distance_radians=0.02)
    mr = np.asarray(subduction_kinematics.df['migr_rate'], dtype=np.float64)
    mo = np.asarray(subduction_kinematics.df['migr_obliq'], dtype=np.float64)
    subduction_kinematics.df['ortho_migr_rate'] = pd.Series(mr*np.sin(np.radians(np.abs(mo))), index=subduction_kinematics.df.index)

    reconstructed_plates = reconstruction_model.plate_snapshot(reconstruction_time)


    velocity_domain = PointDistributionOnSphere(distribution_type='healpix', N=32)
    velocity_field = reconstructed_plates.velocity_field(velocity_domain_features=[velocity_domain.meshnode_feature])
    velgrd = velocity_field.to_grid(spacing='0.25d')

    velocity_domain = PointDistributionOnSphere(distribution_type='healpix', N=8)
    velocity_field = reconstructed_plates.velocity_field(velocity_domain_features=[velocity_domain.meshnode_feature])


    scaling=10
    tmp = np.vstack((subduction_kinematics.df.lon,
                     subduction_kinematics.df.lat,
                     subduction_kinematics.df.arc_azimuth+180,
                     subduction_kinematics.df.ortho_migr_rate/scaling)).T

    fig.basemap(region=region, projection=projection, frame=frame)

    pygmt.makecpt(cmap='magma', series=[0,100,10], background='o', reverse=False, transparency=50)
    fig.grdimage(grid=velgrd, cmap=True, transparency=75)

    reconstruction_model.polygon_snapshot('continents', reconstruction_time).plot(fig, transparency=70, pen=None, color='gray20')

    reconstructed_plates.plot_mid_ocean_ridges(fig, pen='1p,white', transparency=10)
    reconstructed_plates.plot_other_boundaries(fig, pen='0.7p,gray50', transparency=40)

    velocity_field.plot(fig, scaling=300., pen="0.7p,gray40", color='gray40',
                        style='V0.15c+e+a45+ggray40+n', transparency=30)

    if legend:
        with pygmt.config(FONT_ANNOT_PRIMARY='8p', FONT_LABEL='10p'):
            fig.colorbar(position='JBL+jBL+o0.0c+w3.5c/0.3c+h', transparency=20, frame=['x+lPlate Velocity [mm/yr]'])#, box='+gwhite+p0.5p')

    fig.plot(x=subduction_kinematics.df.lon,
             y=subduction_kinematics.df.lat,
             style='c0.33c', color='gray20', transparency=30)

    pygmt.makecpt(cmap='polar', series=[-60,60,10], background='o', reverse=True)
    fig.plot(x=subduction_kinematics.df.lon,
             y=subduction_kinematics.df.lat,
             color=subduction_kinematics.df.ortho_migr_rate*10,
             style='c0.25c', cmap=True, transparency=50)
    if legend:
        with pygmt.config(FONT_ANNOT_PRIMARY='8p', FONT_LABEL='10p'):
            fig.colorbar(position='JBR+jBR+o0.0c+w3.5c/0.3c+h', frame=['x+lTrench Migration [mm/yr]'])#, box='+gwhite+p0.5p')

    reconstructed_plates.plot_subduction_zones(fig, pen='1p,black', gap=8, size=3)

    
    
def polygons_colored_by_data_density(fig, polygons, lookup_table,
                                     reconstruction_model, reconstruction_time, anchor_plate_id=1):
    
    tmp = polygons[polygons.PLATEID1.isin(list(lookup_table.keys()))].reset_index(drop=True)
    tmp = reconstruction_model.reconstruct(tmp, reconstruction_time, anchor_plate_id=anchor_plate_id)

    tmp['counts'] = list(lookup_table[tmp.PLATEID1])

    tmp = tmp[['geometry', 'counts']]
    pygplates.FeatureCollection(gdf2gpml(tmp)).write('tmp.gmt')
    fig.plot(data='tmp.gmt', aspatial='Z=counts', color='+z', cmap=True, transparency=50)
    


#for reconstruction_time in np.arange(500,0,-animation_time_step):
def apwp_reconstruction(fig, reconstruction_model, reconstruction_time, vgps_gdf, 
                        anchor_plate_id=1, reference_plate_id=701, color_polygons=None,
                        apwp_min_time=0, apwp_max_time=540, apwp_time_step=10,
                        region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):
    
    path_times=np.arange(apwp_min_time, apwp_max_time+apwp_time_step, apwp_time_step)
    rGAPWaP = MotionPathFeature(seed_points=(-90, 0), 
                               path_times=path_times, 
                               reconstruction_plate_id=anchor_plate_id, 
                               relative_plate_id=reference_plate_id).reconstruct_motion_path(reconstruction_model,
                                                                              reconstruction_time=reconstruction_time,
                                                                              anchor_plate_id=anchor_plate_id)
    if len(rGAPWaP)>0:
        rGAPWaP = rGAPWaP[0]
        rGAPWaP = gpd.GeoDataFrame(pd.DataFrame(data={'Age':path_times[::-1][:rGAPWaP.shape[0]]}), 
                                   geometry=gpd.points_from_xy(rGAPWaP[:,1], rGAPWaP[:,0]), crs=4326)
    else:
        rGAPWaP = None


    rpoles = pmag.rotate_to_common_reference(vgps_gdf.copy(), reconstruction_model, reference_plate_id=reference_plate_id)
    rpoles['_PLATEID1'] = rpoles['PLATEID1'] # presernve the original plateids
    rpoles['PLATEID1'] = reference_plate_id    # set the 'in use' plateid to be 701 
    rpoles = rpoles.set_geometry(gpd.points_from_xy(rpoles.PoleLongitude, rpoles.PoleLatitude))

    rrpoles = reconstruction_model.reconstruct(rpoles.query('AverageAge >= @reconstruction_time'), 
                                               reconstruction_time=reconstruction_time,
                                               anchor_plate_id=anchor_plate_id)

    
    fig.basemap(region=region, projection=projection, frame=frame)
    reconstruction_model.polygon_snapshot(polygon_type='continents', 
                                          reconstruction_time=reconstruction_time,
                                          anchor_plate_id=anchor_plate_id).plot(fig, color='grey90')
    
    
    pygmt.makecpt(cmap='oslo', series=[0,600,25], reverse=True)
    fig.plot(x=rrpoles.geometry.x, y=rrpoles.geometry.y, size=rrpoles.PoleA95*111.,
             style='E-', 
             #pen='1.0p,red', 
             color=rrpoles.AverageAge, 
             cmap=True, 
             transparency=85)
    
    rrpoles = rrpoles.query('AverageAge <= @reconstruction_time+@apwp_time_step')
    if len(rrpoles)>0:
        #print(rrpoles._PLATEID1.value_counts())
        fig.plot(x=rrpoles.geometry.x, y=rrpoles.geometry.y, size=rrpoles.PoleA95*111.,
                 style='E-', 
                 pen='5.0p,blue', 
                 color=rrpoles.AverageAge, 
                 cmap=True, 
                 transparency=80)
        if color_polygons is not None:
            counts = rrpoles._PLATEID1.value_counts()
            pygmt.makecpt(cmap='yellow,darkorange,darkred', series=[1,3], background='o')
            polygons_colored_by_data_density(fig, color_polygons, counts, 
                                             reconstruction_model, reconstruction_time, 
                                             anchor_plate_id=anchor_plate_id)

    if rGAPWaP is not None:
        pygmt.makecpt(cmap='oslo', series=[0,600,25], reverse=True)
        fig.plot(x=rGAPWaP.geometry.x, y=rGAPWaP.geometry.y, pen='1p,darkblue', label='Modelled Polar Wander Path for Southern Africa')
        fig.plot(x=rGAPWaP.geometry.x, y=rGAPWaP.geometry.y, style='c0.25c', color=rGAPWaP.Age, cmap=True, pen='0.25p,gray20')

    
    
    
def add_depth_slice(fig, slice_depth, cmapstring, projection, region, contour=0.25):

    depth_slice = '/Users/simon/Data/SeismicTomography/DETOX_MODELS/DETOX-P3/grid_nc4/DETOX-P3_{:0.0f}.00.grd'.format(slice_depth)
    depth_slice_filt = pygmt.grdfilter(grid=depth_slice, distance=2, filter='g250', coltypes='g')
    depth_slice_filt.data[depth_slice_filt.data<contour] = np.nan
    depth_slice_filt.data[depth_slice_filt.data>contour] = contour
    pygmt.makecpt(cmap=cmapstring, series=[contour,contour*2,contour], reverse=True)#, background='o')
    fig.grdimage(grid=depth_slice_filt, 
                 transparency=10, region=region, projection=projection, cmap=True, nan_transparent=True)
    

def tomotectonic_model(fig, reconstruction_model, reconstruction_time,
                       region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME, 
                       colorbar=False):
    
    # ???????DETOX models use a simple depth (in km) = time (in Myr) * 10
    slice_depth = np.round(reconstruction_time/5)*5*10.
    print(slice_depth)
    
    #pygmt.makecpt(cmap='devon', series=[0,0.6,0.2], background='o', reverse=True)
    #fig.grdimage(grid='/Users/simon/Data/SeismicTomography/DETOX_MODELS/DETOX-P3/grid_nc4/DETOX-P3_{:0.0f}.00.grd'.format(slice_depth),
    #             cmap=True)
    fig.basemap(region=region, projection=projection, frame=frame)

    contour = 0.3
    depth_slice_interval=150
    add_depth_slice(fig, slice_depth+(depth_slice_interval*3), 'purple3,purple3', projection, region, contour)
    add_depth_slice(fig, slice_depth+(depth_slice_interval*2), 'mediumblue,mediumblue', projection, region, contour)
    add_depth_slice(fig, slice_depth+depth_slice_interval, 'royalblue,royalblue', projection, region, contour)
    add_depth_slice(fig, slice_depth, 'slateblue1,slateblue1', projection, region, contour) 
    
    reconstruction_model.polygon_snapshot('coastlines', reconstruction_time).plot(fig, color='gray50', 
                                                                                  pen='1p,grey70',
                                                                                  transparency=70,
                                                                                  projection=projection, region=region)
    reconstruction_model.plate_snapshot(
        reconstruction_time).plot_subduction_zones(fig, pen='1.25p,black', gap=10, size=3, color='black', transparency=20,
                                                   projection=projection, region=region)

    if colorbar:
        pygmt.makecpt(categorical=True, cmap='slateblue1,royalblue,mediumblue,purple3', 
                      color_model='+c{:0.0f},{:0.0f},{:0.0f},{:0.0f}'.format(slice_depth,
                                                                             slice_depth+depth_slice_interval,
                                                                             slice_depth+(depth_slice_interval*2),
                                                                             slice_depth+(depth_slice_interval*3)))
        with pygmt.config(FONT_ANNOT_PRIMARY='9p,white', FONT_LABEL='10p', MAP_TICK_LENGTH_PRIMARY='2p', MAP_FRAME_PEN='1p,black'):
            fig.colorbar(position='JBR+jBR+o0.2c/-0.2c+w4.5c/0.33c+h', frame=['+n','xa1+lDepth of Tomography Slice [km]'])
        with pygmt.config(FONT_LABEL='10p', MAP_TICK_LENGTH_PRIMARY='2p', MAP_FRAME_PEN='1p,black'):
            fig.colorbar(position='JBR+jBR+o0.2c/-0.2c+w4.5c/0.33c+h')
        
        
        
def climate_model(fig, dataset, variable, reconstruction_time, background_grids=None,
                  region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):
    
    reconstruction_time = find_nearest(dataset.simulation.data*10., reconstruction_time)
    tmp = dataset.sel(simulation=54-int(reconstruction_time/10)).mean(dim='month')[variable].drop_vars(['simulation'])
    tmp.to_netcdf('tmp.nc')

    if background_grids is not None:
        pygmt.makecpt(cmap='gray', series=[-5000,5000], background='o')
        fig.grdimage(grid=pygmt.grdclip(background_grids[reconstruction_time], below='-2500/-2000'), 
                     region=region, projection=projection, cmap=True, shading='+d')
    if variable=='T':
        pygmt.makecpt(cmap='polar', series=[-10,35,5], background='o')
    elif variable=='P':
        pygmt.makecpt(cmap='viridis', series=[0,300,50], background='o', reverse=True)
    else:
        print('Unknown Variable {:s}'.format(variable))
        return
    fig.grdimage(grid='tmp.nc', cmap=True, region=region, projection=projection, transparency=50)
    
    fig.basemap(region=region, projection=projection, frame=frame)

        
def mantle_temperature(fig, reconstruction_time, reconstruction_model, path_to_grid_files, 
                       slice_depth=2840, 
                       region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):

    pygmt.makecpt(cmap='polar', series=[-500,500,100], background='o')
    fig.grdimage(grid='{:s}/OPT1-temp-{:03d}Ma-{:04d}km_mean_removed_Dim.grd'
                 .format(path_to_grid_files, int(np.round(reconstruction_time/20)*20), int(slice_depth)), 
                 cmap=True, region=region, projection=projection)
    plate_snapshot = reconstruction_model.plate_snapshot(reconstruction_time).plot_boundaries(fig)
    fig.basemap(region=region, projection=projection, frame=frame)
    
    
def plate_hierarchy(fig, tree_object, reconstruction_model, reconstruction_time, 
                    tree_polygons='static', show_topologies=True, legend=False, 
                    region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):
    
    fig.basemap(region=region, projection=projection, frame=frame)
    
    reconstruction_model.polygon_snapshot('continents', reconstruction_time).plot(fig, color='gray60', 
                                                                                  pen='0.2p,gray80', transparency=30)
    if show_topologies:
        reconstructed_plates = reconstruction_model.plate_snapshot(reconstruction_time).plot_polygons(fig, pen='1p,gray50', 
                                                                                                      transparency=60)

    tree_object.plot_gmt(fig, reconstruction_time, polygons=tree_polygons,
                         node_color='dodgerblue', root_node_style='a0.5c')

    if legend:
        fig.legend(transparency=20, position='JBR+jBR+o-0.1c/-0.3c', box='+gwhite+p0.5p')
    
    
    
def LLSVP(fig, reconstruction_model, reconstruction_time, LIP_centroids=None, kimberlites=None,
          region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):
    
    fig.basemap(region=region, projection=projection, frame=frame)

    lekic_filt = pygmt.grdfilter(grid='../data/Clustering_Lekic.nc', distance=2, filter='g800', spacing='0.5d')
    pygmt.makecpt(cmap='bilbao', series=[-0.5,5.5,1], background=False) #, color_model="+c" + ",".join(['0','1','2','3','4','5']))
    fig.grdimage(grid=lekic_filt.astype(int), cmap=True, transparency=0, interpolation='l')

    reconstructed_continents = reconstruction_model.polygon_snapshot('continents', 
                                                                      reconstruction_time).plot(fig, 
                                                                                                color='grey20',
                                                                                                transparency=80)

    if LIP_centroids is not None:
        LIPs_at_birth_location = reconstruction_model.reconstruct_to_time_of_appearance(LIP_centroids.clone(), ReconstructTime='MidTime')
        LIPs_at_birth_location = gpml2gdf(LIPs_at_birth_location).query('FROMAGE-5 > @reconstruction_time')
        rLIPs = reconstruction_model.reconstruct(LIP_centroids, reconstruction_time)

        if len(LIPs_at_birth_location)>0:
            fig.plot(data=LIPs_at_birth_location, style='d0.3c', pen='0.7p,black', color='dodgerblue', 
                     transparency=60, label='Reconstructed LIPs')
        if len(rLIPs)>0:
            fig.plot(data=gpml2gdf(rLIPs), style='d0.45c', pen='0.7p,black', color='dodgerblue', 
                     transparency=15, label='Reconstructed LIPs')
        
    if kimberlites is not None:
        kimberlites_at_birth_location = reconstruction_model.reconstruct_to_time_of_appearance(kimberlites.clone(), ReconstructTime='MidTime')
        kimberlites_at_birth_location = gpml2gdf(kimberlites_at_birth_location).query('FROMAGE-5 > @reconstruction_time')
        rkimberlites = reconstruction_model.reconstruct(kimberlites, reconstruction_time)

        if len(kimberlites_at_birth_location)>0:
            fig.plot(data=kimberlites_at_birth_location, style='c0.1c', pen='0.7p,black', color='slateblue', 
                     transparency=90, label='Reconstructed kimberlites')
        if len(rkimberlites)>0:
            fig.plot(data=gpml2gdf(rkimberlites), style='c0.15c', pen='0.7p,black', color='slateblue', 
                     transparency=15, label='Reconstructed kimberlites')
            
            
            
def seafloor_age(fig, reconstruction_model, reconstruction_time, grid_filename_template, cmap='inferno',
                 region=DEFAULT_REGION, projection=DEFAULT_PROJECTION, frame=DEFAULT_FRAME):
    
    if reconstruction_time==0:
        reconstruction_time=0.001
        
    reconstructed_coastlines = reconstruction_model.polygon_snapshot('coastlines', reconstruction_time)
    reconstructed_plates = reconstruction_model.plate_snapshot(reconstruction_time)

    fig.basemap(region=region, projection=projection, frame=frame)

    pygmt.makecpt(cmap=cmap, series=[0,200,10], background='o', reverse=True)
    
    fig.grdimage(grid=grid_filename_template.format(reconstruction_time), cmap=True, interpolation='b')

    reconstructed_coastlines.plot(fig, color='gray40', pen='0p,gray40')
    reconstructed_plates.plot_subduction_zones(fig)
    reconstructed_plates.plot_mid_ocean_ridges(fig)
    reconstructed_plates.plot_other_boundaries(fig)

    