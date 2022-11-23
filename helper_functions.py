import pygplates
import pygmt
import numpy as np
from gprm import PointDistributionOnSphere
import geopandas as gpd
import pandas as pd

from gprm.utils.deformation import get_crustal_thickness_points, raster_topological_reconstruction



def get_topological_model(reconstruction_model, anchor_plate_id=0):
    topological_model = pygplates.TopologicalModel(
        reconstruction_model.dynamic_polygons,
        reconstruction_model.rotation_model,
        anchor_plate_id=anchor_plate_id,
        # Enable strain rate clamping to better control crustal stretching factors...
        default_resolve_topology_parameters=pygplates.ResolveTopologyParameters(enable_strain_rate_clamping=True))
    
    return topological_model


def get_time_spans(topological_model, initial_time, oldest_time, youngest_time, N=128):
    
    pts = PointDistributionOnSphere(distribution_type='healpix', N=N)
    pts_list = [p for p in pts.multipoint.get_points()]
    reconstructed_crustal_thickness = get_crustal_thickness_points(pts_list)

    time_spans = topological_model.reconstruct_geometry(
        pts_list,
        initial_time=initial_time,
        oldest_time=oldest_time,
        youngest_time=youngest_time,
        initial_scalars={pygplates.ScalarType.gpml_crustal_thickness : reconstructed_crustal_thickness},
        # All our points are on continental crust so we keep them active through time (ie, never deactivate them)...
        deactivate_points = pygplates.ReconstructedGeometryTimeSpan.DefaultDeactivatePoints(
            deactivate_points_that_fall_outside_a_network = False))
    
    return time_spans


def get_crustal_thickness(time_spans, reconstruction_time, return_inactive_points=False,
                          region='d', spacing='0.5d', search_radius='5d', gridding='surface'):
    
    d0 = get_scalars(time_spans, reconstruction_time, return_inactive_points=return_inactive_points)
    d0 = d0.dropna()

    if gridding=='nearneighbor':
        res = pygmt.nearneighbor(
            pygmt.blockmean(data=d0[['x','y','ct']], 
                            region=region, 
                            spacing='{:f}d'.format(float(spacing[:-1])/2)), #TODO this spacing should be smaller than final grid spacing?? 
            region=region, 
            spacing=spacing, 
            search_radius=search_radius)
    elif gridding=='surface':
        res = pygmt.surface(
            pygmt.blockmean(data=d0[['x','y','ct']], 
                            region=region, 
                            spacing='{:f}d'.format(float(spacing[:-1])/2)), #TODO this spacing should be smaller than final grid spacing?? 
            region=region, 
            spacing=spacing, 
            T=0.4, M=search_radius)
    
    return res
    
    
def get_scalars(time_spans, reconstruction_time, return_inactive_points=False):
    # given a time_spans object, return a dataframe of reconstructed crustal thickness
    
    # Note that here we would need to keep return_inactive_points as True to ensure that the length of points remains the same
    # between the two times
    reconstructed_points = time_spans.get_geometry_points(reconstruction_time, return_inactive_points=return_inactive_points)
    crustal_thickness = time_spans.get_scalar_values(reconstruction_time, 
                                                     return_inactive_points=return_inactive_points)[pygplates.ScalarType.gpml_crustal_thickness]

    # convert reconstructing crustal thickness points to dataframe, invalid points are preserved as nans 
    d = [[rp.to_lat_lon()[1], rp.to_lat_lon()[0], ct] if rp is not None else [-999,-999,-999] for rp,ct in zip(reconstructed_points, crustal_thickness)]

    df = pd.DataFrame(data=np.vstack(d), columns=['x', 'y', 'ct'])
    df = df.replace({-999: np.nan})
    
    return df


def get_deltas(time_spans, reconstruction_time, delta_time=-1, region='d', spacing='0.5d', search_radius='5d'):
    
    d0 = get_scalars(time_spans, reconstruction_time, return_inactive_points=True)
    d1 = get_scalars(time_spans, reconstruction_time+delta_time, return_inactive_points=True)

    d0['delta'] = d0.ct-d1.ct
    d0 = d0[d0['delta'].abs()>1]
    
    res = pygmt.nearneighbor(
        pygmt.blockmean(data=d0[['x','y','delta']], 
                        region=region, 
                        spacing='{:f}d'.format(float(spacing[:-1])/2)), #TODO this spacing should be smaller than final grid spacing?? 
        region=region, 
        spacing=spacing, 
        search_radius=search_radius)
    
    return res
    
    