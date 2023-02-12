# Deconstructing Plate Tectonic Reconstructions

This repository contains code to reproduce figures and animations associated with the manuscript **Deconstructing plate tectonic reconstructions**, authors Seton, Williams, Domeier, Collins & Sigloch, to appear in *Nature Reviews Earth & Environment*

### Citation
```
@article{Seton++2023DeconstructingPlateTectonicReconstructions,
  author  = {Maria Seton and Simon Williams and Mathew Domeier and Karin Sigloch and Alan Collins},
  title   = {Deconstructing Plate Tectonic Reconstructions},
  journal = {Nature Reviews Earth and Environment},
  year    = {2023},
  volume  = {},
  number  = {},
  pages   = {},
  doi     = {TBA},
}
```

## Requirements
- all code is based on python3, and run through jupyter notebooks. 
- Most python dependencies are installable via conda or pip; additional requirements [pygplates](https://www.gplates.org/docs/pygplates/), [PlateTectonicsTools](https://github.com/EarthByte/PlateTectonicTools) and [GPlatesReconstructionModel](https://github.com/siwill22/GPlatesReconstructionModel) must currently be configured separately - follow the links for instructions.
- The code relies on a wide range of independently compiled data sources and plate tectonic reconstructions. Many of these data sources are automatically downloaded and saved to a local cache when the code is run for the first time. Other data are contained and described in the [Data folder](../main/data).


## Notes
- The illustrations in the paper were created using either [GPlates](http://www.gplates.org) or through python code. This repository contains the code for the latter of these groups.
- The code is organised into folders corresponding to each figure. 
- Figures produced by this code will be similar to, but not identical to, the figures in the paper due to formatting changes made by the publishers.

