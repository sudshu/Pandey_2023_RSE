# Methane Plume Detection 
Code for methane plume detection and quantification using Sentinel-2, Sentinel-3, and Landsat instruments

## Methane_MSI

The Methane_MSI folder includes code for both Sentinel-2 and Landsat. This code retrieves satellite data from Google Earth Engine, which requires configuring an Earth Engine Python API.
https://developers.google.com/earth-engine/tutorials/community/intro-to-python-api

- `runner_notebook.ipynb` shows how to run a single-day retrieval, a time series, and the lookup-table tools.
- `scripts/MS_XCH4_retrieval.py`: Sentinel-2 / Landsat retrieval. `ch4ret.optimizer_v2` converts the delR image to a methane column enhancement (mol/m2) with the lookup table `scripts/data_files/test_srf_210104_delr_to_omega.pkl`.
- `scripts/Lookup_tables_generator.py`: absorption cross sections, radiative transfer, and lookup-table generation.
- `scripts/data_files/`: HITRAN line lists, atmosphere profiles (midlatitude winter/summer, subarctic winter/summer, tropical, US standard, Europe background), solar spectra, absorption cross-section csv files, and lookup tables. See `scripts/data_files/README.md` for the format and units of each file.

Both scripts find `data_files/` relative to their own location, so they can be run or imported from any working directory.

## Sentinel-3 
The Sentinel-3 folder handles the downloading and processing of Sentinel-3 SWIR data. To obtain Sentinel-3 SLSTR observations, the Sentinelsat API https://scihub.copernicus.eu/dhus/#/home must be configured locally by the user.




## Cite as: 
Pandey, S., van Nistelrooij, M., Maasakkers, J. D., Sutar, P., Houweling, S., Varon, D. J., Tol, P., Gains, D., Worden, J., and Aben, I.: Daily detection and quantification of methane leaks using Sentinel-3: a tiered satellite observation approach with Sentinel-2 and Sentinel-5p, Remote Sens. Environ., 296, 113716, https://doi.org/10.1016/j.rse.2023.113716, 2023.



