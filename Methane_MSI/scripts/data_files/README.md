# data_files

Input and output files used by `Lookup_tables_generator.py` and `MS_XCH4_retrieval.py`.
Both scripts locate this folder relative to their own location, so they can be run from any working directory.

## Inputs

| File | Content |
| ---- | ------- |
| `linelist2016_1600nm.dat`, `linelist2016_2100nm.dat` | HITRAN 2016 line parameters for H2O, CO2, N2O, CO and CH4 in the SWIR1 (~1.6 µm) and SWIR2 (~2.2 µm) regions. Column meanings are listed in the file header and in `ch4ret.abCalc`. |
| `atmosphere_*.dat` | Atmosphere profiles. One row per layer, starting at the surface. Columns: layer index, pressure [hPa], temperature [K], then layer column number densities [molecules cm^-2] of H2O, CO2, O3, N2O, CO, CH4 and O2 (see file header). |
| `solar_irradiance_1600nm_highres_extended_sparse.dat`, `solar_irradiance_2100nm_highres_sparse.dat` | High-resolution solar spectrum: wavenumber [cm^-1] and solar irradiance on a 0.05 cm^-1 grid. |
| `srf_s2_l8_s2.pkl` | SWIR1 and SWIR2 spectral response functions for S2, L8 and S3 (not read by the current scripts). |

Available atmosphere profiles: `atmosphere_midlatitudewinter.dat` (default), `atmosphere_midlatitudesummer.dat`,
`atmosphere_subarcticsummer.dat`, `atmosphere_subarcticwinter.dat`, `atmosphere_tropical.dat`, `atmosphere_standard.dat`
(AFGL / US Standard, 24 layers) and `atmosphere_europebackground.dat` (CAMELOT Europe background, 23 layers).
Select one with the `atmfile` argument of `ch4ret` or `CreateAbsorptionCrossSections`.

## Absorption cross sections: `<sat>_absorption_cs_<gas>_<band>.csv`

Written by `ch4ret.abCalc` (called from `CreateAbsorptionCrossSections`) and read by `radianceCalc`.

- `<sat>`: `S2`, `S3` or `L8`; `<gas>`: `H2O`, `CO2`, `N2O`, `CO`, `CH4`; `<band>`: `SWIR1` or `SWIR2`.
- Plain comma-separated matrix, no header.
- **Rows**: spectral grid points, in order of increasing wavenumber (decreasing wavelength), spacing 0.05 cm^-1, covering the band limits set in `abCalc` for that satellite. Row `i` is at wavenumber `w[i] = round(1e7/lambda_u, 1) + 0.05*i` cm^-1, i.e. wavelength `1e7/w[i]` nm.
- **Columns**: atmosphere layers, in the same order as the rows of the atmosphere file used (column 0 = surface layer).
- **Values**: absorption cross section in cm^2 per molecule, at the pressure and temperature of that layer.

The vertical optical depth is the cross section times the layer column density, summed over layers and gases:
`tau_vert = sum_gas sigma_gas @ n_gas`. The files shipped in this repository are for S3 with 24-layer profiles.
The file names do not include the atmosphere profile, so regenerate them (`CreateAbsorptionCrossSections(sat, atmfile)`)
with the same profile you use in `radianceCalc`.

Quick look:

```python
import numpy as np, matplotlib.pyplot as plt
sigma = np.genfromtxt("S3_absorption_cs_CH4_SWIR2.csv", delimiter=",")   # (n_wavenumbers, n_layers)
w = np.round(1e7/(2255.70 + 50.15/2), 1) + 0.05*np.arange(sigma.shape[0])  # S3 SWIR2 grid [cm^-1]
plt.semilogy(1e7/w, sigma[:, 0]); plt.xlabel("Wavelength (nm)"); plt.ylabel("CH4 cross section, surface layer (cm$^2$/molecule)")
```

`Lookup_tables_generator.plotOpticalDepths(sat="S3", band="SWIR2")` plots the resulting optical depths of H2O, CO2 and CH4.

## Lookup tables

| File | Content |
| ---- | ------- |
| `test_srf_210104_delr_to_omega.pkl` | Used by `MS_XCH4_retrieval.fullMBMP2Omega` and `ch4ret.optimizer_v2`. Dictionary with `omegas` (methane enhancements, mol m^-2) and, for `S2`, `L8` and `S3`, band-integrated `SWIR1` and `SWIR2` radiances for each total air mass factor (keys `'2.0'` to `'4.9'`). |
| `<sat>_full_mdata_poly_10_delr_to_omega.pkl` | Written by `createLookupTables`: `data[sat]['MBMP'][amf]` holds degree-10 polynomial coefficients mapping `log(delR + 1)` to methane enhancement (mol m^-2). |
| `L8_210104_delr_to_omega.pkl` | Same format as `test_srf_210104_delr_to_omega.pkl`, L8 only. (`Sentinel-3/test_srf_210104_delr_to_omega.pkl` is an identical copy of the S2/L8/S3 table.) |
