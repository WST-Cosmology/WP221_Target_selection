# Target Selection Strategies for the Wide-Field Spectroscopic Telescope (WST)

This repository compiles and documents the different **target selection strategies** for the main cosmological tracers:

- Bright Galaxy Survey (BGS)  
- Emission Line Galaxies (ELG)  
- Luminous Red Galaxies (LRG)  
- Quasi-Stellar Objects (QSO)  
- Lyman-Break Galaxies (LBG)  

in the context of the Wide-Field Spectroscopic Telescope (WST).

The scientific framework and reference forecasts are based on:

> *The Wide-field Spectroscopic Telescope (WST) Science White Paper*  
> arXiv:2403.05398 (2024)

The repository covers:
- Photometric catalog preparation  
- Target selection implementation  
- Redshift distribution modeling  
- Forecast inputs (bias, clustering, efficiency)  
- Survey design tools  

---

# Repository Structure

---

## photometric_catalogs/

Preparation and manipulation of photometric datasets used for target selection.

- `compute_COSMOS_XMM_surface.py` — Surface density calculations.
- `degrade_photometry.py` — Simulate photometric depth degradation.
- `example_degrading_CLAUDS-r_band.ipynb` — Example degradation workflow.
- `explore_COSMOS.ipynb`, `explore_COSMOS_2.ipynb` — Catalog exploration.
- `mask_catalogs.ipynb` — Mask construction and application.
- `match_CLAUDSHSC_to_allWISE.ipynb` — Cross-matching workflow.
- `matching_catalog.py` — Catalog matching utilities.

---

## target_selection/

Core implementation of target selection algorithms.

### General utilities
- `catalog_infos.py` — Catalog metadata and configuration.
- `conversion.py` — Photometric conversions.
- `redshift_distribution_format.py` — Formatting tools for redshift distributions.

### bg_elg_lrg/
Placeholder directory for BGS, ELG, and LRG selection implementations.

### lyman_break_galaxies/
Placeholder directory for LBG selection implementations.

### quasars/
Placeholder for QSO target selection implementation.

### MagMax/
Placeholder directory for MagMax selection implementations.

---

### photom_redshift_distribution/
Placeholder directory for magnitude-photometric redshift distributions

## survey_design/

Survey strategy and observing configuration tools.

---

## forecasts/

Tools and data products used for cosmological forecasting and tracer validation.

---

## (OLD) redshift_distribution_white_paper/

Reference redshift distributions used in WST Science White Paper forecasts.

- `nz_wst_bg_bright.txt`
- `nz_wst_bg_faint.txt`
- `nz_wst_elg.txt`
- `nz_wst_lrg.txt`
- `nz_wst_qso.txt`
- `nz_wst_lbg-*.txt` (dropout selections)
- `plot_tracer_density.ipynb` — Plot tracer densities and compare distributions.

These files define the fiducial \( n(z) \) used in cosmological forecasts.