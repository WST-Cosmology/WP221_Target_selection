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

---

# Repository Structure

---

## photometric_catalogs/

Preparation and manipulation of photometric datasets used for target selection.

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

Reference redshift distributions used in WST Science White Paper forecasts, and summarized in ["WST 221 Technical note"](https://docs.google.com/document/d/1RkDFxtx-Jb-XMQDaSXwRJVphXaIqcfKyYUigVJ7xyJw/edit?usp=sharing). 

- `nz_wst_bg_bright.txt`
- `nz_wst_bg_faint.txt`
- `nz_wst_elg.txt`
- `nz_wst_lrg.txt`
- `nz_wst_qso.txt`
- `nz_wst_lbg-*.txt` (dropout selections)
- `plot_tracer_density.ipynb` — Plot tracer densities and compare distributions.

These files define the fiducial \( n(z) \) used in cosmological forecasts.