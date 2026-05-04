"""
Utility module to compute and save 2D redshift–magnitude histograms.

This module defines:
- Default binning for redshift (z) and magnitude (mag)
- Midpoint arrays for bins
- A helper function `save_targets` to compute a 2D histogram
  and store the result as a compressed NumPy (.npz) file
"""

import os
import numpy as np


# -----------------------------------------------------------------------------
# Default bin definitions
# -----------------------------------------------------------------------------
Z_EDGES = np.linspace(0, 6, 60)
MAG_EDGES = np.linspace(17, 28, 50)

Z_MID = 0.5 * (Z_EDGES[:-1] + Z_EDGES[1:])
MAG_MID = 0.5 * (MAG_EDGES[:-1] + MAG_EDGES[1:])

BINS = [Z_EDGES, MAG_EDGES]


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def save_targets(
    z_sample,
    mag_sample,
    surface,
    where_to_save="../photom_redshift_distribution/",
    field=None,
    tracer=None,
    selection=None,
    add_name_selection="",
    info_about_the_sample="",
):
    """
    Compute a 2D histogram in (z, magnitude) and save it to disk.

    Parameters
    ----------
    z_sample : array-like
        Redshift values.
    mag_sample : array-like
        Magnitude values.
    surface : float
        Survey surface in deg^2.
    where_to_save : str, optional
        Directory where the output file will be stored.
    field : str
        Field name (used in output filename).
    tracer : str
        Tracer name (used in output filename).
    selection : str
        Selection name (used in output filename).
    add_name_selection : str, optional
        Additional suffix appended to the filename.
    info_about_the_sample : str, optional
        Free-text metadata stored in the output file.

    Returns
    -------
    None
        Saves a compressed .npz file containing histogram and metadata.
    """

    if field is None or tracer is None or selection is None:
        raise ValueError("'field', 'tracer', and 'selection' must be provided.")

    # Ensure output directory exists
    os.makedirs(where_to_save, exist_ok=True)

    # Compute 2D histogram
    counts, z_edges, mag_edges = np.histogram2d(
        z_sample,
        mag_sample,
        bins=BINS,
    )

    # Build filename
    name_parts = [field, tracer, selection]
    if add_name_selection:
        name_parts.append(add_name_selection)

    filename = "_".join(name_parts)
    filepath = os.path.join(where_to_save, f"{filename}.npz")

    # Save data
    np.savez(
        filepath,
        z_center=Z_MID,
        mag_center=MAG_MID,
        z_edges=z_edges,
        mag_edges=mag_edges,
        object_count=counts,
        surface_deg2=surface,
        info_about_the_sample=info_about_the_sample,
    )

    print(f"Saved in: {where_to_save}")
    print(f"Filename: {filename}.npz")

    return None

             

