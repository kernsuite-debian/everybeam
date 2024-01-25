#!/usr/bin/env python3
# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
    Script converts OVRO-LWA simulation data to a coefficient file, by fitting spherical harmonics to
    the simulated data. In order to run this script, please make sure
    that the everybeam shared library is on your LD_LIBRARY_PATH and on your PYTHONPATH, See below for
    example export commands:
    - export LD_LIBRARY_PATH=~/opt/everybeam/lib:$LD_LIBRARY_PATH
    - export PYTHONPATH=/home~/opt/everybeam/lib/python3.6/site-packages:$PYTHONPATH

    You can run this script as follows:
    python3 scripts/coeff_scripts/convert_lwa.py  --simulation_path LWA_SIMULATION.h5 --destination_path LWA_COEFFICIENTS.h5
"""
import argparse
import numpy as np
import h5py
from convert_lobes import fit_and_write_lobes_coeffs
import os


def read_lwa_simulation_file(path):
    """
    Read simulation files, where all frequencies are contained in a single HDF5 file

    Parameters
    ----------
    path : str
        Path to directory containing simulated results

    Returns
    -------
    dict
        Dictionary, containing the following fields:
        - theta: zenith angles of simulated results [rad]
        - phi: elevation angle of simulated results [rad]
        - freq, np.1darray of frequencies [MHz]
        - v_pol1, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
        - v_pol2, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
    """

    h5 = h5py.File(path, "r")

    freqs = np.array(h5["Freq(Hz)"])

    # Read polarizations from FEKO simulation file
    # For compatibility with the LOBES coefficients, add one extra dimension for the stations.
    # Since the data is the same for all stations, the dimension of this new axis is 1.
    v_pol1 = np.expand_dims(
        np.stack(
            (
                np.array(h5["X-pol_Efields/etheta"]),
                np.array(h5["X-pol_Efields/ephi"]),
            )
        ),
        axis=0,
    )
    v_pol2 = np.expand_dims(
        np.stack(
            (
                np.array(h5["Y-pol_Efields/etheta"]),
                np.array(h5["Y-pol_Efields/ephi"]),
            )
        ),
        axis=0,
    )

    # reshape from shape (#antenna elements, 2, #freqs, #theta, #phi)
    #           to shape (#antenna elements, 2, #theta, #phi, #freqs)
    v_pol1 = np.transpose(v_pol1, (0, 1, 3, 4, 2))
    v_pol2 = np.transpose(v_pol2, (0, 1, 3, 4, 2))

    # Theta and phi are already in radians
    theta = np.array(h5["theta_pts"])
    phi = np.array(h5["phi_pts"])

    if np.any(freqs[:-1] > freqs[1:]):
        # Then we have to reorder frequencies - and simulated data -
        # ascendingly
        idcs = np.argsort(freqs)
        freqs = freqs[idcs]
        v_pol1 = v_pol1[..., idcs]
        v_pol2 = v_pol2[..., idcs]

    theta = theta.flatten()
    phi = phi.flatten()

    # Return as dictionary
    # Note: the phi axis runs from [0, 2pi> (not including 2 pi) by cutting off
    # the last element. This avoids 0 and 2pi to be weighted twice.
    return {
        "theta": theta,
        "phi": phi[:-1],
        "v_pol1": v_pol1[..., :-1, :],
        "v_pol2": v_pol2[..., :-1, :],
        "freq": freqs,
    }


def main(simulation_path, destination_path):
    """
    This script fits spherical harmonics coefficients to a FEKO simulation for the OVRO-LWA telescope.
    The simulation file is provided as an HDF5 file.
    """

    # Order of Legendre polynomial to be fitted. Adjust to suit your needs
    nmax = 21

    """
    Simulated results for each frequency stored in separate file. Station coordinates
    stored in separate HDF5 file.
    """

    # Read simulated data from file
    simdata = read_lwa_simulation_file(simulation_path)

    fit_and_write_lobes_coeffs(
        simdata["theta"],
        simdata["phi"],
        simdata["freq"],
        None,
        simdata["v_pol1"],
        simdata["v_pol2"],
        nmax,
        destination_path,
        apply_weighting=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to generate coeffiecients file for the OVRO-LWA telescope, based on a simulation file"
    )
    parser.add_argument(
        "--simulation_path", help="Path to LWA simulation file."
    )
    parser.add_argument(
        "--destination_path", help="Path to the output coefficients file."
    )
    args = parser.parse_args()
    main(args.simulation_path, args.destination_path)
