#!/usr/bin/env python3
"""
Bandshape Equalisation for pulsar FITS archives.

Author : Avinash Kumar Paladi <avinashkumarpaladi@gmail.com>

Usage:
    python3 bandshape_equalisation.py <FITS-files>... [--sbin SBIN] [--ebin EBIN] [--sphase SPHASE] [--ephase EPHASE] [--output OUTPUT]

Arguments:
    fits_files       : Input PSRCHIVE FITS archive(s) (required)
    --sbin SBIN      : Start phase bin of on-pulse window (optional, default: 0)
    --ebin EBIN      : End phase bin of on-pulse window (optional, default: nbin)
    --sphase SPHASE  : Start phase of on-pulse window (0.0 to 1.0, overrides --sbin)
    --ephase EPHASE  : End phase of on-pulse window (0.0 to 1.0, overrides --ebin)
    --output OUTPUT  : Output filename (optional, default: <input>.beq.fits)

Examples:
    python3 bandshape_equalisation.py obs.fits
    python3 bandshape_equalisation.py obs.fits --sbin 100 --ebin 200
    python3 bandshape_equalisation.py obs.fits --sbin 100 --ebin 200 --output J0437.Band3.beq.fits
"""

__author__  = "Avinash Kumar Paladi"
__email__   = "avinashkumarpaladi@gmail.com"

import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')   # non-interactive backend; change to TkAgg/Qt5Agg for pop-up windows
import matplotlib.pyplot as plt

try:
    import psrchive
except ImportError:
    sys.exit("ERROR: psrchive Python bindings not found. "
             "Install PSRCHIVE with Python bindings and try again.")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_gmrt_band(arch):
    """
    Determines the uGMRT band (2, 3, 4, or 5) for a PSRCHIVE archive
    observed with GMRT, based on the central frequency.

    uGMRT Band Ranges:
        Band 2: 120  - 250  MHz
        Band 3: 250  - 500  MHz
        Band 4: 550  - 850  MHz
        Band 5: 1000 - 1460 MHz

    Parameters
    ----------
    arch : psrchive.Archive
        Loaded PSRCHIVE archive object.

    Returns
    -------
    x : int or None
        Band number (2, 3, 4, or 5), or None if telescope is not GMRT
        or frequency does not match any known band.
    """
    telescope = arch.get_telescope().strip().upper()
    if telescope not in ("GMRT", "UGMRT"):
        print(f"Telescope is '{telescope}', not GMRT. Returning None.")
        return None

    centre_freq = arch.get_centre_frequency()

    GMRT_BANDS = {
        2: (120,  250),
        3: (250,  500),
        4: (550,  850),
        5: (1000, 1460),
    }

    x = None
    for band, (f_low, f_high) in GMRT_BANDS.items():
        if f_low <= centre_freq <= f_high:
            x = band
            break

    if x is None:
        print(f"Warning: Central frequency {centre_freq:.2f} MHz does not "
              f"fall within any known uGMRT band.")
    return x


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def beq(filename, sbin=None, ebin=None, sphase=None, ephase=None, outfile=None):
    """
    Perform bandshape equalisation on a PSRCHIVE FITS file.

    Parameters
    ----------
    filename : str
        Path to the input FITS archive.
    sbin : int, optional
        First phase bin of the on-pulse window.
    ebin : int, optional
        Last phase bin (exclusive) of the on-pulse window.
    sphase : float, optional
        First phase of the on-pulse window (0.0 to 1.0). Overrides sbin.
    ephase : float, optional
        Last phase of the on-pulse window (0.0 to 1.0). Overrides ebin.
    outfile : str, optional
        Output filename. Defaults to <filename>.beq.fits.

    Returns
    -------
    str or None
        Path to the output file, or None if filename is falsy.
    """
    if not filename:
        return None

    arch = psrchive.Archive_load(filename)
    arch.dedisperse()
    nbin = arch.get_nbin()

    if sphase is not None:
        sbin = int(sphase * nbin)
    elif sbin is None:
        sbin = 0

    if ephase is not None:
        ebin = int(ephase * nbin)
    elif ebin is None:
        ebin = nbin

    # Ensure bounds
    sbin = max(0, min(sbin, nbin))
    ebin = max(0, min(ebin, nbin))

    data = arch.get_data()          # shape: (nsub, npol, nchan, nbin)
    arch.dededisperse()

    # Full-band pulse profile (summed over sub-integrations, polarisations, channels)
    pulseprofile = np.sum(data[:, :, :, :], axis=(0, 1, 2))

    # Per-channel bandshape within the selected pulse-window bins
    prof = np.sum(data[:, :, :, sbin:ebin], axis=(0, 1, 3))

    # Compute equalisation weights
    wts = np.max(prof) / prof
    wts[0]  = 0                                   # zero-weight edge channels
    wts[-1] = 0
    wts[prof < 0.01 * np.max(prof)] = 0           # mask very faint channels

    # Apply weights channel-by-channel
    nchan = arch.get_nchan()
    for ichan in range(nchan):
        profile = arch.get_Integration(0).get_Profile(0, ichan)
        profile.get_amps()[:] *= wts[ichan]

    # Determine output filename
    if outfile is None:
        outfile = f"{filename.removesuffix('.fits')}.beq.fits"

    arch.unload(outfile)
    print(f"[beq] Saved equalised archive -> {outfile}")

    # ---- Diagnostic plots ------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].plot(range(len(pulseprofile)), pulseprofile)
    axes[0].plot(
        range(len(pulseprofile))[sbin:ebin],
        pulseprofile[sbin:ebin],
        label="Selected Bins",
        color="orange",
    )
    axes[0].set_title("Pulse Profile")
    axes[0].set_xlabel("Phase bin")
    axes[0].set_ylabel("Intensity")
    axes[0].legend()

    axes[1].plot(range(len(prof)), prof)
    axes[1].set_title("Bandshape")
    axes[1].set_xlabel("Channel")
    axes[1].set_ylabel("Intensity")

    axes[2].plot(range(len(wts)), wts)
    axes[2].set_title("Weights")
    axes[2].set_xlabel("Channel")
    axes[2].set_ylabel("Weight")

    plt.tight_layout()
    diag_plot = f"{outfile.removesuffix('.fits')}.diagnostic.png"
    plt.savefig(diag_plot, dpi=150)
    plt.close(fig)
    print(f"[beq] Diagnostic plot saved -> {diag_plot}")

    # ---- Before / after comparison plot ----------------------------------
    plot(filename, outfile)

    return outfile


def plot(b3_filename, b5_filename=None):
    """
    Side-by-side frequency-phase plot of two archives (before / after).

    Parameters
    ----------
    b3_filename : str
        Path to the original archive.
    b5_filename : str, optional
        Path to the equalised archive.
        If None, both panels show the same file.
    """
    if not b5_filename:
        b5_filename = b3_filename

    arch_b3 = psrchive.Archive_load(b3_filename)
    arch_b5 = psrchive.Archive_load(b5_filename)

    arch_b3.dedisperse()
    arch_b3.remove_baseline()
    data_b3 = arch_b3.get_data()

    arch_b5.dedisperse()
    arch_b5.remove_baseline()
    data_b5 = arch_b5.get_data()

    b3_band = get_gmrt_band(arch_b3)
    b5_band = get_gmrt_band(arch_b5)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    im1 = axes[0].imshow(data_b3[0, 0, :, :], aspect='auto', cmap='hot', origin='lower')
    axes[0].set_title(
        f"Original Band {b3_band} - "
        f"{arch_b3.get_centre_frequency():.0f} MHz "
        f"{arch_b3.get_bandwidth():.0f} MHz",
        fontsize=13,
    )
    axes[0].set_xlabel("Phase Bins")
    axes[0].set_ylabel("Frequency Channels")
    plt.colorbar(im1, ax=axes[0], label="Intensity")

    im2 = axes[1].imshow(data_b5[0, 0, :, :], aspect='auto', cmap='hot', origin='lower')
    axes[1].set_title(
        f"Bandequalised Band {b5_band} - "
        f"{arch_b5.get_centre_frequency():.0f} MHz "
        f"{arch_b5.get_bandwidth():.0f} MHz",
        fontsize=13,
    )
    axes[1].set_xlabel("Phase Bins")
    axes[1].set_ylabel("Frequency Channels")
    plt.colorbar(im2, ax=axes[1], label="Intensity")

    fig.suptitle(f"Pulsar: {arch_b3.get_source()}", fontsize=15, fontweight='bold')
    plt.tight_layout()

    comparison_plot = f"{b5_filename.removesuffix('.fits')}.comparison.png"
    plt.savefig(comparison_plot, dpi=150)
    plt.close(fig)
    print(f"[plot] Comparison plot saved -> {comparison_plot}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        prog="bandshape_equalisation.py",
        description=(
            "Bandshape equalisation for uGMRT pulsar FITS archives.\n"
            "Computes per-channel weights from the on-pulse bandshape and\n"
            "applies them to flatten the frequency response of the archive.\n\n"
            f"Author : {__author__} <{__email__}>"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 bandshape_equalisation.py obs.fits\n"
            "  python3 bandshape_equalisation.py *.fits --sbin 100 --ebin 200\n"
            "  python3 bandshape_equalisation.py obs.fits --sbin 100 --ebin 200 "
            "--output J0437.Band3.beq.fits"
        ),
    )

    parser.add_argument(
        "fits_files",
        nargs="+",
        metavar="FITS-files",
        help="Input PSRCHIVE FITS archive(s).",
    )
    parser.add_argument(
        "--sbin",
        type=int,
        default=None,
        metavar="SBIN",
        help="Start phase bin of the on-pulse window (default: 0).",
    )
    parser.add_argument(
        "--ebin",
        type=int,
        default=None,
        metavar="EBIN",
        help="End phase bin of the on-pulse window (default: nbin).",
    )
    parser.add_argument(
        "--sphase",
        type=float,
        default=None,
        metavar="SPHASE",
        help="Start phase of the on-pulse window (0.0 to 1.0). Overrides --sbin.",
    )
    parser.add_argument(
        "--ephase",
        type=float,
        default=None,
        metavar="EPHASE",
        help="End phase of the on-pulse window (0.0 to 1.0). Overrides --ebin.",
    )
    parser.add_argument(
        "--output",
        default=None,
        metavar="OUTPUT",
        help="Output FITS filename (default: <input>.beq.fits). Ignored if multiple files are provided.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Bandshape Equalisation  |  {__author__} <{__email__}>")
    print(f"  Inputs : {len(args.fits_files)} file(s)")
    if args.sphase is not None:
        print(f"  sphase : {args.sphase}")
    else:
        print(f"  sbin   : {args.sbin}")
    if args.ephase is not None:
        print(f"  ephase : {args.ephase}")
    else:
        print(f"  ebin   : {args.ebin}")
    if len(args.fits_files) > 1 and args.output:
        print("  Warning: --output is ignored when multiple input files are provided.")
    else:
        print(f"  Output : {args.output or '<auto>'}")
    print()

    for f in args.fits_files:
        print(f"Processing: {f}")
        outfile = args.output if len(args.fits_files) == 1 else None
        result = beq(f, sbin=args.sbin, ebin=args.ebin, sphase=args.sphase, ephase=args.ephase, outfile=outfile)

        if result:
            print(f"Done. Output archive: {result}\n")
        else:
            print(f"No output produced for {f}\n")


if __name__ == "__main__":
    main()
