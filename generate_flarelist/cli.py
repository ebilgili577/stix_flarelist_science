#!/usr/bin/env python3
"""
STIX Flare List Fetching CLI.

Entry point for the stix-fetch command.
Fetches and processes flare data from STIX.
"""

from astropy.time import Time

from .flarelist_generate import get_flares

# Default configuration
RAW_FITS_PATH = 'data/raw_fits'


def main():
    """Main entry point for stix-fetch command."""
    # TODO: Add argparse for configurable date ranges and paths
    tstart = Time('2025-09-01', format="iso")
    tend = Time('2026-01-22', format="iso")
    
    flares = get_flares(tstart, tend, RAW_FITS_PATH)
    print(f"Fetched and processed {len(flares)} flares")


if __name__ == "__main__":
    main()
