#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# __init__.py - GWIscan: genome-wide identification and annotation of gene families.               #
#                                                                                                  #
# Exposes one console command (`gwiscan`) whose subcommands are the individual pipeline stages     #
# (see gwiscan.cli). Each stage is a small module with a run(cfg) function operating on a shared   #
# gwiscan.config.Config, so the stages share path handling, the hit schema, and I/O helpers.       #
#                                                                                                  #
####################################################################################################
"""

from importlib.metadata import PackageNotFoundError, version

try:
    # Single source of truth: the version declared in pyproject.toml, read from the
    # installed package metadata. Never hardcode it here -- that drifts on releases.
    __version__ = version("gwiscan")
except PackageNotFoundError:          # imported from a source tree without an install
    __version__ = "0.0.0+unknown"
