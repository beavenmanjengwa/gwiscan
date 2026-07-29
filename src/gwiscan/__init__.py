#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# __init__.py - GWIscan package: genome-wide gene-family annotation pipeline.                      #
#                                                                                                  #
# Exposes one console command (`gwiscan`) whose subcommands are the individual pipeline stages     #
# (see gwiscan.cli). Each stage is a small module with a run(cfg) function operating on a shared   #
# gwiscan.config.Config, so the stages share path handling, the hit schema, and I/O helpers.       #
#                                                                                                  #
####################################################################################################
"""

__version__ = "1.0.0"
