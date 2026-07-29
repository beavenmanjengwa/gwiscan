#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# __main__.py - Enable `python -m gwiscan` as an alias for the `gwiscan` command.                  #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
