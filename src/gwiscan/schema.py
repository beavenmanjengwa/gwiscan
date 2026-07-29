#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# schema.py - Shared hit-table schema.                                                             #
#                                                                                                  #
# The hit schema is produced by both search stages (search-hmm and search-diamond) and consumed    #
# by merge and extract-domains. Defining it once here lets the two hit tables concatenate cleanly. #
# family always holds the curated Family from the family table (e.g. Legume, GNA, EUL); accession, #
# evalue, bitscore, start and end are method-neutral -- for an hmmscan row they come from the HMM   #
# domain hit, for a DIAMOND row from the BLAST alignment (accession is then "-"; DIAMOND has none). #
#                                                                                                  #
####################################################################################################
"""

HIT_HEADER = [
    "protein_id",
    "family",
    "accession",
    "evalue",
    "bitscore",
    "start",
    "end",
    "method",
]
