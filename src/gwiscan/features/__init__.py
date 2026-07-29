"""Per-candidate annotation stages: physicochemical properties, signal/transit
peptides (TargetP), transmembrane topology (DeepTMHMM), subcellular localization
(DeepLoc), and domain/GO annotation (InterProScan). Each only depends on
intermediate/candidates.fasta, so they form independent parallel branches.
"""
