# Supplementary data for "Ribosomes are covered by a coat of flexible protein fragments"

Authors of the publication: Hugo McGrath, Rudolf Kvasňovský, Michal Kolář

Preprint: https://www.biorxiv.org/content/10.64898/2026.06.18.733196v1

All authors are affiliated with the Department of physical chemistry at the University of Chemistry and Technology, Prague, Czechia.

The repository contains input FASTA and MMCIF files. Furthermore, the repository includes output data and Python scripts necessary to reproduce all figures in the manuscript.

Project structure:
```
.
├── SCRIPTS                                 # jupyter notebooks to make plots from result data
├── FASTA                                   # raw .fasta files from PDB
├── FASTA_CLEAN                             # .fasta of only ribosomal proteins, one .fasta per ribosome ({pdb_id}.fasta)
├── fpf_info_full.json                      # example .json structure:
                                            # {
                                            #   "3J7P": [ # pdb_id
                                            #       {
                                            #           "chain_name": "A",
                                            #           "protein_name": "Ribosomal protein uL2",
                                            #           "fragment_id": 0, # either 0 or 1, depends if C or N terminus
                                            #           "sequence": "LRGTKTVQEKEN", # sequence of FPT
                                            #           "length": 12, # length of FPT
                                            #           "full_sequence": "MGRVIRGQRKGAGSVFRAHVKHRKGAARLRAVDFAERHGYIKGIVKDIIHDPGRGAPLAKVVFRDPYRFKKRTELFIAAEGIHTGQFVYCGKKAQLNIGNVLPVGTMPEGTIVCCLEEKPGDRGKLARASGNYATVISHNPETKKTRVKLPSGSKKVISSANRAVVGVVAGGGRIDKPILKAGRAYHKYKAKRNCWPRVRGVAMNPVEHPFGGGNHQHIGKPSTIRRDAPAGRKVGLIAARRTGRLRGTKTVQEKEN", # full protein sequence
                                            #           "fpf_position": [
                                            #               245,
                                            #               257
                                            #           ] # index range of FPT in full protein sequence, python indexing
                                            #       },
                                            #       ...
                                            #   ],
                                            #   ...
                                            # }
├── fpf_info.json                           # example .json structure:
                                            # {
                                            #   "3J7P": [ # pdb_id
                                            #       {
                                            #           "chain_name": "A",
                                            #           "protein_name": "Ribosomal protein uL2",
                                            #           "fragment_id": 0, # either 0 or 1, depends if C or N terminus
                                            #           "sequence": "LRGTKTVQEKEN", # sequence of FPT
                                            #           "length": 12 # length of FPT
                                            #       },
                                            #       ...
                                            #   ],
                                            #   ...
                                            # }
├── fpf_info_presequence_corrected.json     # same as fpf_info_full.json, except "sequence" has targeting sequence removed
├── MMCIF                                   # raw .mmcif from PDB
├── PLOTS
├── PREDICTED_PDBS                          # .pdb structures of proteins including FPTs predicted by ESMFold
├── RESULTS_AA_CONTENT                      # {pdb_id}/{chain_name}_{fragment_id}_{aa_code}_content.npy
├── RESULTS_DSSP                            # {pdb_id}/{chain_name}.csv, secondary structure prediction
├── RESULTS_MITOFATES                       # {pdb_id}.tsv, prediction of mitoribosomal protein targeting sequences
├── RESULTS_PLDDT                           # {pdb_id}/{chain_name}.npy, ESMFold pLDDT prediction
└── ribosome_info.json                      # example .json structure:
                                            # [
                                            #   {
                                            #       "organism": "S. scrofa",
                                            #       "pdb_id": "3J7P",
                                            #       "resolution": 3.5,
                                            #       "ribosome_class": "eukaryotes"
                                            #   },
                                            #   ...
                                            # ]
```