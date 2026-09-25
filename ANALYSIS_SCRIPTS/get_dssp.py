import os
import csv
import pydssp
from Bio.PDB import PDBParser


BASE_DIR = ".."
PREDICTED_PDBS_DIR = f"{BASE_DIR}/PREDICTED_PDBS"
RESULTS_DSSP_DIR = f"{BASE_DIR}/RESULTS_DSSP"

DSSP_COLUMNS = ["residue_number", "residue_name", "secondary_structure"]


def process_chain(pdb_id, chain_id, pdb_path, output_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(chain_id, pdb_path)
    model = structure[0]
    residues = [r for r in model.get_residues() if r.id[0] == " "]

    with open(pdb_path, "r") as f:
        pdb_text = f.read()

    coords = pydssp.read_pdbtext(pdb_text)
    dssp_codes = pydssp.assign(coords)

    rows = [
        {
            "residue_number": res.id[1],
            "residue_name": res.resname,
            "secondary_structure": code,
        }
        for res, code in zip(residues, dssp_codes)
    ]

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=DSSP_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main():
    os.makedirs(RESULTS_DSSP_DIR, exist_ok=True)

    for pdb_id in os.listdir(PREDICTED_PDBS_DIR):
        pdb_subdir = os.path.join(PREDICTED_PDBS_DIR, pdb_id)
        if not os.path.isdir(pdb_subdir):
            continue
        for fname in os.listdir(pdb_subdir):
            if not fname.endswith(".pdb"):
                continue
            chain_id = fname[:-4]
            pdb_path = os.path.join(pdb_subdir, fname)
            output_path = os.path.join(RESULTS_DSSP_DIR, pdb_id, f"{chain_id}.csv")
            if not os.path.exists(output_path):
                process_chain(pdb_id, chain_id, pdb_path, output_path)


if __name__ == "__main__":
    main()