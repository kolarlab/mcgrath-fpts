import json
import numpy as np
from pathlib import Path
import shutil

BASE_DIR = ".."
RESULTS_AA_CONTENT = f"{BASE_DIR}/RESULTS_AA_CONTENT"

MIN_LENGTH = 20

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def compute_metrics(sequence):
    n = len(sequence)
    metrics = {aa.lower() + "_content": sequence.count(aa) / n for aa in AMINO_ACIDS}
    metrics["aromatic_content"] = (sequence.count("Y") + sequence.count("F") + sequence.count("W") + sequence.count("H")) / n
    metrics["sequence_length"] = float(n)
    return metrics


def main():
    results_dir = Path(RESULTS_AA_CONTENT)

    if results_dir.exists():
        shutil.rmtree(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    with open(Path(BASE_DIR) / "fpf_info.json") as f:
        fpf_info = json.load(f)

    for pdb_id, fragments in sorted(fpf_info.items()):
        for frag in fragments:
            sequence = frag["sequence"]
            if len(sequence) < MIN_LENGTH:
                continue

            chain_id = frag["chain_name"]
            fragment_id = frag["fragment_id"]

            metrics = compute_metrics(sequence)

            out_dir = results_dir / pdb_id
            out_dir.mkdir(parents=True, exist_ok=True)

            for analysis_name, value in metrics.items():
                out_path = out_dir / f"{chain_id}_{fragment_id}_{analysis_name}.npy"
                np.save(out_path, np.float64(value))


if __name__ == "__main__":
    main()