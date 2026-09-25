import json
import subprocess
from pathlib import Path


BASE_DIR = Path("..")
FASTA_CLEAN = BASE_DIR / "FASTA_CLEAN"
RESULTS_MITOFATES = BASE_DIR / "RESULTS_MITOFATES"

MITOFATES_PL = Path().home() / "MitoFates/script/MitoFates.pl"

ORGANISM_MAP = {
    "A. thaliana/B. oleracea": "plant",
    "H. sapiens": "metazoa",
    "N. crassa": "fungi",
    "S. scrofa": "metazoa",
}


def classify_organism(organism):
    if organism in ORGANISM_MAP:
        return ORGANISM_MAP[organism]
    return None


with open(BASE_DIR / "ribosome_info.json") as f:
    ribosome_info = json.load(f)

mito_entries = [e for e in ribosome_info if e["ribosome_class"] == "mitochondria"]

pdb_to_flag = {}

for entry in mito_entries:
    flag = classify_organism(entry["organism"])
    pdb_to_flag[entry["pdb_id"]] = flag

RESULTS_MITOFATES.mkdir(exist_ok=True)

for pdb_id, organism_flag in pdb_to_flag.items():
    fasta_path = FASTA_CLEAN / f"{pdb_id}.fasta"
    if not fasta_path.exists():
        print(f"fasta file [{fasta_path}] not found")
        continue

    out_path = RESULTS_MITOFATES / f"{pdb_id}.tsv"
    print(f"Running MitoFates on {pdb_id} ({organism_flag})")

    with open(out_path, "w") as out:
        result = subprocess.run(
            ["perl", str(MITOFATES_PL), str(fasta_path), organism_flag],
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
        )

    if result.returncode != 0:
        print(f"  MitoFates failed for {pdb_id}: {result.stderr.strip()}")