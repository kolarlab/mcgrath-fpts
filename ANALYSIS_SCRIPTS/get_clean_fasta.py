import json
import re
from pathlib import Path


BASE_DIR = Path("..")
FASTA = BASE_DIR / "FASTA"
FASTA_CLEAN = BASE_DIR / "FASTA_CLEAN"

if FASTA_CLEAN.exists():
    for f in FASTA_CLEAN.glob("*.fasta"):
        f.unlink()
else:
    FASTA_CLEAN.mkdir()

with open(BASE_DIR / "fpf_info.json") as f:
    fpf_info = json.load(f)

with open(BASE_DIR / "ribosome_info.json") as r:
    ribosome_info = json.load(r)


def parse_fasta(path):
    records = {}
    header = None
    seq_lines = []

    def store(header, sequence):
        fields = header[1:].split("|")
        if len(fields) >= 3:
            records.setdefault(fields[2], []).append((header, sequence))

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    store(header, "".join(seq_lines))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            store(header, "".join(seq_lines))
    return records


def header_chains(header):
    fields = header[1:].split("|")
    if len(fields) < 2:
        return None
    field = fields[1]
    matches = re.findall(r"\[auth\s+([^\]]+)\]", field)
    if matches:
        last = matches[-1].split(",")[0].strip()
        return [last]
    chain_field = re.sub(r"^Chains?\s+", "", field)
    return [c.strip() for c in chain_field.split(",")]


for entry in ribosome_info:
    pdb_id = entry["pdb_id"]
    if pdb_id not in fpf_info:
        continue

    fasta_path = FASTA / f"{pdb_id}.fasta"
    if not fasta_path.exists():
        print(f"fasta file [{fasta_path}] not found")
        continue

    fasta_records = parse_fasta(fasta_path)
    written_headers = set()

    for protein in fpf_info[pdb_id]:
        if protein["fragment_id"] != 0:
            continue

        protein_name = protein["protein_name"]
        chain_name = protein["chain_name"]

        candidates = fasta_records.get(protein_name)
        if not candidates:
            print(f"protein [{protein_name}] not found")
            continue

        match = None
        for header, sequence in candidates:
            chains = header_chains(header)
            if chains is not None and chain_name in chains:
                match = (header, sequence)
                break

        if match is None:
            print(
                f"chain [{chain_name}] for protein [{protein_name}] in {pdb_id} "
                f"not found among {len(candidates)} candidate record(s)"
            )
            continue

        header, sequence = match
        if header in written_headers:
            continue
        written_headers.add(header)

        out_path = FASTA_CLEAN / f"{pdb_id}.fasta"
        with open(out_path, "a") as out:
            out.write(f"{header}\n{sequence}\n")