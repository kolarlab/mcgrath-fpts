import json
import re
from pathlib import Path


BASE_DIR = Path("..")
FASTA_CLEAN = BASE_DIR / "FASTA_CLEAN"

with open(BASE_DIR / "fpf_info.json") as f:
    fpf_info = json.load(f)


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


fpf_info_full = {}

for pdb_id, proteins in fpf_info.items():
    fasta_path = FASTA_CLEAN / f"{pdb_id}.fasta"
    if not fasta_path.exists():
        print(f"fasta file [{fasta_path}] not found")
        fpf_info_full[pdb_id] = proteins
        continue

    fasta_records = parse_fasta(fasta_path)

    updated_proteins = []
    for protein in proteins:
        protein_name = protein["protein_name"]
        fragment_sequence = protein["sequence"]
        chain_name = protein["chain_name"]

        candidates = fasta_records.get(protein_name)
        if not candidates:
            print(f"protein [{protein_name}] not found in {fasta_path}")
            updated_proteins.append(protein)
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
            updated_proteins.append(protein)
            continue

        header, full_sequence = match
        start = full_sequence.find(fragment_sequence)

        if start == -1:
            print(f"fragment sequence for [{protein_name}] in {pdb_id} not found in full sequence")
            updated_proteins.append(protein)
            continue

        occurrences = full_sequence.count(fragment_sequence)
        if occurrences > 1:
            print(f"multiple matches for [{protein_name}] in {pdb_id}, using first occurrence")

        end = start + len(fragment_sequence)

        new_protein = dict(protein)
        new_protein["full_sequence"] = full_sequence
        new_protein["fpf_position"] = [start, end]
        updated_proteins.append(new_protein)

    fpf_info_full[pdb_id] = updated_proteins

with open(BASE_DIR / "fpf_info_full.json", "w") as out:
    json.dump(fpf_info_full, out, indent=4)