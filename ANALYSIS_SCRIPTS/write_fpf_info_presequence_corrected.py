import csv
import json
import re
from pathlib import Path


BASE_DIR = Path("..")
FASTA_CLEAN = BASE_DIR / "FASTA_CLEAN"
RESULTS_MITOFATES = BASE_DIR / "RESULTS_MITOFATES"

PRESEQUENCE_PROBABILITY_THRESHOLD = 0.385

with open(BASE_DIR / "fpf_info_full.json") as f:
    fpf_info_full = json.load(f)

with open(BASE_DIR / "ribosome_info.json") as r:
    ribosome_info = json.load(r)

ribosome_class = {entry["pdb_id"]: entry["ribosome_class"] for entry in ribosome_info}


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


def load_mitofates_results(path):
    results = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq_id = row["Sequence ID"]
            probability = float(row["Probability of presequence"])
            cleavage_field = row["Cleavage site (processing enzyme)"]
            match = re.search(r"(\d+)\(MPP\)", cleavage_field)
            mpp_position = int(match.group(1)) if match else None
            results[seq_id] = {
                "probability": probability,
                "mpp_position": mpp_position,
            }
    return results


fpf_info_corrected = {}

for pdb_id, proteins in fpf_info_full.items():
    if ribosome_class.get(pdb_id) != "mitochondria":
        fpf_info_corrected[pdb_id] = proteins
        continue

    fasta_path = FASTA_CLEAN / f"{pdb_id}.fasta"
    if not fasta_path.exists():
        print(f"fasta file [{fasta_path}] not found, leaving {pdb_id} unmodified")
        fpf_info_corrected[pdb_id] = proteins
        continue

    mitofates_path = RESULTS_MITOFATES / f"{pdb_id}.tsv"
    if not mitofates_path.exists():
        print(f"mitofates results [{mitofates_path}] not found, leaving {pdb_id} unmodified")
        fpf_info_corrected[pdb_id] = proteins
        continue

    fasta_records = parse_fasta(fasta_path)
    mitofates_results = load_mitofates_results(mitofates_path)
    header_cache = {}

    updated_proteins = []
    for protein in proteins:
        protein_name = protein["protein_name"]
        chain_name = protein["chain_name"]

        cache_key = (protein_name, chain_name)
        if cache_key not in header_cache:
            candidates = fasta_records.get(protein_name, [])
            match_header = None
            for header, _ in candidates:
                chains = header_chains(header)
                if chains is not None and chain_name in chains:
                    match_header = header
                    break
            header_cache[cache_key] = match_header

        header = header_cache[cache_key]
        if header is None:
            print(
                f"chain [{chain_name}] for protein [{protein_name}] in {pdb_id} "
                f"not found in {fasta_path}, leaving entry unmodified"
            )
            updated_proteins.append(protein)
            continue

        seq_id = header[1:]
        if seq_id not in mitofates_results:
            print(
                f"no mitofates result for [{seq_id}] in {pdb_id}, leaving entry unmodified"
            )
            updated_proteins.append(protein)
            continue

        result = mitofates_results[seq_id]
        probability = result["probability"]

        if probability < PRESEQUENCE_PROBABILITY_THRESHOLD:
            new_protein = dict(protein)
            new_protein["presequence_probability"] = probability
            updated_proteins.append(new_protein)
            continue

        if result["mpp_position"] is None:
            print(
                f"probability above threshold but no MPP position parsed for [{seq_id}] "
                f"in {pdb_id}, leaving entry unmodified apart from probability"
            )
            new_protein = dict(protein)
            new_protein["presequence_probability"] = probability
            updated_proteins.append(new_protein)
            continue

        mpp_position = result["mpp_position"]
        presequence_range = [0, mpp_position]

        start, end = protein["fpf_position"]

        if mpp_position <= start:
            new_protein = dict(protein)
            new_protein["presequence_probability"] = probability
            new_protein["presequence_range"] = presequence_range
            updated_proteins.append(new_protein)
            continue

        if mpp_position >= end:
            print(
                f"fragment [{chain_name}/{protein_name}, fragment_id={protein['fragment_id']}] "
                f"in {pdb_id} fully overlaps predicted presequence "
                f"(presequence 0-{mpp_position}, fragment {start}-{end}), dropping entry"
            )
            continue

        removed = mpp_position - start
        new_sequence = protein["sequence"][removed:]
        new_protein = dict(protein)
        new_protein["sequence"] = new_sequence
        new_protein["length"] = len(new_sequence)
        new_protein["fpf_position"] = [mpp_position, end]
        new_protein["presequence_probability"] = probability
        new_protein["presequence_range"] = presequence_range
        updated_proteins.append(new_protein)
        print(
            f"trimmed [{chain_name}/{protein_name}, fragment_id={protein['fragment_id']}] "
            f"in {pdb_id}: removed {removed} residue(s) overlapping predicted presequence"
        )

    fpf_info_corrected[pdb_id] = updated_proteins

with open(BASE_DIR / "fpf_info_presequence_corrected.json", "w") as out:
    json.dump(fpf_info_corrected, out, indent=4)