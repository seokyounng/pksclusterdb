#!/usr/bin/env python3
import csv
import sys
import json
import argparse
from pathlib import Path
from Bio import SeqIO
import re

csv.field_size_limit(sys.maxsize)

# ---------------------------
# Helper functions
# ---------------------------
def parse_coords(coord_str):
    try:
        start, end = map(int, coord_str.split('-'))
        return start, end
    except Exception:
        return None, None

def overlap(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or b_end < a_start)

def parse_mibig_entry(gbk_file: Path, json_file: Path):
    mibig_id = gbk_file.stem.upper()
    accession, product, organism = "", "", ""
    start, end = None, None

    try:
        with open(json_file) as jf:
            js = json.load(jf)
            if "loci" in js and js["loci"]:
                accession = js["loci"][0].get("accession", "").split(".")[0]
            if "compounds" in js:
                product = "; ".join([c.get("name", "") for c in js["compounds"] if "name" in c])
            if "taxonomy" in js and "name" in js["taxonomy"]:
                organism = js["taxonomy"]["name"].strip()
    except Exception as e:
        print(f"[WARN] Could not parse JSON {json_file}: {e}")

    try:
        gbk_text = gbk_file.read_text()
        m_header = re.search(r"ORGANISM\s+.*?\n(.*?)\nFEATURES", gbk_text, re.DOTALL)
        header_section = m_header.group(1) if m_header else ""
        m_start = re.search(r"Orig\. start\s*::\s*(\d+)", header_section)
        m_end   = re.search(r"Orig\. end\s*::\s*(\d+)", header_section)
        if m_start and m_end:
            start, end = int(m_start.group(1)), int(m_end.group(1))

        if start is None or end is None:
            record = SeqIO.read(str(gbk_file), "genbank")
            for feature in record.features:
                if feature.type == "region":
                    start, end = int(feature.location.start)+1, int(feature.location.end)
                    break
            if start is None or end is None:
                for feature in record.features:
                    if feature.type == "misc_feature":
                        start, end = int(feature.location.start)+1, int(feature.location.end)
                        break
            if start is None or end is None:
                for feature in record.features:
                    if feature.type == "source":
                        start, end = int(feature.location.start)+1, int(feature.location.end)
                        break
    except Exception as e:
        print(f"[WARN] Could not parse coordinates from {gbk_file}: {e}")

    if not accession or accession.upper() == "MIBIG":
        return None

    return {"accession": accession, "start": start, "end": end,
            "mibig_id": mibig_id, "product": product, "organism": organism}

def parse_mibig_dir(mibig_gbk_dir: Path, mibig_json_dir: Path):
    entries = []
    for gbk_file in mibig_gbk_dir.glob("*.gbk"):
        json_file = mibig_json_dir / (gbk_file.stem + ".json")
        if not json_file.exists():
            print(f"[WARN] Missing JSON for {gbk_file}, skipping")
            continue
        entry = parse_mibig_entry(gbk_file, json_file)
        if entry:
            entries.append(entry)
    return entries

def load_or_parse_mibig(mibig_gbk_dir: Path, mibig_json_dir: Path, mibig_table: Path):
    entries = []
    if mibig_table.exists():
        with open(mibig_table) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                row["start"] = int(row["start"])
                row["end"] = int(row["end"])
                entries.append(row)
    else:
        entries = parse_mibig_dir(mibig_gbk_dir, mibig_json_dir)
        entries.sort(key=lambda e: e["mibig_id"])
        with open(mibig_table, "w", newline="") as f:
            writer = csv.DictWriter(f, delimiter="\t",
                                    fieldnames=["accession", "start", "end",
                                                "mibig_id", "product", "organism"])
            writer.writeheader()
            for e in entries:
                writer.writerow(e)

    mibig_index = {}
    for e in entries:
        mibig_index.setdefault(e["accession"], []).append(e)
    return mibig_index

# ---------------------------
# Main function
# ---------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Annotate a final PKS cluster summary with overlapping MiBIG entries."
    )
    parser.add_argument("--summary", required=True, help="Input final_summary.tsv from cleanup_clusters.py")
    parser.add_argument("--mibig_gbk_dir", default="data/mibig/mibig_gbk_4.0")
    parser.add_argument("--mibig_json_dir", default="data/mibig/mibig_json_4.0")
    parser.add_argument("--mibig_table", default="outputs/mibig_clusters.tsv",
                        help="Cached parsed MiBIG table to read/write")
    parser.add_argument("--output_full", default="outputs/final_summary_with_mibig.tsv")
    parser.add_argument("--output_known", default="outputs/final_summary_known.tsv")
    args = parser.parse_args()

    summary_file = args.summary
    mibig_gbk_dir = Path(args.mibig_gbk_dir)
    mibig_json_dir = Path(args.mibig_json_dir)
    mibig_table = Path(args.mibig_table)
    output_full = args.output_full
    output_known = args.output_known

    mibig_table.parent.mkdir(parents=True, exist_ok=True)
    Path(output_full).parent.mkdir(parents=True, exist_ok=True)
    Path(output_known).parent.mkdir(parents=True, exist_ok=True)

    mibig_entries = load_or_parse_mibig(mibig_gbk_dir, mibig_json_dir, mibig_table)
    print(f"Loaded {sum(len(v) for v in mibig_entries.values())} MiBIG entries.")

    with open(summary_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        summary_rows = list(reader)

    fieldnames = list(summary_rows[0].keys()) + ["known product", "mibig ID/accession"]
    with open(output_full, "w", newline="") as fout_full, open(output_known, "w", newline="") as fout_known:
        writer_full = csv.DictWriter(fout_full, delimiter="\t", fieldnames=fieldnames)
        writer_known = csv.DictWriter(fout_known, delimiter="\t", fieldnames=fieldnames)
        writer_full.writeheader()
        writer_known.writeheader()

        for row in summary_rows:
            acc = row["accession"].split(".")[0]
            start, end = parse_coords(row["coordinates"])
            row["known product"] = ""
            row["mibig ID/accession"] = ""
            found = False

            # 1) Representative -- coordinate-based
            if acc in mibig_entries and start and end:
                for entry in mibig_entries[acc]:
                    if overlap(start, end, entry["start"], entry["end"]):
                        row["known product"] = entry["product"]
                        row["mibig ID/accession"] = f"{entry['mibig_id']}/{entry['accession']}"
                        found = True
                        break

            # 2) Duplicates / redundant (NO coordinate lookup, match only by accession)
            if not found:
                for col in ["identical_members", "homolog_members"]:
                    members = [m.strip() for m in row[col].split(",")] if row[col] else []
                    for m in members:
                        if "_c" not in m:
                            continue
                        dup_acc = m.split("_c")[0]

                        if dup_acc in mibig_entries:
                            entry = mibig_entries[dup_acc][0]  # take first match
                            row["known product"] = entry["product"]
                            row["mibig ID/accession"] = (
                                f"{entry['mibig_id']}/{entry['accession']} | "
                                f"duplicate: {m} | representative: {acc}"
                            )
                            print(f"[INFO] MiBIG match by duplicate/redundant member: {m} "
                                  f"→ {entry['mibig_id']} ({entry['accession']})")
                            found = True
                            break
                    if found:
                        break

            writer_full.writerow(row)
            if row["known product"]:
                writer_known.writerow(row)

    print(f"Done. Full table: {output_full}")
    print(f"Known clusters table: {output_known}")

if __name__ == "__main__":
    main()
