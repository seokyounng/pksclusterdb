#!/usr/bin/env python3
import argparse
import datetime
import json
import re
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO


# ---------------------------
# Parsing helpers
# ---------------------------

def build_accession_desc_dict(desc_file):
    acc_desc = {}
    with open(desc_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                acc_desc[parts[0]] = parts[1]
    return acc_desc


def parse_dates(date_file):
    MONTHS = {
        "JAN": "01", "FEB": "02", "MAR": "03", "APR": "04",
        "MAY": "05", "JUN": "06", "JUL": "07", "AUG": "08",
        "SEP": "09", "OCT": "10", "NOV": "11", "DEC": "12",
    }
    accession_to_date = {}
    with open(date_file) as f:
        for line in f:
            if "\t" not in line:
                continue
            acc, date_str = line.strip().split("\t")
            try:
                day, mon, year = date_str.strip().split("-")
                iso_date = f"{year}-{MONTHS[mon.upper()]}-{day.zfill(2)}"
                accession_to_date[acc] = datetime.datetime.strptime(iso_date, "%Y-%m-%d")
            except Exception:
                accession_to_date[acc] = datetime.datetime.max
    return accession_to_date


# ---------------------------
# Antismash parsing
# ---------------------------

def count_domains_in_gbk(gbk_file: Path, domain_pattern: str) -> int:
    count = 0
    try:
        for record in SeqIO.parse(str(gbk_file), "genbank"):
            for feature in record.features:
                if feature.type == "aSDomain":
                    if "aSDomain" in feature.qualifiers:
                        if any(domain_pattern in v for v in feature.qualifiers["aSDomain"]):
                            count += 1
    except Exception:
        pass
    return count


def get_cluster_coords(js_file: Path, cluster_idx: int, seq_id: str):
    if not js_file.exists():
        return None
    text = js_file.read_text()
    m = re.search(r"var\s+(all_regions|recordData)\s*=\s*(\[.*?\]);", text, re.DOTALL)
    if not m:
        return None
    try:
        data = json.loads(m.group(2))
    except Exception:
        return None
    for record in data:
        if record.get("seq_id") == seq_id:
            for region in record.get("regions", []):
                if region.get("idx") == cluster_idx:
                    return f"{region.get('start')}-{region.get('end')}"
    return None

# ---------------------------
# GC helpers (region GBK primary, genome FASTA fallback)
# ---------------------------

def _gc_percent(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    gc = s.count("G") + s.count("C")
    return 100.0 * gc / len(s)


def _parse_orig_coords_from_region_gbk_text(region_gbk: Path):
    """
    Parse antiSMASH comment lines:
      Orig. start  :: 651602
      Orig. end    :: 692206
    Returns (orig_start, orig_end) as ints, or (None, None).
    """
    try:
        txt = region_gbk.read_text(errors="ignore")
    except Exception:
        return None, None

    m1 = re.search(r"Orig\.\s*start\s*::\s*(\d+)", txt)
    m2 = re.search(r"Orig\.\s*end\s*::\s*(\d+)", txt)
    if not (m1 and m2):
        return None, None
    return int(m1.group(1)), int(m2.group(1))


def gc_from_region_with_fallback(region_gbk: Path, genome_fasta_dir: Path | None = None):
    """
    Returns (length_nt, gc_percent, source)
      source in {"region_gbk", "genome_fasta_slice", "missing"}
    """
    region_record = None

    # Primary: region GBK sequence
    try:
        region_record = next(SeqIO.parse(str(region_gbk), "genbank"))
        region_seq = str(region_record.seq).upper()
        if region_seq:
            return len(region_seq), _gc_percent(region_seq), "region_gbk"
    except Exception:
        region_record = None

    # Fallback: slice original genome FASTA
    if not genome_fasta_dir:
        return 0, 0.0, "missing"

    orig_start, orig_end = _parse_orig_coords_from_region_gbk_text(region_gbk)
    if not (orig_start and orig_end):
        return 0, 0.0, "missing"

    accession = region_gbk.name.split(".region", 1)[0]
    genome_fa = genome_fasta_dir / f"{accession}.fasta"
    if not genome_fa.exists():
        return 0, 0.0, "missing"

    target_id = region_record.id if region_record is not None else None

    chosen = None
    for rec in SeqIO.parse(str(genome_fa), "fasta"):
        if target_id and (rec.id == target_id or rec.name == target_id or rec.description.split()[0] == target_id):
            chosen = rec
            break
        if chosen is None:
            chosen = rec  # fallback to first record if no match

    if chosen is None:
        return 0, 0.0, "missing"

    # 1-based inclusive -> python slice [start-1:end]
    start0 = max(orig_start - 1, 0)
    end0 = orig_end
    sliced = str(chosen.seq[start0:end0]).upper()
    if not sliced:
        return 0, 0.0, "missing"

    return len(sliced), _gc_percent(sliced), "genome_fasta_slice"


# ---------------------------
# Summary table helpers
# ---------------------------

def get_species(desc: str) -> str:
    words = desc.split()

    # Handle "xxx: ..." prefix (e.g., "MAG:")
    if words and words[0].endswith(":"):
        words = words[1:]

    if "uncultured" in words or "Uncultured" in words:
        idx = words.index("uncultured") if "uncultured" in words else words.index("Uncultured")
        words = words[idx+1:]

    if "[Organism:" in words:
        idx = words.index("[Organism:")
        words = words[idx+1:]
        words[1] = words[1].rstrip("]")

    genus = words[0] if len(words) > 0 else ""
    species = words[1] if len(words) > 1 else ""

    # Handle "sp." and "ATCC"
    if species == "sp." and len(words) > 2:
        species += " " + words[2]
        if len(words) > 3 and "ATCC" in words[2]:
            species += " " + words[3]
        return f"{genus} {species}".strip()

    # MAG/uncultured bacterium-style labels
    if species in {"bacterium", "archaeon"}:
        return f"{genus} {species}".strip()

    return f"{genus} {species}".strip()


def parse_accession_and_cluster(rec_id: str):
    """
    Supports both:
      - AM420293.1.region003  (protein_concat.faa style)
      - AM420293.1_c3         (older style)
      - AM420293.1_c003
    Returns (accession, cluster_num_str)
    """
    m = re.search(r"^(?P<acc>.+?)\.region(?P<num>\d+)$", rec_id)
    if m:
        return m.group("acc"), str(int(m.group("num")))  # normalize leading zeros

    m = re.search(r"^(?P<acc>.+?)_c(?P<num>\d+)$", rec_id)
    if m:
        return m.group("acc"), str(int(m.group("num")))

    return rec_id, "unknown"


def get_cluster_type_from_gbk(gbk_file: Path) -> str:
    def clean_products(values):
        products = set()
        for value in values:
            value = re.sub(r"\s+", " ", str(value)).strip()
            if value:
                products.add(value)
        return products

    try:
        candidate_products = set()
        protocluster_products = set()
        legacy_products = set()

        for record in SeqIO.parse(str(gbk_file), "genbank"):
            for feature in record.features:
                products = clean_products(feature.qualifiers.get("product", []))
                if not products:
                    continue
                if feature.type == "cand_cluster":
                    candidate_products.update(products)
                elif feature.type == "protocluster":
                    protocluster_products.update(products)
                elif feature.type == "cluster":
                    legacy_products.update(products)

        products = candidate_products or protocluster_products or legacy_products
        return ";".join(sorted(products))
    except Exception:
        return ""


def write_final_summary(records, kept_ids, accession_to_desc, accession_to_date,
                        antismash_dir, out_file, genome_fasta_dir=None,
                        redundant_groups=None, homolog_groups=None):
    header = [
        "accession", "description", "date", "species", "cluster_number",
        "coordinates", "length", "gc_content(%)",
        "num_KSs", "num_ATs", "num_transATs", "num_ACPs", "num_KRs",
        "num_DHs", "num_ERs", "num_TEs", "num_Cs", "num_As", "num_Ps",
        "cluster_type", "identical_members", "homolog_members"
    ]

    with open(out_file, "w") as fout:
        fout.write("\t".join(header) + "\n")

        rec_by_id = {r.id: r for r in records}
        missing_count = 0
        for cid in sorted(kept_ids):
            rec = rec_by_id.get(cid)
            if rec is None:
                missing_count += 1
                continue
            accession, cluster_num = parse_accession_and_cluster(cid)
            desc = accession_to_desc.get(accession, "")
            date = accession_to_date.get(accession, datetime.datetime.max)
            date_str = date.strftime("%Y-%m-%d") if date != datetime.datetime.max else "NA"
            species = get_species(desc)

            # Coordinates from regions.js (best-effort; can be blank if seq_id mismatch)
            coords = ""
            if cluster_num != "unknown":
                js_file = antismash_dir / accession / "regions.js"
                coords = get_cluster_coords(js_file, int(cluster_num), accession) or ""

            length = len(rec.seq)

            # Build region GBK path (antiSMASH folder uses .gbk)
            region_gbk = None
            if cluster_num != "unknown":
                region_gbk = antismash_dir / accession / f"{accession}.region{int(cluster_num):03d}.gbk"
                if not region_gbk.exists():
                    region_gbk = antismash_dir / accession / f"{accession}.1.region{int(cluster_num):03d}.gbk"

            # GC% from region GBK (primary) with genome FASTA fallback (optional)
            length_nt, gc_content, gc_source = (0, 0.0, "missing")
            if region_gbk is not None and region_gbk.exists():
                length_nt, gc_content, gc_source = gc_from_region_with_fallback(region_gbk, genome_fasta_dir)

            if region_gbk is not None and region_gbk.exists():
                num_KSs = count_domains_in_gbk(region_gbk, "PKS_KS")
                num_ATs = count_domains_in_gbk(region_gbk, "PKS_AT")
                num_transATs = count_domains_in_gbk(region_gbk, "Trans-AT_docking")
                # ACP/PP-binding variants
                num_ACPs = (
                    count_domains_in_gbk(region_gbk, "ACP") +
                    count_domains_in_gbk(region_gbk, "PKS_PP") +
                    count_domains_in_gbk(region_gbk, "PP_binding") +
                    count_domains_in_gbk(region_gbk, "PP-binding")
                )
                num_KRs = count_domains_in_gbk(region_gbk, "PKS_KR")
                num_DHs = count_domains_in_gbk(region_gbk, "PKS_DH")  # matches PKS_DH2 etc via substring
                num_ERs = count_domains_in_gbk(region_gbk, "PKS_ER")
                num_TEs = count_domains_in_gbk(region_gbk, "Thioesterase")
                num_Cs = count_domains_in_gbk(region_gbk, "Condensation_")
                num_As = count_domains_in_gbk(region_gbk, "AMP-binding")
                num_Ps = count_domains_in_gbk(region_gbk, "PCP")
                cluster_type = get_cluster_type_from_gbk(region_gbk)
            else:
                num_KSs = num_ATs = num_transATs = num_ACPs = num_KRs = 0
                num_DHs = num_ERs = num_TEs = num_Cs = num_As = num_Ps = 0
                cluster_type = ""

            # Identical members (redundant_groups)
            identical = []
            if redundant_groups:
                identical = [m for m in redundant_groups.get(cid, []) if m != cid]

            # Homolog members (homolog_groups)
            homologs = []
            if homolog_groups:
                homologs = [m for m in homolog_groups.get(cid, []) if m != cid and m not in identical]

            fields = [
                accession, desc, date_str, species, cluster_num,
                coords, str(length), f"{gc_content:.2f}",
                str(num_KSs), str(num_ATs), str(num_transATs), str(num_ACPs),
                str(num_KRs), str(num_DHs), str(num_ERs), str(num_TEs),
                str(num_Cs), str(num_As), str(num_Ps),
                cluster_type,
                ",".join(identical),
                ",".join(homologs)
            ]
            fout.write("\t".join(fields) + "\n")

        if missing_count:
            print(f"[WARN] {missing_count} kept IDs were not found in records_faa; not written.", flush=True)
