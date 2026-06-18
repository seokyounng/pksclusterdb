#!/usr/bin/env python3
import argparse
import re
import sys
import shutil
import subprocess
import datetime
from pathlib import Path
from collections import defaultdict
from write_summary import (
    build_accession_desc_dict,
    parse_dates,
    write_final_summary,
)
import csv

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


KS_PATTERNS = ["PKS_KS"]
ACP_PCP_PATTERNS = ["ACP", "PKS_PP", "PP-binding", "PCP"]

# ---------------------------
# Command helpers
# ---------------------------

def run_cmd(cmd):
    """Run a command, echoing it to stderr for provenance."""
    sys.stderr.write("[CMD] " + " ".join(map(str, cmd)) + "\n")
    subprocess.run(cmd, check=True)


def safe_rmdb(prefix: Path):
    """Remove mmseqs DB leftovers (best-effort)."""
    prefix = Path(prefix)
    try:
        subprocess.run(["mmseqs", "rmdb", str(prefix)], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass
    parent = prefix.parent
    if parent.exists():
        for f in parent.iterdir():
            if f.name == prefix.name or f.name.startswith(prefix.name + "_") or f.name.startswith(prefix.name + "."):
                try:
                    shutil.rmtree(f, ignore_errors=True) if f.is_dir() else f.unlink()
                except Exception:
                    pass


# --------------------------------
# Protein sequence concatenation
# --------------------------------

def count_domains(gbk_file: Path, patterns) -> int:
    """
    Count aSDomain features where any qualifier value contains any of the patterns.
    Counts each aSDomain feature at most once.
    """
    patterns = list(patterns)
    count = 0
    try:
        for record in SeqIO.parse(str(gbk_file), "genbank"):
            for feature in record.features:
                if feature.type != "aSDomain":
                    continue
                vals = feature.qualifiers.get("aSDomain", [])
                if not vals:
                    continue
                hit = False
                for v in vals:
                    for p in patterns:
                        if p in v:
                            hit = True
                            break
                    if hit:
                        break
                if hit:
                    count += 1
    except Exception:
        return 0
    return count


def concat_region_proteins(region_gbk: Path, reverse_minus_blocks: bool = False):
    """
    Extract all CDS translations from an antiSMASH region GBK and concatenate them.

    Returns:
      concat_seq (str): concatenated protein sequence
      ordered_cds (list[dict]): metadata for each CDS in the concatenation order
    """
    it = SeqIO.parse(str(region_gbk), "genbank")
    record = next(it, None)
    if record is None:
        sys.stderr.write(f"[WARN] No GenBank records parsed from: {region_gbk}\n")
        return "", []

    cds_feats = []
    for f in record.features:
        if f.type != "CDS":
            continue
        tr = f.qualifiers.get("translation")
        if not tr:
            continue
        aa = tr[0].replace(" ", "").replace("\n", "")
        if not aa:
            continue

        start = int(f.location.start)
        end = int(f.location.end)
        strand = int(f.location.strand or 0)

        if "locus_tag" in f.qualifiers:
            locus = f.qualifiers["locus_tag"][0]
        elif "gene" in f.qualifiers:
            locus = f.qualifiers["gene"][0]
        elif "protein_id" in f.qualifiers:
            locus = f.qualifiers["protein_id"][0]
        else:
            locus = "unknown"

        cds_feats.append({
            "locus_tag": locus,
            "start": start,
            "end": end,
            "strand": strand,
            "aa": aa,
        })

    # Sort by genomic coordinate (increasing)
    cds_feats.sort(key=lambda x: (min(x["start"], x["end"]), max(x["start"], x["end"])))

    if reverse_minus_blocks and cds_feats:
        # Split into contiguous same-strand blocks in genomic order,
        # and reverse only the '-' blocks.
        blocks = []
        cur = [cds_feats[0]]
        for item in cds_feats[1:]:
            if item["strand"] == cur[-1]["strand"]:
                cur.append(item)
            else:
                blocks.append(cur)
                cur = [item]
        blocks.append(cur)

        new_order = []
        for blk in blocks:
            if blk and blk[0]["strand"] == -1:
                new_order.extend(list(reversed(blk)))
            else:
                new_order.extend(blk)
        cds_feats = new_order

    concat_seq = "".join(x["aa"] for x in cds_feats)
    return concat_seq, cds_feats


def generate_protein_concat_faa(
    antismash_dir: Path,
    out_faa: Path,
    reverse_minus_blocks: bool = False,
    min_ks: int = 3,
    min_acp_pcp: int = 3,
):
    """
    Create one combined protein_concat.faa across all accession subdirectories under antismash_dir.
    Each record id: <accession>.region###
    Only include regions with >= min_ks PKS_KS domains and >= min_acp_pcp ACP/PCP-like domains.
    """
    tmp_records = []
    skipped_domain_filter = 0
    skipped_empty_or_bad = 0
    kept = 0

    for acc_dir in sorted(p for p in antismash_dir.iterdir() if p.is_dir()):
        accession = acc_dir.name
        region_gbks = sorted(acc_dir.glob("*.region*.gbk"))
        for region_gbk in region_gbks:
            m = re.search(r"\.region(\d+)\.gbk$", region_gbk.name)
            if not m:
                continue
            region_num = int(m.group(1))
            rec_id = f"{accession}.region{region_num:03d}"

            # --- Domain filter (BEFORE concat) ---
            ks_count = count_domains(region_gbk, KS_PATTERNS)
            acp_count = count_domains(region_gbk, ACP_PCP_PATTERNS)
            if ks_count < min_ks or acp_count < min_acp_pcp:
                skipped_domain_filter += 1
                continue

            # --- Concat translations ---
            concat_seq, cds_feats = concat_region_proteins(region_gbk, reverse_minus_blocks=reverse_minus_blocks)
            if not concat_seq:
                skipped_empty_or_bad += 1
                continue

            header = f"{rec_id} locus_tags=" + ",".join(x["locus_tag"] for x in cds_feats)
            tmp_records.append(SeqRecord(Seq(concat_seq), id=rec_id, description=header))
            kept += 1

    out_faa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_faa, "w") as fh:
        SeqIO.write(tmp_records, fh, "fasta")

    sys.stderr.write(
        f"[INFO] protein_concat.faa: kept={kept}, skipped_domain_filter={skipped_domain_filter}, "
        f"skipped_empty_or_bad={skipped_empty_or_bad}\n"
    )


# ---------------------------
# mmseqs clustering helpers
# ---------------------------

def parse_mmseqs_clusters(cluster_tsv: Path):
    """
    Parse mmseqs *_cluster.tsv into rep -> members dict.
    Expected format: representative<TAB>member per line.
    Note: representatives may or may not appear as their own member; we treat rep as a member of its cluster implicitly.
    """
    groups = defaultdict(set)
    with open(cluster_tsv) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, mem = parts[0], parts[1]
            if rep == mem:
                continue
            groups[rep].add(mem)
    return groups


def _mmseqs_expected_cluster_tsv(out_prefix: Path) -> Path:
    """
    MMseqs easy-cluster/easy-linclust typically writes '<out_prefix>_cluster.tsv'.
    Return that path if present; otherwise search nearby for '*cluster.tsv'.
    """
    out_prefix = Path(out_prefix)
    cand = Path(str(out_prefix) + "_cluster.tsv")
    if cand.exists():
        return cand
    # fallback search
    parent = out_prefix.parent if out_prefix.parent.exists() else Path(".")
    hits = sorted(parent.glob(out_prefix.name + "*cluster.tsv"))
    if hits:
        return hits[0]
    # last resort
    hits = sorted(parent.glob("*cluster.tsv"))
    if hits:
        return hits[0]
    raise FileNotFoundError(f"Could not find mmseqs cluster TSV for prefix: {out_prefix}")


def sort_fasta_by_length_desc(in_faa: Path, out_faa: Path) -> None:
    """Write FASTA records sorted by decreasing sequence length (stable by id for ties)."""
    recs = list(SeqIO.parse(str(in_faa), "fasta"))
    recs.sort(key=lambda r: (-len(r.seq), r.id))
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_faa, "w") as f:
        SeqIO.write(recs, f, "fasta")


def run_mmseqs_easy_cluster(
    in_faa: Path,
    out_prefix: Path,
    tmpdir: Path,
    *,
    min_seq_id: float,
    coverage: float,
    cov_mode: int,
    threads: int,
):
    """
    Run mmseqs easy-cluster.
      - cov-mode (0/1/2 are coverage-based; 3/4/5 are length-ratio modes).
    """
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "mmseqs", "easy-cluster",
        str(in_faa), str(out_prefix), str(tmpdir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--alignment-mode", "2",
        "--cluster-mode", "2",     # greedy by sequence length (CDHIT-like)
        "-s", "4.0",
        "--max-seqs", "300",
        "--threads", str(threads),
        "-v", "3",
    ]
    run_cmd(cmd)
    return _mmseqs_expected_cluster_tsv(out_prefix)


def run_mmseqs_easy_linclust(
    in_faa: Path,
    out_prefix: Path,
    tmpdir: Path,
    *,
    min_seq_id: float,
    coverage: float,
    cov_mode: int,
    threads: int,
):
    """
    Run mmseqs easy-linclust (linear-time clustering for large sets).
    """
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "mmseqs", "easy-linclust",
        str(in_faa), str(out_prefix), str(tmpdir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--alignment-mode", "3",
        "--cluster-mode", "2",     # greedy by sequence length (CDHIT-like)
        "--threads", str(threads),
        "-v", "3",
    ]
    run_cmd(cmd)
    return _mmseqs_expected_cluster_tsv(out_prefix)


class UnionFind:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x):
        p = self.parent.get(x, x)
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent.get(x, x)

    def union(self, a, b):
        self.add(a)
        self.add(b)
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def _record_accession(rec_id: str) -> str:
    # record ids are <accession>.region###
    return rec_id.split(".region", 1)[0]


def _dominates(a: str, b: str, lengths: dict, acc_dates: dict) -> bool:
    """True if a should be the representative over b (length desc, date asc, id asc)."""
    la, lb = lengths.get(a, 0), lengths.get(b, 0)
    if la != lb:
        return la > lb
    da = acc_dates.get(_record_accession(a), datetime.datetime.max)
    db = acc_dates.get(_record_accession(b), datetime.datetime.max)
    if da != db:
        return da < db
    return a < b


def stage1_dedup(
    records: list,
    records_faa: Path,
    out_dir: Path,
    *,
    accession_date: dict,
    threads: int,
    min_seq_id: float,
    coverage: float,
):
    """
    Stage 1: remove exact/near-exact duplicates and absorb obvious fragments.

    Implementation:
      - Sort FASTA by length descending.
      - Run mmseqs easy-linclust with target-coverage filtering (cov-mode=1) so
        fragments can be absorbed into longer representatives.
      - Post-process cluster TSV to ensure representative selection is stable:
          choose rep by _dominates() among each connected component.
    """
    mmseqs_dir = out_dir / "mmseqs_results"
    mmseqs_dir.mkdir(parents=True, exist_ok=True)

    sorted_faa = mmseqs_dir / "protein_concat.sorted.faa"
    sort_fasta_by_length_desc(records_faa, sorted_faa)

    tmpdir = out_dir / "mmseqs_tmp_dedup"
    out_prefix = mmseqs_dir / "dedup"
    cluster_tsv = run_mmseqs_easy_linclust(
        sorted_faa,
        out_prefix,
        tmpdir,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=1,  # target coverage
        threads=threads,
    )

    lengths = {r.id: len(r.seq) for r in records}
    uf = UnionFind()
    for rid in lengths:
        uf.add(rid)

    # Build components from rep-member edges.
    groups = parse_mmseqs_clusters(cluster_tsv)
    for rep, mems in groups.items():
        if rep in lengths:
            uf.add(rep)
        for m in mems:
            if m in lengths:
                uf.union(rep, m)

    comp = defaultdict(set)
    for rid in lengths:
        comp[uf.find(rid)].add(rid)

    kept = set()
    dedup_map = defaultdict(set)
    for _, members in comp.items():
        if len(members) == 1:
            m = next(iter(members))
            kept.add(m)
            continue
        rep = None
        for m in members:
            if rep is None or _dominates(m, rep, lengths, accession_date):
                rep = m
        kept.add(rep)
        for m in members:
            if m != rep:
                dedup_map[rep].add(m)

    # Artifacts
    with open(mmseqs_dir / "dedup_map.tsv", "w") as f:
        for rep, mems in sorted(dedup_map.items()):
            for m in sorted(mems):
                f.write(f"{rep}\t{m}\n")
    with open(mmseqs_dir / "dedup_kept_ids.txt", "w") as f:
        for rid in sorted(kept):
            f.write(rid + "\n")

    return kept, dedup_map


# ---------------------------
# PKS-only concatenation (for stage 2)
# ---------------------------

PKS_DOMAINS = [
    "PKS_KS", "PKS_AT", "Trans-AT_docking", "ACP", "PKS_PP",
    "PKS_KR", "PKS_DH", "PKS_DH2", "PKS_DHt", "PKS_ER", "Thioesterase"
]


def _intervals_overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return not (a1 <= b0 or b1 <= a0)


def _alphabetic_gene_key(locus_tag: str):
    """
    Return (prefix, suffix) for characterized gene names like fosA/herB.

    These names do not encode genomic coordinates, so within a same-prefix
    family their biological order should be A, B, C, ... rather than the
    coordinate order reported in the region file.
    """
    m = re.match(r"^([A-Za-z]+)([A-Z])$", locus_tag)
    if not m:
        return None
    return m.group(1), m.group(2)


def sort_alphabetic_gene_blocks(cds_feats: list[dict]) -> list[dict]:
    if not cds_feats:
        return cds_feats

    groups = defaultdict(list)
    for idx, item in enumerate(cds_feats):
        key = _alphabetic_gene_key(item["locus_tag"])
        if key is None:
            continue
        prefix, suffix = key
        groups[prefix].append((idx, suffix, item))

    reordered = list(cds_feats)
    for items in groups.values():
        if len(items) < 2:
            continue
        idxs = [x[0] for x in items]
        sorted_items = [x[2] for x in sorted(items, key=lambda x: x[1])]
        for idx, item in zip(sorted(idxs), sorted_items):
            reordered[idx] = item

    return reordered


def concat_region_pks_proteins(region_gbk: Path, reverse_minus_blocks: bool = False):
    """
    Extract CDS translations that overlap any antiSMASH aSDomain whose aSDomain qualifier
    contains one of PKS_DOMAINS, then concatenate those CDS translations in genomic order.

    Returns:
      concat_seq (str), ordered_cds (list[dict])  # same schema as concat_region_proteins
    """
    it = SeqIO.parse(str(region_gbk), "genbank")
    record = next(it, None)
    if record is None:
        sys.stderr.write(f"[WARN] No GenBank records parsed from: {region_gbk}\n")
        return "", []

    # Collect CDS features
    cds_feats = []
    for f in record.features:
        if f.type != "CDS":
            continue
        tr = f.qualifiers.get("translation")
        if not tr:
            continue
        aa = tr[0].replace(" ", "").replace("\n", "")
        if not aa:
            continue

        start = int(f.location.start)
        end = int(f.location.end)
        strand = int(f.location.strand or 0)

        if "locus_tag" in f.qualifiers:
            locus = f.qualifiers["locus_tag"][0]
        elif "gene" in f.qualifiers:
            locus = f.qualifiers["gene"][0]
        elif "protein_id" in f.qualifiers:
            locus = f.qualifiers["protein_id"][0]
        else:
            locus = "unknown"

        cds_feats.append({"locus_tag": locus, "start": start, "end": end, "strand": strand, "aa": aa})

    if not cds_feats:
        return "", []

    # Collect PKS domain intervals
    dom_intervals = []
    for f in record.features:
        if f.type != "aSDomain":
            continue
        vals = f.qualifiers.get("aSDomain", [])
        if not vals:
            continue
        hit = any(any(tok in v for tok in PKS_DOMAINS) for v in vals)
        if not hit:
            continue
        dom_intervals.append((int(f.location.start), int(f.location.end)))

    if not dom_intervals:
        return "", []

    # Select CDS that overlap any PKS domain interval
    keep = []
    for cds in cds_feats:
        s0, s1 = min(cds["start"], cds["end"]), max(cds["start"], cds["end"])
        if any(_intervals_overlap(s0, s1, d0, d1) for d0, d1 in dom_intervals):
            keep.append(cds)

    if not keep:
        return "", []

    # Sort by genomic coordinate (increasing)
    keep.sort(key=lambda x: (min(x["start"], x["end"]), max(x["start"], x["end"])))

    if reverse_minus_blocks and keep:
        blocks = []
        cur = [keep[0]]
        for item in keep[1:]:
            if item["strand"] == cur[-1]["strand"]:
                cur.append(item)
            else:
                blocks.append(cur)
                cur = [item]
        blocks.append(cur)

        new_order = []
        for blk in blocks:
            if blk and blk[0]["strand"] == -1:
                new_order.extend(list(reversed(blk)))
            else:
                new_order.extend(blk)
        keep = new_order

    keep = sort_alphabetic_gene_blocks(keep)

    concat_seq = "".join(x["aa"] for x in keep)
    return concat_seq, keep


def _id_to_region_gbk(antismash_dir: Path, rec_id: str) -> Path:
    """
    rec_id: <accession>.region###
    returns: path to <antismash_dir>/<accession>/*.region###.gbk
    """
    accession, region_part = rec_id.split(".region", 1)
    region_num = int(region_part)
    acc_dir = antismash_dir / accession
    # antiSMASH naming: <accession>.region###.gbk (may have extra prefixes); use glob
    pats = [
        f"*.region{region_num}.gbk",
        f"*.region{region_num:03d}.gbk",
        f"{accession}.region{region_num}.gbk",
        f"{accession}.region{region_num:03d}.gbk",
    ]
    for pat in pats:
        hits = sorted(acc_dir.glob(pat)) if acc_dir.exists() else []
        if hits:
            return hits[0]
    # fallback: search any region gbk with that number
    hits = sorted(acc_dir.glob(f"*.region{region_num}*.gbk")) if acc_dir.exists() else []
    if hits:
        return hits[0]
    return Path("")


def generate_pks_concat_faa_for_ids(
    antismash_dir: Path,
    rec_ids: set,
    out_faa: Path,
    reverse_minus_blocks: bool = False,
):
    """
    Build a FASTA where each record is a concatenation of PKS-related CDS translations
    (as detected by overlapping aSDomain features with PKS_DOMAINS).
    Only IDs present in rec_ids are attempted.
    Returns: (written_ids_set)
    """
    tmp_records = []
    written = set()
    missing = 0
    empty = 0

    for rec_id in sorted(rec_ids):
        region_gbk = _id_to_region_gbk(antismash_dir, rec_id)
        if not region_gbk or not region_gbk.exists():
            missing += 1
            continue
        concat_seq, cds_feats = concat_region_pks_proteins(region_gbk, reverse_minus_blocks=reverse_minus_blocks)
        if not concat_seq:
            empty += 1
            continue
        header = f"{rec_id} locus_tags=" + ",".join(x["locus_tag"] for x in cds_feats)
        tmp_records.append(SeqRecord(Seq(concat_seq), id=rec_id, description=header))
        written.add(rec_id)

    out_faa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_faa, "w") as fh:
        SeqIO.write(tmp_records, fh, "fasta")

    sys.stderr.write(
        f"[INFO] pks_concat.faa: written={len(written)}, missing_gbk={missing}, empty_pks_concat={empty}\n"
    )
    return written


def write_sorted_subset_faa(in_faa: Path, keep_ids: set[str], out_sorted_faa: Path) -> None:
    """Write only keep_ids from in_faa to out_sorted_faa, sorted by length desc."""
    keep_ids = set(keep_ids)
    recs = [r for r in SeqIO.parse(str(in_faa), "fasta") if r.id in keep_ids]
    recs.sort(key=lambda r: (-len(r.seq), r.id))
    out_sorted_faa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_sorted_faa, "w") as f:
        SeqIO.write(recs, f, "fasta")


def stage2_collapse_homologs(
    records: list,
    pks_faa: Path,
    out_dir: Path,
    *,
    universe_ids: set,
    accession_date: dict,
    # 2a (fragment absorption)
    a_min_seq_id: float = 0.95,
    a_coverage: float = 0.90,
    # 2b (full-length collapse)
    b_min_seq_id: float = 0.95,
    b_coverage: float = 0.80,
    threads: int,
):
    """
    Stage 2 (two-step):

      2a) Fragment absorption on PKS concatenations:
          - cov-mode=1 (target coverage), higher coverage (e.g., 0.90)
          - absorbs partials into longer reps if they match well.

      2b) Full-length homolog collapse among 2a-kept:
          - cov-mode=0 (bidirectional), coverage ~0.80
          - collapses near-identical full-length homologs.

    Inputs:
      - pks_faa should already be built from dedup_kept (your pipeline does that).
      - universe_ids should usually be dedup_kept (or dedup_kept ∩ written_pks_ids).
    """
    out_dir = Path(out_dir)
    mmseqs_dir = out_dir / "mmseqs_results"
    mmseqs_dir.mkdir(parents=True, exist_ok=True)

    universe_ids = set(universe_ids)

    # Full-length sizes (for rep choice) restricted to universe only
    full_lengths = {r.id: len(r.seq) for r in records if r.id in universe_ids}

    # PKS ids actually present
    pks_ids_all = {r.id for r in SeqIO.parse(str(pks_faa), "fasta")}
    pks_ids = pks_ids_all & universe_ids

    if not pks_ids:
        sys.stderr.write("[WARN] Stage2: no PKS sequences in universe; skipping.\n")
        return set(universe_ids), defaultdict(set)

    # --- Make the canonical sorted input for Stage 2a from *pks_faa* (already dedup-restricted) ---
    pks_sorted = mmseqs_dir / "pks_concat.sorted.faa"
    # Only sort once, but ensure we don't carry non-universe IDs:
    write_sorted_subset_faa(pks_faa, pks_ids, pks_sorted)

    def _clusters_to_components(cluster_tsv: Path, allowed_ids: set[str]) -> list[set[str]]:
        """Build connected components from rep-member edges, but only within allowed_ids."""
        uf = UnionFind()
        for rid in allowed_ids:
            uf.add(rid)

        groups = parse_mmseqs_clusters(cluster_tsv)
        for rep, mems in groups.items():
            if rep not in allowed_ids:
                continue
            for m in mems:
                if m in allowed_ids:
                    uf.union(rep, m)

        comp = defaultdict(set)
        for rid in allowed_ids:
            comp[uf.find(rid)].add(rid)
        return list(comp.values())

    def _pick_reps_and_map(components: list[set[str]]) -> tuple[set[str], dict]:
        kept = set()
        repmap = defaultdict(set)
        for members in components:
            if len(members) == 1:
                kept.add(next(iter(members)))
                continue
            rep = None
            for m in members:
                if rep is None or _dominates(m, rep, full_lengths, accession_date):
                    rep = m
            kept.add(rep)
            for m in members:
                if m != rep:
                    repmap[rep].add(m)
        return kept, repmap

    # =====================
    # 2a) Fragment absorption
    # =====================
    tmp_a = out_dir / "mmseqs_tmp_homolog_a"
    out_a = mmseqs_dir / "homolog_a_fragabsorb"
    tsv_a = run_mmseqs_easy_cluster(
        pks_sorted,
        out_a,
        tmp_a,
        min_seq_id=a_min_seq_id,
        coverage=a_coverage,
        cov_mode=1,
        threads=threads,
    )
    comps_a = _clusters_to_components(tsv_a, pks_ids)
    kept_a, map_a = _pick_reps_and_map(comps_a)

    print(f"[DEBUG] Stage 2a: {len(kept_a)} kept after fragment absorption, {sum(len(mems) for mems in map_a.values())} absorbed as fragments.", flush=True)

    # =====================
    # 2b) Full-length collapse among 2a-kept
    # =====================
    pks_a_sorted = mmseqs_dir / "pks_concat.a_kept.sorted.faa"
    write_sorted_subset_faa(pks_faa, kept_a, pks_a_sorted)

    tmp_b = out_dir / "mmseqs_tmp_homolog_b"
    out_b = mmseqs_dir / "homolog_b_fullcollapse"
    tsv_b = run_mmseqs_easy_cluster(
        pks_a_sorted,
        out_b,
        tmp_b,
        min_seq_id=b_min_seq_id,
        coverage=b_coverage,
        cov_mode=0,
        threads=threads,
    )
    comps_b = _clusters_to_components(tsv_b, kept_a)
    kept_b, map_b = _pick_reps_and_map(comps_b)

    print(f"[DEBUG] Stage 2b: {len(kept_b)} kept after full-length collapse, {sum(len(mems) for mems in map_b.values())} collapsed as homologs.", flush=True)

    # =====================
    # Merge 2a + 2b maps transitively
    # =====================
    final_map = defaultdict(set)
    # Start with 2a absorption
    for rep, mems in map_a.items():
        final_map[rep].update(mems)

    # Then apply 2b collapse; if 2b absorbs a rep that already absorbed fragments in 2a,
    # carry those members over to the 2b winner.
    for rep, mems in map_b.items():
        for m in mems:
            final_map[rep].add(m)
            final_map[rep].update(map_a.get(m, set()))
        final_map[rep].update(map_a.get(rep, set()))

    # Artifacts
    with open(mmseqs_dir / "homolog_map.tsv", "w") as f:
        for rep, mems in sorted(final_map.items()):
            for m in sorted(mems):
                f.write(f"{rep}\t{m}\n")
    with open(mmseqs_dir / "homolog_kept_ids.txt", "w") as f:
        for rid in sorted(kept_b):
            f.write(rid + "\n")

    return kept_b, final_map


def remap_dedup_to_final_reps(dedup_map, homolog_map, kept_ids):
    """
    Keep only stage1 dedup members that belong directly to final surviving reps.

    If a stage1 rep was later removed in stage2, its dedup members should not be
    reclassified as identical_members of the final rep; they should instead travel
    with that removed rep in homolog_members.
    """
    final_dedup = defaultdict(set)

    for rep in kept_ids:
        final_dedup[rep].update(dedup_map.get(rep, set()))

    return final_dedup


def merge_nz_pairs_in_final_summary(final_tsv: Path):
    """
    Post-process final_summary.tsv:
      - If accession 'NZ_xx' exists and non-NZ accession 'xx' also exists
      - and |cluster_number difference| <= 2
      - and gc_content(%) is identical
    then merge the NZ row into the non-NZ row's identical_members
    and remove the NZ row from the table.

    Important: direction is NZ_xx -> xx, not the reverse lookup.
    """
    csv.field_size_limit(sys.maxsize)

    with open(final_tsv, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
        fieldnames = rows[0].keys() if rows else []

    # index rows by accession
    rows_by_accession = defaultdict(list)
    for row in rows:
        rows_by_accession[row["accession"]].append(row)

    # rows to remove after merging
    remove_keys = set()

    def norm_gc(x):
        # exact-ish numeric comparison; safer than raw string compare
        try:
            return f"{float(x):.6f}"
        except Exception:
            return str(x).strip()

    def make_cluster_id(row):
        # Match your summary's cluster identifier style as needed.
        # If identical_members already uses accession.region###, keep this.
        try:
            n = int(row["cluster_number"])
            return f'{row["accession"]}.region{n:03d}'
        except Exception:
            return f'{row["accession"]}.region{row["cluster_number"]}'

    for nz_acc, nz_rows in rows_by_accession.items():
        if not nz_acc.startswith("NZ_"):
            continue

        base_acc = nz_acc[3:]   # NZ_xx -> xx
        if base_acc not in rows_by_accession:
            continue

        base_rows = rows_by_accession[base_acc]

        for nz_row in nz_rows:
            try:
                nz_cluster = int(nz_row["cluster_number"])
            except Exception:
                continue
            nz_gc = norm_gc(nz_row["gc_content(%)"])

            matched_base = None
            for base_row in base_rows:
                try:
                    base_cluster = int(base_row["cluster_number"])
                except Exception:
                    continue
                base_gc = norm_gc(base_row["gc_content(%)"])

                if abs(base_cluster - nz_cluster) <= 2 and base_gc == nz_gc:
                    matched_base = base_row
                    break

            if matched_base is None:
                continue

            nz_cluster_id = make_cluster_id(nz_row)

            existing = [x for x in (matched_base.get("identical_members") or "").split(",") if x]
            if nz_cluster_id not in existing:
                existing.append(nz_cluster_id)
            matched_base["identical_members"] = ",".join(existing)

            remove_keys.add((nz_row["accession"], nz_row["cluster_number"]))

    # keep original row order, except merged-away NZ rows
    kept_rows = [
        row for row in rows
        if (row["accession"], row["cluster_number"]) not in remove_keys
    ]

    with open(final_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept_rows)

# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_ks", type=int, default=3)
    parser.add_argument("--min_acp_pcp", type=int, default=3)
    parser.add_argument("--concat_faa", default=None,
                        help="FASTA of AA sequences; if omitted, protein_concat.faa is generated")
    parser.add_argument("--pks_concat_faa", default=None,
                    help="FASTA of AA sequences; if omitted, pks_concat.faa is generated")
    parser.add_argument("--desc", required=True)
    parser.add_argument("--dates", required=True)
    parser.add_argument("--antismash_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--genome_fasta_dir", default=None, help="Optional: directory of whole-genome FASTAs for GC fallback")
    # Optional: provide precomputed mapping TSVs (rep<TAB>member). If omitted, we run mmseqs easy-search.
    parser.add_argument("--dedup_clusters", default=None, help="Optional: rep<TAB>member map to use for dedup step")
    parser.add_argument("--homolog_clusters", default=None, help="Optional: rep<TAB>member map to use for homolog collapse")

    # Stage 1 (dedup + fragment absorption) parameters
    parser.add_argument("--dedup_threads", type=int, default=8)
    parser.add_argument("--dedup_min_seq_id", type=float, default=0.99,
                        help="mmseqs --min-seq-id for stage1 (dedup/fragment absorption). Typical: 0.98-1.0")
    parser.add_argument("--dedup_coverage", type=float, default=0.90,
                        help="mmseqs -c coverage threshold for stage1 with --cov-mode 1 (target coverage). Typical: 0.90-0.95")

    # Stage 2 (homolog collapse) parameters
    parser.add_argument("--homolog_min_seq_id", type=float, default=0.95)
    parser.add_argument("--homolog_coverage", type=float, default=0.80,
                        help="mmseqs -c coverage threshold for stage2 with --cov-mode 0 (bidirectional). Typical: 0.70-0.85")
    parser.add_argument("--homolog_threads", type=int, default=8)
    args = parser.parse_args()

    accession_desc = build_accession_desc_dict(args.desc)
    accession_date = parse_dates(args.dates)
    antismash_dir = Path(args.antismash_dir)
    output_dir = Path(args.out_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    if args.concat_faa:
        records_faa = Path(args.concat_faa)
    else:
        print("Generating protein_concat.faa...", flush=True)
        records_faa = output_dir / "protein_concat.faa"
        generate_protein_concat_faa(
            antismash_dir,
            records_faa,
            reverse_minus_blocks=True,
            min_ks=args.min_ks,
            min_acp_pcp=args.min_acp_pcp,
        )

    records = list(SeqIO.parse(records_faa, "fasta"))

    genome_fasta_dir = Path(args.genome_fasta_dir) if args.genome_fasta_dir else None

    # ---------------------------
    # Stage 1: Deduplication
    # ---------------------------
    print("Stage 1: Deduplication...", flush=True)

    if args.dedup_clusters:
        mmseqs_dir = Path(args.dedup_clusters)
        if not mmseqs_dir.is_dir():
            raise ValueError(
                "--dedup_clusters must be the mmseqs_results directory "
                "containing dedup_map.tsv and dedup_kept_ids.txt"
            )

        map_file = mmseqs_dir / "dedup_map.tsv"
        kept_file = mmseqs_dir / "dedup_kept_ids.txt"
        if not map_file.exists():
            raise FileNotFoundError(f"Missing: {map_file}")
        if not kept_file.exists():
            raise FileNotFoundError(f"Missing: {kept_file}")

        # Load rep->member mapping
        dedup_map = parse_mmseqs_clusters(map_file)
        # Load true kept representatives
        with open(kept_file) as f:
            dedup_kept = {line.strip() for line in f if line.strip()}

    else:
        dedup_kept, dedup_map = stage1_dedup(
            records,
            records_faa,
            output_dir,
            accession_date=accession_date,
            threads=args.dedup_threads,
            min_seq_id=args.dedup_min_seq_id,
            coverage=args.dedup_coverage,
        )

    print(f"Deduplicated from {len(records)} to {len(dedup_kept)} sequences.", flush=True)

    print("Stage 2: Collapsing homologs...", flush=True)

    if args.homolog_clusters:
        mmseqs_dir = Path(args.homolog_clusters)
        if not mmseqs_dir.is_dir():
            raise ValueError(
                "--homolog_clusters must be the mmseqs_results directory "
                "containing homolog_map.tsv and homolog_kept_ids.txt"
            )

        map_file = mmseqs_dir / "homolog_map.tsv"
        kept_file = mmseqs_dir / "homolog_kept_ids.txt"
        if not map_file.exists():
            raise FileNotFoundError(f"Missing: {map_file}")
        if not kept_file.exists():
            raise FileNotFoundError(f"Missing: {kept_file}")

        homolog_map = parse_mmseqs_clusters(map_file)
        with open(kept_file) as f:
            homolog_kept = {line.strip() for line in f if line.strip()}

    else:
        # Build PKS-only concatenations for deduplicated IDs (used for homolog collapsing).
        if args.pks_concat_faa:
            pks_faa = Path(args.pks_concat_faa)
            written_pks_ids = set(r.id for r in SeqIO.parse(str(pks_faa), "fasta"))
        else:
            print("Generating pks_concat.faa...", flush=True)
            pks_faa = output_dir / "pks_concat.faa"
            written_pks_ids = generate_pks_concat_faa_for_ids(
                antismash_dir,
                dedup_kept,
                pks_faa,
                reverse_minus_blocks=True,
            )

        if not written_pks_ids:
            sys.stderr.write("[WARN] No PKS-only concatenations were generated; skipping stage2 clustering.\n")
            homolog_kept, homolog_map = set(dedup_kept), defaultdict(set)
        else:
            homolog_kept, homolog_map = stage2_collapse_homologs(
                records,
                pks_faa,
                output_dir,
                universe_ids=dedup_kept,
                accession_date=accession_date,
                b_min_seq_id=args.homolog_min_seq_id,
                b_coverage=args.homolog_coverage,
                threads=args.homolog_threads,
            )

    kept_ids = set(homolog_kept)

    redundant_groups = remap_dedup_to_final_reps(dedup_map, homolog_map, kept_ids)

    homolog_groups = defaultdict(set)
    for final_rep, mems in homolog_map.items():
        homolog_groups[final_rep].update(mems)
        for rid in mems:
            homolog_groups[final_rep].update(dedup_map.get(rid, set()))

    final_tsv = output_dir / "final_summary.tsv"
    write_final_summary(
        records,
        kept_ids,
        accession_desc,
        accession_date,
        antismash_dir,
        final_tsv,
        genome_fasta_dir=genome_fasta_dir,
        redundant_groups=redundant_groups,
        homolog_groups=homolog_groups,
    )
    merge_nz_pairs_in_final_summary(final_tsv)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
