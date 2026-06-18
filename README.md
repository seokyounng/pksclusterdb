# PKSClusterDB

This repository contains the construction workflow for **PKSClusterDB**, a curated catalogue of non-redundant assembly-line polyketide synthase (PKS) biosynthetic gene clusters. The workflow was developed for the review article **“Genomic Scale Analysis of Assembly-Line Polyketide Synthase Diversity and Evolution”** by Seokyoung Lee and Chaitan Khosla. A publication link will be added here when available.

PKSClusterDB was designed to organize the rapidly expanding sequence space of bacterial assembly-line PKSs, including both chemically characterized systems and orphan biosynthetic gene clusters. The workflow starts from KS-domain homology searches against NCBI BLAST/WGS resources, filters candidate loci by clustered KS-domain hits (HSPs), annotates candidate accessions with antiSMASH, removes redundancy with a staged MMseqs2-based procedure, and optionally annotates overlap with known MiBIG biosynthetic gene clusters.

In the manuscript-associated build, this workflow yielded **16,633 non-redundant assembly-line PKSs**. The resulting catalogue was used as the sequence foundation for family-resolved comparative analyses, including anchor-window-based interrogation of conserved multimodular regions. This repository covers the **database construction workflow**; downstream evolutionary or chemical analyses may require additional scripts or manual curation steps described in the paper.

<br>

## Scope and reproducibility notes

This workflow is intended for large-scale reconstruction of the PKSClusterDB build. It uses public NCBI databases, WGS/SRA resources, antiSMASH annotations, MMseqs2 clustering, and MiBIG matching. Because NCBI records and WGS databases change over time, an exact rerun may not produce identical accession counts unless the same database snapshots are used.

The workflow includes automated steps and curation-oriented outputs. Users should inspect intermediate files and final summaries, especially when applying the pipeline to newer database snapshots, alternative antiSMASH versions, or different redundancy thresholds.

<br>

## Citation

If you use this workflow or the resulting PKSClusterDB catalogue, please cite:

> Lee, S. and Khosla, C. **Genomic Scale Analysis of Assembly-Line Polyketide Synthase Diversity and Evolution.** *Natural Product Reports*. Submitted.

Please also cite the external tools used in the workflow, including antiSMASH, MMseqs2, NCBI BLAST+, SRA Toolkit, and MiBIG, as appropriate.

<br>

## License

This code is released under the MIT License. See `LICENSE` for details.

<br>

## Requirements

Install dependencies manually with conda, pip, modules, or your HPC software stack. The scripts do not install software automatically.

Python:

- Python >= 3.9
- biopython
- pandas

External tools:

- NCBI BLAST+ 2.16.0+
- SRA Toolkit 3.2.0 (`prefetch`, `tblastn_vdb`)
- GNU parallel 20200822
- antiSMASH 7.1
- MMseqs2 13.45111
- Perl 5
- Perl modules: `LWP::UserAgent`, `LWP::Protocol::https`

The shell examples use `bash` and do not require SLURM. On a cluster, submit the same commands with your local scheduler settings if needed.

<br>

## Inputs

Main required inputs:

- `queries/KSSignatureConsensusPKSDB.fasta`: consensus KS-domain query sequence used for initial `tblastn`/`tblastn_vdb` searches.
- NCBI BLAST databases, downloaded locally or available through an existing institutional installation.
- WGS/SRA accessions retrieved through SRA Toolkit for bacteria, archaea, and eukaryota WGS searches.
- MiBIG JSON and GBK archives for optional annotation of known biosynthetic gene clusters.

Several scripts require NCBI Entrez access. Set `NCBI_EMAIL` to a valid email address before fetching GenBank or FASTA records.

<br>

## Estimated storage

Storage requirements depend on the current NCBI/WGS database state and how many hits pass filtering. For the July 2025 database downloads used in this workflow:

- Initial NCBI BLAST databases plus WGS/SRA downloads: about 10 TB
- Fetched GenBank records for HSP-filtered entries: about 160 GB
- antiSMASH output for eligible candidate clusters: about 660 GB

The GenBank records fetched from `ref_euk_rep_genomes` and `ref_prok_rep_genomes` may not contain full sequence or annotated CDS features, so this workflow also downloads FASTA records for candidate accessions. `run_antismash.sh` uses the FASTA files as a fallback when a GenBank record is missing usable sequence.

Recommended project layout:

```text
pksclusterdb/
  queries/
    KSSignatureConsensusPKSDB.fasta
  data/
    ncbi_blastdbs/
    genbank_gb/
    genbank_fasta/
    mibig/
  outputs/
    tblastn/
    ks_filter/
    antismash/
    cleaned_clusters/
  log/
  scripts/
```

<br>

## Workflow overview

The workflow has two search inputs: preformatted NCBI BLAST databases and WGS/SRA accessions. Both write BLAST-tabular output under `outputs/tblastn/`, then merge at the KS-filtering step.

The major stages are:

1. Search NCBI BLAST databases and WGS/SRA accessions with a consensus KS-domain query.
2. Merge search results and filter loci with clustered KS-domain hits (HSPs).
3. Fetch candidate GenBank and FASTA records.
4. Annotate candidate accessions with antiSMASH.
5. Filter antiSMASH regions and collapse redundancy with MMseqs2.
6. Annotate overlap with known MiBIG clusters.

The commands below can be run manually. For convenience, the repository also includes an optional `Makefile` that wraps the same steps into named targets; see **Optional Makefile shortcuts** below.

<br>

## Optional Makefile shortcuts

The `Makefile` provides step-based shortcuts for running the same workflow commands shown below. Workflow variables can be overridden at runtime.


```bash
make ncbi-download
make ncbi-search
make wgs-bacteria
make wgs-archaea
make wgs-eukaryota
make ks-filter
make fetch NCBI_EMAIL=name@example.org
make metadata
make antismash CONDA_ENV=antismash THREADS_PER_JOB=8 JOBS_AT_ONCE=4
make cleanup
make mibig
```

<br>

## 1. Search NCBI BLAST databases and WGS/SRA accessions

### 1.1 Download and extract NCBI BLAST databases

Skip this step if your BLAST databases already exist elsewhere. By default this writes to `data/ncbi_blastdbs/<db_name>/`.

```bash
bash scripts/download_ncbi_blastdbs.sh \
  ref_prok_rep_genomes:00-24 \
  ref_euk_rep_genomes:000-161 \
  env_nt \
  patnt:00-07 \
  tsa_nt:00-03 \
  nt:000-223
```

Use `name:start-end` only at the download step for sharded NCBI BLAST archive files. Use just `name` for single-archive databases.

### 1.2 Run `tblastn` against the downloaded NCBI databases

```bash
QUERY=queries/KSSignatureConsensusPKSDB.fasta \
DB_ROOT=data/ncbi_blastdbs \
OUT_DIR=outputs/tblastn \
DB_NAMES="ref_prok_rep_genomes ref_euk_rep_genomes env_nt patnt tsa_nt nt" \
bash ncbi_tblastn.sh
```

If BLAST+ is not already on `PATH`, set `BLAST_BIN=/path/to/ncbi-blast/bin`.

### 1.3 Fetch WGS accessions and run `tblastn_vdb`

The full database construction uses both preformatted NCBI BLAST databases and WGS searches. Use `wgs_tblastn.sh` for bacteria, archaea, and eukaryota WGS searches. `prefetch` and `tblastn_vdb` are installed with NCBI SRA Toolkit; `prefetch` downloads the WGS/SRA accessions listed in the `.nvl` alias file, then `tblastn_vdb` searches those accessions. By default, downloaded WGS/SRA accessions go under `data/wgs_<group>/`, and `.tblastn.out` files go under `outputs/tblastn/wgs_<group>/`.

```bash
for group_taxid in bacteria:2 archaea:2157 eukaryota:2759; do
  GROUP=${group_taxid%:*}
  TAXID=${group_taxid#*:}
  CREATE_ALIAS=1 \
  GROUP="$GROUP" \
  TAXID="$TAXID" \
  QUERY=queries/KSSignatureConsensusPKSDB.fasta \
  bash wgs_tblastn.sh
done
```

TaxIDs used in this workflow:

- bacteria: `2`
- archaea: `2157`
- eukaryota: `2759`

For very large `.nvl` files, `wgs_tblastn.sh` splits the alias file internally before `prefetch`. You can also split an alias file manually:

```bash
CHUNK_SIZE=45000 GROUP=bacteria bash scripts/split_nvl.sh bacteria-wgs.nvl
```

To resume an incomplete WGS download/search directory by excluding accessions already present in the current directory:

```bash
bash scripts/make_remaining_nvl.sh bacteria-wgs.nvl bacteria-wgs-remaining.nvl
```

<br>

## 2. Merge search outputs and filter KS clusters

### 2.1 Filter all BLAST/WGS hits for clustered KS domains

At this point, regular `tblastn` outputs and WGS `tblastn_vdb` outputs should both be under `outputs/tblastn/`. This step filters every `.out` file found there and merges candidate loci into one file for the downstream GenBank/FASTA fetch.

```bash
bash scripts/filter_ks_clusters.sh
```

`ks_filter.py` reports subjects with at least 3 merged KS hits, using a 20 kb cluster distance and 3 kb HSP merge distance. This step is intended to enrich for assembly-line PKS-like loci while excluding isolated KS-like hits and many non-assembly-line PKS architectures.

<br>

## 3. Fetch candidate sequences

### 3.1 Fetch GenBank and FASTA records for candidate accessions

Fetch both formats. GenBank records provide metadata and annotations when available, while FASTA records are needed as a fallback for candidates whose GenBank files do not contain full sequence or annotated CDS features.

```bash
NCBI_EMAIL=name@example.org \
python get_genbank.py outputs/ks_filter/clusters_all.txt genbank data/genbank_gb

NCBI_EMAIL=name@example.org \
python get_genbank.py outputs/ks_filter/clusters_all.txt fasta data/genbank_fasta
```

<br>

## 4. Annotate candidate accessions with antiSMASH

### 4.1 Run antiSMASH

```bash
THREADS_PER_JOB=8 JOBS_AT_ONCE=4 CONDA_ENV=antismash \
bash run_antismash.sh data/genbank_gb data/genbank_fasta outputs/antismash antismash
```

If you are using SLURM, submit the same command with your site-specific resource options, for example:

```bash
sbatch --cpus-per-task=32 --mem=128G --time=2-00:00:00 \
  --wrap="THREADS_PER_JOB=8 JOBS_AT_ONCE=4 CONDA_ENV=antismash bash run_antismash.sh data/genbank_gb data/genbank_fasta outputs/antismash antismash"
```

Useful `run_antismash.sh` environment variables:

- `CONDA_ENV`: conda environment to activate, for example `antismash`
- `CONDA_EXE`: full path to `conda` if it is not on `PATH`
- `MODULE_SETUP`: optional module command, for example `module load system parallel/20200822`
- `FORCE_RERUN=1`: rerun all accessions
- `FORCE_IDS=A,B`: rerun selected accessions
- `PY_FETCH_SCRIPT=get_genbank.py`: FASTA fallback fetch helper

<br>

## 5. Filter antiSMASH regions and collapse redundancy

### 5.1 Build accession metadata tables

```bash
python utils.py build_desc data/genbank_gb outputs/accession_desc.txt
python utils.py get_dates data/genbank_gb outputs/accession_dates.txt
```

### 5.2 Clean and cluster antiSMASH regions

This step parses antiSMASH outputs, retains candidate assembly-line PKS regions, generates concatenated protein and PKS-specific FASTA files, collapses redundant entries using MMseqs2, and writes the final summary table.

```bash
python cleanup_clusters.py \
  --desc outputs/accession_desc.txt \
  --dates outputs/accession_dates.txt \
  --antismash_dir outputs/antismash \
  --out_dir outputs/cleaned_clusters \
  --genome_fasta_dir data/genbank_fasta \
  --dedup_min_seq_id 0.90 \
  --dedup_coverage 0.10 \
  --dedup_threads 32 \
  --homolog_min_seq_id 0.80 \
  --homolog_coverage 0.60 \
  --homolog_threads 32
```

Main outputs under `outputs/cleaned_clusters/` include:

- `protein_concat.faa`: concatenated protein sequences for all ORFs in each retained antiSMASH region.
- `pks_concat.faa`: concatenated PKS-relevant proteins/domains for each retained region.
- `mmseqs_results/`: intermediate MMseqs2 clustering outputs.
- `final_summary.tsv`: non-redundant PKSClusterDB summary table.

The staged redundancy strategy is important. The all-ORF concatenation helps identify highly similar or duplicate biosynthetic regions, while the PKS-restricted concatenation prioritizes redundancy removal based on the assembly-line PKS itself rather than unrelated neighboring genes.

<br>

## 6. Annotate known MiBIG clusters

### 6.1 Download MiBIG data and annotate known clusters

Download and extract the MiBIG JSON and GBK archives into `data/mibig/`.

```bash
MIBIG_VERSION=4.0 bash scripts/download_mibig.sh
```

```bash
python find_mibig.py \
  --summary outputs/cleaned_clusters/final_summary.tsv \
  --mibig_gbk_dir data/mibig/mibig_gbk_4.0 \
  --mibig_json_dir data/mibig/mibig_json_4.0 \
  --mibig_table outputs/mibig_clusters.tsv \
  --output_full outputs/final_summary_with_mibig.tsv \
  --output_known outputs/final_summary_known.tsv
```

This step adds known-cluster annotations to the final summary and produces a subset of entries associated with known MiBIG products.

<br>

## Main outputs

The most important final files are:

- `outputs/cleaned_clusters/final_summary.tsv`: primary non-redundant PKSClusterDB table.
- `outputs/final_summary_with_mibig.tsv`: final summary with MiBIG overlap annotations.
- `outputs/final_summary_known.tsv`: subset of PKSClusterDB entries annotated as known clusters through MiBIG matching.
- `outputs/mibig_clusters.tsv`: processed MiBIG reference table used for matching.

Typical columns in `final_summary.tsv` include accession, description, sequence date, species, cluster number, genomic coordinates, cluster length, GC content, domain counts, cluster type, identical-member count, and homolog-member count.

<br>

## Script summary

- `Makefile`: optional workflow wrapper that exposes named targets for the major construction steps.
- `scripts/download_ncbi_blastdbs.sh`: downloads and extracts NCBI BLAST database archives.
- `ncbi_tblastn.sh`: BLAST+ search wrapper for downloaded NCBI BLAST databases.
- `wgs_tblastn.sh`: configurable WGS/SRA `prefetch` and `tblastn_vdb` workflow.
- `scripts/taxid2wgs.pl`: NCBI helper for TaxID-to-WGS alias files.
- `scripts/split_nvl.sh`: splits a WGS `.nvl` alias file into smaller chunks.
- `scripts/make_remaining_nvl.sh`: creates a resume `.nvl` by excluding accessions already present in the working directory.
- `scripts/filter_ks_clusters.sh`: filters all BLAST/WGS result files and writes one merged candidate list.
- `ks_filter.py`: filters BLAST outfmt 6 hits into candidate KS clusters.
- `get_genbank.py`: fetches GenBank or FASTA records from NCBI Entrez.
- `run_antismash.sh`: runs antiSMASH over GenBank inputs, with FASTA fallback.
- `utils.py`: writes accession description and sequence date tables.
- `cleanup_clusters.py`: filters antiSMASH regions, clusters proteins with MMseqs2, and writes final summaries.
- `write_summary.py`: summary-table helper functions used by `cleanup_clusters.py`.
- `scripts/download_mibig.sh`: downloads and extracts MiBIG JSON/GBK archives.
- `find_mibig.py`: adds MiBIG overlap annotations to the final summary.
