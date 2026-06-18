.PHONY: help dirs ncbi-download ncbi-search wgs-bacteria wgs-archaea wgs-eukaryota \
	ks-filter fetch-genbank fetch-fasta fetch metadata antismash cleanup \
	mibig-download mibig-annotate mibig

PYTHON ?= python
QUERY ?= queries/KSSignatureConsensusPKSDB.fasta

DB_ROOT ?= data/ncbi_blastdbs
DB_SPECS ?= ref_prok_rep_genomes:00-24 ref_euk_rep_genomes:000-161 env_nt patnt:00-07 tsa_nt:00-03 nt:000-223
DB_NAMES ?= ref_prok_rep_genomes ref_euk_rep_genomes env_nt patnt tsa_nt nt

TBLASTN_DIR ?= outputs/tblastn
KS_FILTER_OUT ?= outputs/ks_filter
CLUSTERS ?= $(KS_FILTER_OUT)/clusters_all.txt

GB_DIR ?= data/genbank_gb
FASTA_DIR ?= data/genbank_fasta
NCBI_EMAIL ?=

ANTISMASH_OUT ?= outputs/antismash
ANTISMASH_RUN ?= antismash
THREADS_PER_JOB ?= 8
JOBS_AT_ONCE ?= 4
CONDA_ENV ?=

CLEAN_OUT ?= outputs/cleaned_clusters
DEDUP_MIN_SEQ_ID ?= 0.90
DEDUP_COVERAGE ?= 0.10
DEDUP_THREADS ?= 32
HOMOLOG_MIN_SEQ_ID ?= 0.80
HOMOLOG_COVERAGE ?= 0.60
HOMOLOG_THREADS ?= 32

MIBIG_VERSION ?= 4.0
MIBIG_DIR ?= data/mibig

help:
	@echo "PKSClusterDB workflow targets"
	@echo ""
	@echo "  make ncbi-download       Download/extract NCBI BLAST databases"
	@echo "  make ncbi-search         Run tblastn against downloaded NCBI databases"
	@echo "  make wgs-bacteria        Fetch/search bacterial WGS accessions"
	@echo "  make wgs-archaea         Fetch/search archaeal WGS accessions"
	@echo "  make wgs-eukaryota       Fetch/search eukaryotic WGS accessions"
	@echo "  make ks-filter           Merge BLAST/WGS hits and filter KS clusters"
	@echo "  make fetch               Fetch GenBank and FASTA records; set NCBI_EMAIL=..."
	@echo "  make metadata            Build accession description/date tables"
	@echo "  make antismash           Run antiSMASH"
	@echo "  make cleanup             Cluster antiSMASH regions and write final_summary.tsv"
	@echo "  make mibig               Download MiBIG and annotate known clusters"
	@echo ""
	@echo "Override variables as needed, for example:"
	@echo "  make fetch NCBI_EMAIL=name@example.org"
	@echo "  sbatch --wrap=\"make antismash CONDA_ENV=antismash THREADS_PER_JOB=8 JOBS_AT_ONCE=4\""

dirs:
	mkdir -p data outputs log

ncbi-download: dirs
	DB_ROOT="$(DB_ROOT)" bash scripts/download_ncbi_blastdbs.sh $(DB_SPECS)

ncbi-search: dirs
	QUERY="$(QUERY)" DB_ROOT="$(DB_ROOT)" OUT_DIR="$(TBLASTN_DIR)" DB_NAMES="$(DB_NAMES)" bash ncbi_tblastn.sh

wgs-bacteria: dirs
	CREATE_ALIAS=1 GROUP=bacteria TAXID=2 QUERY="$(QUERY)" bash wgs_tblastn.sh

wgs-archaea: dirs
	CREATE_ALIAS=1 GROUP=archaea TAXID=2157 QUERY="$(QUERY)" bash wgs_tblastn.sh

wgs-eukaryota: dirs
	CREATE_ALIAS=1 GROUP=eukaryota TAXID=2759 QUERY="$(QUERY)" bash wgs_tblastn.sh

ks-filter: dirs
	TBLASTN_DIR="$(TBLASTN_DIR)" OUT_DIR="$(KS_FILTER_OUT)" bash scripts/filter_ks_clusters.sh

fetch: fetch-genbank fetch-fasta

fetch-genbank: dirs
	test -n "$(NCBI_EMAIL)" || { echo "Set NCBI_EMAIL=name@example.org"; exit 2; }
	NCBI_EMAIL="$(NCBI_EMAIL)" $(PYTHON) get_genbank.py "$(CLUSTERS)" genbank "$(GB_DIR)"

fetch-fasta: dirs
	test -n "$(NCBI_EMAIL)" || { echo "Set NCBI_EMAIL=name@example.org"; exit 2; }
	NCBI_EMAIL="$(NCBI_EMAIL)" $(PYTHON) get_genbank.py "$(CLUSTERS)" fasta "$(FASTA_DIR)"

metadata: dirs
	$(PYTHON) utils.py build_desc "$(GB_DIR)" outputs/accession_desc.txt
	$(PYTHON) utils.py get_dates "$(GB_DIR)" outputs/accession_dates.txt

antismash: dirs
	THREADS_PER_JOB="$(THREADS_PER_JOB)" JOBS_AT_ONCE="$(JOBS_AT_ONCE)" CONDA_ENV="$(CONDA_ENV)" \
		bash run_antismash.sh "$(GB_DIR)" "$(FASTA_DIR)" "$(ANTISMASH_OUT)" "$(ANTISMASH_RUN)"

cleanup: dirs
	$(PYTHON) cleanup_clusters.py \
		--desc outputs/accession_desc.txt \
		--dates outputs/accession_dates.txt \
		--antismash_dir "$(ANTISMASH_OUT)" \
		--out_dir "$(CLEAN_OUT)" \
		--genome_fasta_dir "$(FASTA_DIR)" \
		--dedup_min_seq_id "$(DEDUP_MIN_SEQ_ID)" \
		--dedup_coverage "$(DEDUP_COVERAGE)" \
		--dedup_threads "$(DEDUP_THREADS)" \
		--homolog_min_seq_id "$(HOMOLOG_MIN_SEQ_ID)" \
		--homolog_coverage "$(HOMOLOG_COVERAGE)" \
		--homolog_threads "$(HOMOLOG_THREADS)"

mibig: mibig-download mibig-annotate

mibig-download: dirs
	MIBIG_VERSION="$(MIBIG_VERSION)" MIBIG_DIR="$(MIBIG_DIR)" bash scripts/download_mibig.sh

mibig-annotate: dirs
	$(PYTHON) find_mibig.py \
		--summary "$(CLEAN_OUT)/final_summary.tsv" \
		--mibig_gbk_dir "$(MIBIG_DIR)/mibig_gbk_$(MIBIG_VERSION)" \
		--mibig_json_dir "$(MIBIG_DIR)/mibig_json_$(MIBIG_VERSION)" \
		--mibig_table outputs/mibig_clusters.tsv \
		--output_full outputs/final_summary_with_mibig.tsv \
		--output_known outputs/final_summary_known.tsv
