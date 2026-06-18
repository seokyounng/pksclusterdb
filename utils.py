#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import SeqIO
import re

def build_accession_desc(infile, outfile):
    """
    If infile is a GenBank file or folder of GenBank files:
       extract accession and better description (DEFINITION + ORGANISM).
    Else:
       read tab-delimited file with accession and description columns.

    Writes accession and description to output file.
    """

    if Path(infile).is_dir():
        # Process GenBank files in folder
        genbank_path = Path(infile)
        with open(outfile, "w") as fout:
            for gb_file in genbank_path.glob("*.gb"):
                try:
                    for record in SeqIO.parse(str(gb_file), "genbank"):
                        accession = record.id
                        desc = record.description.strip() if record.description else ""

                        # Append organism info if not present
                        organism = record.annotations.get("organism", "").strip()
                        if organism and organism not in desc:
                            desc += f" [Organism: {organism}]"

                        fout.write(f"{accession}\t{desc}\n")
                except Exception as e:
                    print(f"Error processing {gb_file}: {e}", file=sys.stderr)
    else:
        # Assume tab-delimited file
        with open(infile) as fin, open(outfile, "w") as fout:
            for line in fin:
                parts = line.strip().split("\t")
                if not parts:
                    continue
                accession = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                fout.write(f"{accession}\t{desc}\n")


def get_sequence_dates(genbank_folder, outfile=None):
    """
    Reads all *.gb files in genbank_folder and writes accession/date pairs.
    Tries Biopython's record.annotations['date'] first.
    Falls back to parsing LOCUS line manually if needed.
    """
    genbank_path = Path(genbank_folder)
    out_handle = open(outfile, "w") if outfile else sys.stdout
    try:
        for gb_file in genbank_path.glob("*.gb"):
            try:
                record = next(SeqIO.parse(str(gb_file), "genbank"))
                accession = record.id

                # Try Biopython first
                date = record.annotations.get("date", "").strip()

                # If date is missing or empty, fallback to LOCUS line
                if not date:
                    with open(gb_file) as f:
                        for line in f:
                            if line.startswith("LOCUS"):
                                match = re.search(r"\d{2}-[A-Z]{3}-\d{4}", line)
                                if match:
                                    date = match.group(0)
                                break

                print(f"{accession}\t{date}", file=out_handle)
            except Exception as e:
                print(f"Error processing {gb_file}: {e}", file=sys.stderr)
    finally:
        if outfile:
            out_handle.close()


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Build accession desc: python utils.py build_desc <input_cluster_file> <output_desc_file>")
        print("  Get sequence dates:  python utils.py get_dates <genbank_folder> <output_dates_file>")
        sys.exit(1)

    command = sys.argv[1]

    if command == "build_desc":
        if len(sys.argv) != 4:
            print("Usage: python utils.py build_desc <input_cluster_file> <output_desc_file>")
            sys.exit(1)
        build_accession_desc(sys.argv[2], sys.argv[3])

    elif command == "get_dates":
        if len(sys.argv) not in (3, 4):
            print("Usage: python utils.py get_dates <genbank_folder> [output_dates_file]")
            sys.exit(1)
        outfile = sys.argv[3] if len(sys.argv) == 4 else None
        get_sequence_dates(sys.argv[2], outfile)

    else:
        print(f"Unknown command '{command}'")
        sys.exit(1)


if __name__ == "__main__":
    main()
