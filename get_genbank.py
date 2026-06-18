#!/usr/bin/env python3
import sys
import os
import re
import time
from Bio import Entrez, SeqIO
from http.client import IncompleteRead

def get_genbank_ids(infile):
    ids = set()
    with open(infile) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            genbank = line.split("\t")[0]

            if "|" in genbank:
                parts = genbank.split("|")
                if len(parts) >= 2:
                    genbank_id = parts[1]
                else:
                    genbank_id = genbank
            else:
                genbank_id = genbank

            genbank_id = re.sub(r'^(ref|dbj|emb|gb|tpe)', '', genbank_id)
            genbank_id = genbank_id.replace('|', '')

            ids.add(genbank_id)

    return sorted(ids)


def fetch_batch(batch_ids, seq_format, email, temp_outfile, retries=5, delay=5):
    Entrez.email = email
    for attempt in range(1, retries + 1):
        try:
            with Entrez.efetch(
                db="nucleotide",
                id=",".join(batch_ids),
                rettype=seq_format,
                retmode="text",
                post=True
            ) as handle:
                content = handle.read()

            # If content is bytes, decode
            if isinstance(content, bytes):
                content = content.decode("utf-8")

            with open(temp_outfile, "w") as out_handle:
                out_handle.write(content)

            print(f"Successfully fetched {len(batch_ids)} records to {temp_outfile}")
            return  # success, exit function

        except IncompleteRead as e:
            print(f"IncompleteRead encountered on attempt {attempt}/{retries}. Retrying in {delay}s...")
            time.sleep(delay)

    # If all retries fail
    raise Exception(f"Failed to fetch batch after {retries} attempts due to incomplete read")


def split_sequences(temp_outfile, seq_format, outfolder):
    os.makedirs(outfolder, exist_ok=True)
    count = 0

    in_format = seq_format
    out_format = "gb" if seq_format == "genbank" else seq_format

    with open(temp_outfile) as handle:
        for record in SeqIO.parse(handle, in_format):
            rec_id = record.id

            # Try to extract accession from gi|...|
            m = re.match(r'gi\|(.+?)\|', rec_id)
            if m:
                accession = m.group(1)
            else:
                accession = rec_id

            out_file = os.path.join(outfolder, f"{accession}.{out_format}")
            with open(out_file, "w") as out_handle:
                SeqIO.write(record, out_handle, in_format)
            count += 1

    print(f"Total sequences written: {count}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python get_genbank.py <infile> <format> <outfolder>")
        sys.exit(1)

    infile = sys.argv[1]
    seq_format = sys.argv[2].lower()
    outfolder = sys.argv[3]

    if seq_format not in ["genbank", "fasta", "gb"]:
        print(f"Unsupported format: {seq_format}. Use 'genbank', 'fasta', or 'gb'.")
        sys.exit(1)

    email = os.environ.get("NCBI_EMAIL", "your.email@example.com")
    if email == "your.email@example.com":
        print("Set NCBI_EMAIL for Entrez requests, e.g. NCBI_EMAIL=name@example.org", file=sys.stderr)
    ids = get_genbank_ids(infile)
    print(f"Num of unique IDs: {len(ids)}")

    batch_size = 100  # safe size for long requests
    out_format = "gb" if seq_format == "genbank" else seq_format

    for i in range(0, len(ids), batch_size):
        batch = []
        for genbank_id in ids[i:i+batch_size]:
            expected_file = os.path.join(outfolder, f"{genbank_id}.{out_format}")
            if not os.path.exists(expected_file):
                batch.append(genbank_id)
        if not batch:
            continue

        temp_outfile = f"temp_batch_{i//batch_size+1}.{seq_format}"
        print(f"Fetching batch {i // batch_size + 1} with {len(batch)} IDs...")
        fetch_batch(batch, seq_format, email, temp_outfile, retries=3, delay=10)
        split_sequences(temp_outfile, seq_format, outfolder)
        os.remove(temp_outfile)
        time.sleep(2.5)  # polite pause

    print(f"All done! Split files in: {outfolder}")
