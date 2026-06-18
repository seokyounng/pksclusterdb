#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path

# Distance thresholds
DIST_CUTOFF = 20000  # cluster nearby KSs if <= 20kb apart
TOO_CLOSE = 3000     # merge HSPs if <3kb apart

def main():
    infile = sys.argv[1]
    if Path(infile).stat().st_size == 0:
        return

    # Load BLAST tabular: outfmt 6
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    try:
        df = pd.read_csv(infile, sep="\t", names=cols)
    except pd.errors.EmptyDataError:
        return
    if df.empty:
        return

    # Keep subject & hit coords
    df['start'] = df[['sstart', 'send']].min(axis=1)
    df['end']   = df[['sstart', 'send']].max(axis=1)

    # For each subject (contig/scaffold), cluster HSPs by distance
    for sseqid, subdf in df.groupby('sseqid'):
        subdf = subdf[['start', 'end']].sort_values('start').reset_index(drop=True)

        clusters = []
        cluster = [ subdf.loc[0].tolist() ]
        last_end = subdf.loc[0, 'end']

        for idx in range(1, len(subdf)):
            row = subdf.loc[idx]
            distance = row['start'] - last_end

            if distance <= DIST_CUTOFF:
                cluster.append(row.tolist())
                last_end = max(last_end, row['end'])
            else:
                clusters.append(cluster)
                cluster = [ row.tolist() ]
                last_end = row['end']

        clusters.append(cluster)

        # Merge HSPs within each cluster that are <3kb apart
        for cluster_idx, cluster in enumerate(clusters):
            cluster.sort()
            merged = []
            current = list(cluster[0])

            for hsp in cluster[1:]:
                distance = hsp[0] - current[1]
                if distance < TOO_CLOSE:
                    current[1] = max(current[1], hsp[1])
                else:
                    merged.append(tuple(current))
                    current = list(hsp)
            merged.append(tuple(current))

            if len(merged) >= 3:
                print(f"{sseqid}\tCluster{cluster_idx+1}\tNum_KS:{len(merged)}\tCoords:{merged}")

if __name__ == '__main__':
    main()
