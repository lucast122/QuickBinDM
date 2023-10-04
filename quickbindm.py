#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from Bio import SeqIO


class DynamicDatabase:
    def __init__(self, skani_file, ani_threshold=95, query_coverage_threshold=80):
        self.skani_file = skani_file
        self.ani_threshold = ani_threshold
        self.query_coverage_threshold = query_coverage_threshold

    def filter_skani_results(self):
        df = pd.read_csv(self.skani_file, sep='\t')
        filtered_df = df[
            (df['ANI'] >= self.ani_threshold) &
            (df['Align_fraction_query'] >= self.query_coverage_threshold)
            ]
        return filtered_df['Ref_name'].unique()


class DIAMOND:
    def __init__(self, database, query_file, threads=1):
        self.database = database
        self.query_file = query_file
        self.threads = threads

    def align(self):
        pass


class MEGAN:
    def __init__(self, diamond_output):
        self.diamond_output = diamond_output

    def bin_and_correct(self):
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Your pipeline description here.')

    # Mandatory arguments
    parser.add_argument('-i', '--input_contigs', required=True,
                        help='Input contigs after assembly and possibly assembly correction.')
    parser.add_argument('-r', '--reference_seqs', required=True,
                        help='Reference sequences for creating the database. Should be a FASTA or FASTQ file.')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='Path to the output folder that will contain the results.')

    # Optional arguments
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use.'
                                                                     ' Default is 1.')

    args = parser.parse_args()

    # Initialize dynamic database class with skani file and thresholds
    dyn_db = DynamicDatabase(args.reference_seqs)
    remaining_genomes = dyn_db.filter_skani_results()
    print(remaining_genomes)

    # Create DIAMOND database based on remaining_genomes
    # ...

    # Run DIAMOND alignment
    diamond = DIAMOND("path/to/database", args.input_contigs, args.threads)
    diamond.align()

    # Run MEGAN binning and frame-shift correction
    megan = MEGAN("path/to/diamond/output")
    megan.bin_and_correct()

    print("Pipeline completed.")
