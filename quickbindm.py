#!/usr/bin/env python3
import glob
import os
import argparse
import pandas as pd
import zipfile
import requests
import subprocess
import time
from Bio import Entrez, SeqIO
from concurrent.futures import ThreadPoolExecutor




class DynamicDatabase:
    def __init__(self, contigs_path, reference_seq_path , ani_threshold, query_coverage_threshold, ref_coverage_threshold, output_folder,threads):
        self.contigs_path = contigs_path
        self.reference_seq_path = reference_seq_path
        self.ani_threshold = ani_threshold
        self.query_coverage_threshold = query_coverage_threshold
        self.ref_coverage_threshold = ref_coverage_threshold
        self.output_folder = output_folder
        self.threads = threads
        self.refseq_ids = []
        # Initialize pipeline
        self.initialize_pipeline()

    def initialize_pipeline(self):
        # Filter remaining genomes

        print(f"Calculating ANI between contigs and reference genomes...")

        if os.path.exists(f"{self.output_folder}/contigs_vs_refseq.txt"):
            print("Found existing skani result, continue using existing result...")
        else:
            self.calculate_ani()


        remaining_genomes = self.filter_skani_results()
        print(f"Numbers of genomes in dynamic DB: {len(remaining_genomes)}")


        start_time = time.time()

        remaining_genome_ids = [genome.split(" ")[0] for genome in remaining_genomes]

        print(f"Fetching {len(remaining_genome_ids)} refseq ids...")

        start_time = time.time()
        refseq_ids = self.get_refseq_assembly_ids(remaining_genome_ids,num_threads=10,batch_size=100)
        self.refseq_ids = refseq_ids
        print(refseq_ids)
        elapsed_time_1 = time.time() - start_time

        print(f"Time for fetching refseq ids:  {elapsed_time_1}")
        print(f"Bulk downloading genomes using 10 threads to {self.output_folder}...")

        start_time = time.time()
        self.download_proteins_for_refseq_assembly_ids(refseq_ids, self.output_folder, num_threads=10)


        elapsed_time_2 = time.time() - start_time
        print(f"Time for fetching proteins: {elapsed_time_2}")





        # Extract files
        self.extract_files_in_directory(self.output_folder)

        protein_sequences_path = f"{self.output_folder}/proteins/concatenated.fasta"

        # Concatenate fasta
        self.concatenate_fasta(f"{self.output_folder}/proteins/ncbi_dataset/data/", protein_sequences_path)

        # Create Diamond DB
        self.create_diamond_db(protein_sequences_path, f"{self.output_folder}/dynamic_db.dmnd")

    def calculate_ani(self):

        command = ['skani', 'dist', '--qi', '-q', self.contigs_path, '--ri', '-r', self.reference_seq_path, '-t', '64']

        with open(f"{self.output_folder}/contigs_vs_refseq.txt", "w") as outfile:
            subprocess.run(command, stdout=outfile, check=True)



    def filter_skani_results(self):
        df = pd.read_csv(f"{self.output_folder}/contigs_vs_refseq.txt", sep='\t')
        filtered_df = df[
            (df['ANI'] >= self.ani_threshold) &
            (df['Align_fraction_query'] >= self.query_coverage_threshold) &
            (df['Align_fraction_ref'] >= self.ref_coverage_threshold)


            ]
        return filtered_df['Ref_name'].unique()

    def fetch_batch(self,ids_batch, results):
        ids_str = ','.join(ids_batch)
        handle = Entrez.efetch(db="nucleotide", id=ids_str, rettype="gb", retmode="text")
        genome_records = SeqIO.parse(handle, "genbank")
        for record in genome_records:
            results.append(record.dbxrefs[-1].split("Assembly:")[-1])
        handle.close()

    def extract_files_in_directory(self,directory_path):
        zip_files = glob.glob(f"{directory_path}/*.zip")
        print(f"Extracting {len(zip_files)} files...")
        for file in zip_files:
            # print(f"Extracting file {file}")

            try:
                with zipfile.ZipFile(file, 'r') as zip_ref:
                    zip_ref.extractall(f"{directory_path}/proteins/")
            except zipfile.BadZipfile as e:
                print(f"File {file} is not a normal zip file")
                continue

    def concatenate_fasta(self,input_path, out_file_path):
        with open(out_file_path, 'w') as outfile:
            for fname in glob.glob(f"{input_path}/*/*.faa"):
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def download_protein(self,refseq_assembly_id, output_folder):
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{refseq_assembly_id}/download"
        params = {
            "include_annotation_type": "PROT_FASTA",
            "filename": f"{refseq_assembly_id}.zip"
        }
        headers = {
            "Accept": "application/zip"
        }

        response = requests.get(url, params=params, headers=headers)

        if response.status_code == 200:
            with open(f"{output_folder}/{refseq_assembly_id}.zip", "wb") as file:
                file.write(response.content)
        else:
            print(f"Failed to retrieve data for {refseq_assembly_id}: {response.status_code}")

    def download_proteins_for_refseq_assembly_ids(self,refseq_assembly_ids, output_folder, num_threads):
        # Ensure the output directory exists
        os.makedirs(output_folder, exist_ok=True)

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            executor.map(lambda refseq_assembly_id: self.download_protein(refseq_assembly_id, output_folder),
                         refseq_assembly_ids)

    def download_proteins_for_refseq_assembly_id(self,refseq_assembly_id, output_folder):
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{refseq_assembly_id}/download"
        params = {
            "include_annotation_type": "PROT_FASTA",
            f"filename": f"{output_folder}{refseq_assembly_id}.zip"
        }
        headers = {
            "Accept": "application/zip"
        }

        response = requests.get(url, params=params, headers=headers)

        # Check for a valid response
        if response.status_code == 200:
            # Write the content to a file
            with open(f"{output_folder}/{refseq_assembly_id}.zip", "wb") as file:
                file.write(response.content)
        else:
            print(f"Failed to retrieve data: {response.status_code}")

    def get_refseq_assembly_id(self, refseq_genome_id):
        # Fetch genome record
        handle = Entrez.efetch(db="nucleotide", id=refseq_genome_id, rettype="gb", retmode="text")
        genome_record = SeqIO.read(handle, "genbank")
        handle.close()
        # print(dir(genome_record))
        refseq_assembly_id = genome_record.dbxrefs[-1].split("Assembly:")[-1]
        return refseq_assembly_id

    def get_refseq_assembly_ids(self, refseq_genome_ids, batch_size=200, num_threads=5):
        total_ids = len(refseq_genome_ids)
        results = []

        # Split ids into batches
        batches = [refseq_genome_ids[i:i + batch_size] for i in range(0, total_ids, batch_size)]

        # Create a ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            # Use the executor.map method to execute fetch_batch on each batch of ids
            list(executor.map(lambda batch: self.fetch_batch(batch, results), batches))
        # print(results)
        return results

    # wrapper function
    def download_protein_for_refseq_genome_id(self,refseq_genome_id, output_folder):
        download_proteins_for_refseq_assembly_id(get_refseq_assembly_id(refseq_genome_id), output_folder)

    def create_diamond_db(self, input_fasta, db_name):
        command = ['diamond', 'makedb', '--in', input_fasta, '-d', db_name]
        subprocess.run(command, check=True)


class DIAMOND:
    def __init__(self, database, query_file, out, threads):
        self.database = database
        self.query_file = query_file
        self.out = out
        self.threads = threads

    def align(self):
        command = ['diamond', 'blastx', '-q', self.query_file, '-d', self.database, '-p', f"{self.threads}", "-o",
                   f"{self.out}/contigs_vs_dynamic_db.daa", "-f100",  "--long-reads", '-t', '/dev/shm/']
        subprocess.run(command, check=True)


class MEGAN:
    def __init__(self, output_folder, daa_file, mapping_file,threads):
        self.output_folder = output_folder
        self.daa_file = daa_file
        self.mapping_file = mapping_file
        self.threads = threads


    def meganize(self):
        command = ['/home/minion-computer/megan/tools/daa-meganizer', '-i', self.daa_file, '-mdb', self.mapping_file, '-lg', '-t', f"{self.threads}"]
        subprocess.run(command, check=True)


    def bin_and_correct(self):
        command = ['/home/minion-computer/megan/tools/read-extractor', '-i', self.daa_file, '-c', 'Taxonomy', '-o', f'{self.output_folder}/%t.fasta', '-fsc']
        subprocess.run(command, check=True)

class CHECKM:
    def __init__(self,output_folder, num_threads):
        self.num_threads = num_threads
        self.output_folder = output_folder

    def run_checkm(self):
        command = ['checkm', 'lineage_wf', self.output_folder, self.output_folder, '-t', str(self.num_threads), '-x',
                   'fasta', '--pplacer_threads', '64']

        with open(f"{self.output_folder}/checkm.txt", "w") as outfile:
            subprocess.run(command, stdout=outfile, check=True)


if __name__ == "__main__":
    Entrez.email = "timolucas1994@gmail.com"
    print("Parsing arguments...")
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

    parser.add_argument('-a', '--ani', type=int, default=95, help='ANI cutoff.'
                                                                  ' Default is 95.')

    parser.add_argument('-q', '--querycoverage', type=int, default=80,
                        help='Minimum query coverage for ANI calculation.'
                             ' Default is 80')

    parser.add_argument('-c', '--refcoverage', type=int, default=80,
                        help='Minimum reference coverage for ANI calculation.'
                             ' Default is 80')

    args = parser.parse_args()

    # Create output folder if it doesn't exist
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # Initialize dynamic database class with skani file and thresholds
    print("Creating dynamic DB...")

    dyn_db = DynamicDatabase(args.input_contigs, args.reference_seqs, args.ani, args.querycoverage, args.refcoverage, args.output_folder, args.threads)
    # dyn_db.extract_remaining_genomes(args.reference_seqs, args.output_folder)
    # dyn_db.extract_remaining_genomes_parallel(args.reference_seqs, args.output_folder, n_threads=args.threads)
    print("Filtering references using skani ANI results...")






    # Run DIAMOND alignment
    diamond = DIAMOND(f"{args.output_folder}/dynamic_db.dmnd", args.input_contigs, args.output_folder, args.threads)
    diamond.align()

    # Run MEGAN binning and frame-shift correction
    megan = MEGAN(args.output_folder,f"{args.output_folder}/contigs_vs_dynamic_db.daa", "megan-map-Feb2022-ue.db", args.threads)
    megan.meganize()
    megan.bin_and_correct()

    checkm = CHECKM(args.output_folder, args.threads)
    checkm.run_checkm()



    print("Pipeline completed.")
