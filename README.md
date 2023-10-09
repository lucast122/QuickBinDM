 # QuickBinDM
Quick binning and frame-shift correction of long read metagenomic assemblies using DIAMOND+MEGAN pipeline with dynamically created reference DB

## Overview

QuickBinDM is a bioinformatics pipeline designed to efficiently perform metagenomic binning of contigs generated from nanopore sequencing. The pipeline integrates DIAMOND, MEGAN, and Skani to create a dynamic DIAMOND database, which significantly speeds up the binning process.

## Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/your_username/QuickBinDM.git
```

## Requirements

- Python 3.x
- DIAMOND
- MEGAN
- Skani

Please make sure to install the dependencies before running the pipeline.

## Usage

Navigate to the directory containing `quickbindm.py` and execute the following command:

```bash
python quickbindm.py -i <input_contigs> -r <reference_seqs> -o <output_folder> [options]
```

### Parameters

#### Mandatory Parameters

- `-i, --input_contigs`: The path to the input contigs after assembly and possibly assembly correction. Should be a FASTA or FASTQ file.
- `-r, --reference_seqs`: The path to the reference sequences for creating the database. Should be a FASTA or FASTQ file.
- `-o, --output_folder`: The path to the output folder that will contain the results.

#### Optional Parameters

- `-t, --threads`: The number of threads to use. Default is 1.
- `-a, --ani`: The ANI (Average Nucleotide Identity) cutoff. Default is 95.
- `-q, --querycoverage`: The minimum query coverage for ANI calculation. Default is 80.
- `-c, --refcoverage`: The minimum reference coverage for ANI calculation. Default is 80.

## Example

```bash
python quickbindm.py -i input_contigs.fasta -r reference_seqs.fasta -o output/ -t 4 -a 95 -q 80 -c 80
```

## Contributors

- Your Name (your@email.com)

## License

This project is licensed under the MIT License.
