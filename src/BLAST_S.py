import re
import os
import sys
import glob
import pandas as pd
from pathlib import Path
import subprocess
from utils import create_folder

def BLAST_S(input_dir, output_dir):
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parent_dir = os.path.dirname(script_dir)

    database_dir = os.path.join(parent_dir, "database/BLAST")
    # os.chdir('..')
    # database_dir = os.getcwd()

    # parser = argparse.ArgumentParser(description="Run BLASTn against a custom database")
    # parser.add_argument("--bi", type=str, required=True, help="Path to the folder containing query .fna files")
    # args = parser.parse_args()

    # output_dir = Path("blast_results")
    # output_dir.mkdir(parents=True, exist_ok=True)
    create_folder(output_dir)

    merged_fna_file = f"{database_dir}/merged_sequences.fna"
    print("merged_fna_file",merged_fna_file)

    db_name = f"{database_dir}/db/ct_blastdb"
    print("db_name",db_name)

    # query_folder = Path("query_seq/")
    # query_folder = Path(args.query_folder)
    query_folder = f"{input_dir}"
    # if not query_folder.is_dir():
    #     raise NotADirectoryError(f"The query file directory {query_folder} does not exist or is invalid")
    #
    # query_files = list(query_folder.glob("*.fna"))
    # if not query_files:
    #     raise FileNotFoundError("Query file not found (.fna)")
    if not os.path.isdir(query_folder):
        raise NotADirectoryError(f"The query file directory {query_folder} does not exist or is invalid")

    query_files = glob.glob(os.path.join(query_folder, "*.fna"))
    if not query_files:
        raise FileNotFoundError("Query file not found (.fna)")

    blast_columns = [
        "query_id", "subject_id", "identity", "alignment_length",
        "mismatches", "gap_opens", "query_start", "query_end",
        "subject_start", "subject_end", "evalue", "bit_score"
    ]

    # with merged_fna_file.open("r", encoding="utf-8") as f:
    #     fna_lines = [line.strip() for line in f if line.startswith(">")]

    with open(merged_fna_file, "r", encoding="utf-8") as f:
        fna_lines = [line.strip() for line in f if line.startswith(">")]

    source_dict = {
        re.search(r">(\S+)", line).group(1): re.search(r"source=([\w.-]+)", line).group(1)
        for line in fna_lines if re.search(r"source=([\w.-]+)", line)
    }

    # BLASTn
    for query_file in query_files:
        # query_prefix = re.match(r"([^_]*_[^_]*)_", query_file.stem).group(1)
        query_prefix = re.search(r"([^\\\/]+)\.fna$",query_file).group(1)

        output_file = f"{output_dir}/blast_results_{query_prefix}.txt"

        subprocess.run([
            "blastn",
            "-query", str(query_file),
            "-db", str(db_name),
            "-out", str(output_file),
            "-outfmt", "6",
            "-evalue", "1e-5",
            "-max_target_seqs", "20"
        ])

        blast_results = pd.read_csv(output_file, sep="\t", names=blast_columns)

        blast_results["source"] = blast_results["subject_id"].map(source_dict)

        blast_results.to_csv(output_file, sep="\t", index=False)
        print(blast_results)
        print('\n')
        print("--------------------------------------------------------------------------------")
