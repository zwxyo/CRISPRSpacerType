import re
import os
import sys
import glob
import pandas as pd
import subprocess
from utils import create_folder


def BLAST_M(input_dir, user_db, output_dir):

    # parser = argparse.ArgumentParser(description="Run BLASTn against a custom database")
    # parser.add_argument("--bi", type=str, required=True, help="Path to the folder containing query .fna files")
    # args = parser.parse_args()

    # output_dir = Path("blast_results")
    # output_dir.mkdir(parents=True, exist_ok=True)
    create_folder(output_dir)

    # user_path = Path("user_gene_seq/")
    user_path = f"{user_db}"

    # if user_path.is_dir() and any(user_path.iterdir()):
    #     db_folder = user_path
    if os.path.isdir(user_path) and any(os.scandir(user_path)):
        db_folder = user_path

    # db_files = list(db_folder.glob("*.fna"))
    db_files = glob.glob(os.path.join(db_folder, "*.fna"))
    if not db_files:
        raise FileNotFoundError("User database file not found (.fna)")

    # merged_fna_file = Path("user_merged_sequences.fna")
    # merged_fna_file = f"{user_db}/user_merged_sequences.fna"
    # with merged_fna_file.open("w", encoding="utf-8") as merged_file:
    #     for file in db_files:
    #         file_label = file.stem
    #         with file.open("r", encoding="utf-8") as f:
    #             for line in f:
    #                 if line.startswith(">"):
    #                     line = f"{line.strip()} | source={file_label}\n"
    #                 merged_file.write(line)
    merged_fna_file = os.path.join(user_db, "user_merged_sequences.fna")
    with open(merged_fna_file, "w", encoding="utf-8") as merged_file:
        for file in db_files:
            file_label = os.path.splitext(os.path.basename(file))[0]
            with open(file, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith(">"):
                        line = f"{line.strip()} | source={file_label}\n"
                    merged_file.write(line)

    # db_name = Path("userdb/my_custom_blastdb")
    # db_name = f"{user_db}/my_custom_blastdb"
    db_name = os.path.join(user_db, "ct_blastdb")
    subprocess.run([
        "makeblastdb",
        "-in", str(merged_fna_file),
        "-dbtype", "nucl",
        "-out", str(db_name)
    ])

    # if not any((db_name.parent / f"{db_name.name}{ext}").exists() for ext in [".nsq", ".nin", ".nhr"]):
    #     raise RuntimeError("BLAST database creation failed, please check the input file!")
    db_files = [f"{db_name}{ext}" for ext in [".nsq", ".nin", ".nhr"]]
    if not any(os.path.exists(db_file) for db_file in db_files):
        raise RuntimeError("BLAST database creation failed, please check the input file!")

    # query_folder = Path("E:/jupyter-notebook/research_topic/query_seq/")
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
