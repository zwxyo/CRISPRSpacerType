import os
import pandas as pd
import csv
import re
import glob
from pathlib import Path
from collections import defaultdict

# search for eligible GFF files
def get_gff_file(separated_proteins_dir, hit_id):

    hit_id_parts = hit_id.split('_')

    gff_file_prefix = '_'.join(hit_id_parts[:2])
    # print(gff_file_prefix)

    gff_pattern = os.path.join(separated_proteins_dir, f"*{gff_file_prefix}*.gff")
    # print(gff_pattern)

    matched_files = glob.glob(gff_pattern)

    if matched_files:
        return matched_files[0]
    else:
        print(f"Warning: no gff file found for {hit_id}.")
        return None


# read the TSV file
def read_best_solution_tsv(tsv_file):
    try:
        df = pd.read_csv(tsv_file, sep='\t', comment='#')
        print("cas system exist")
        return df
    except pd.errors.EmptyDataError:
        #print(f"Warning: {tsv_file} is empty or contains only comments.")
        return pd.DataFrame()

# determine the system type
def classify_text(text):
    contains_ie = bool(re.search(r'I-E', text))
    contains_if = bool(re.search(r'I-F', text))

    if contains_ie and contains_if:
        return "Contains both I-E and I-F"
    elif contains_ie:
        return "Contains only I-E"
    elif contains_if:
        return "Contains only I-F"
    else:
        return "Contains neither I-E nor I-F"

def classify_file(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', comment='#')
    text = ' '.join(df.astype(str).values.flatten())
    return classify_text(text)

# obtain the loci from the GFF file
def extract_coordinates(gff_file, target_id):
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                if len(columns) > 8:
                    attributes = columns[8]
                    match = re.search(r'ID=' + re.escape(target_id) + r'\b', attributes)
                    if match:
                        start, end = columns[3], columns[4]
                        return int(start), int(end)
    return None

# process best_solution.tsv
def process_best_solution(folder):
    analyzed_folder = os.path.join(folder, "cas").replace("\\", "/")

    for subfolder in os.listdir(analyzed_folder):
        subdir = os.path.join(analyzed_folder, subfolder).replace("\\", "/")

        if not os.path.isdir(subdir):
            continue

        print("Processing subdirectory:", subdir)
        all_results = []
        num_systems = 1

        for sub_subdir, dirs, files in os.walk(subdir):
            for file in files:
                if file == "best_solution.tsv":
                    tsv_path = os.path.join(sub_subdir, file).replace("\\", "/")
                    print(f"Processing {tsv_path}...\n")

                    df = read_best_solution_tsv(tsv_path)

                    if df.empty:
                        print(f"  No valid data found in {tsv_path}. Skipping...\n")
                        continue

                    replicon = df["replicon"].iloc[0]
                    path = Path(sub_subdir)
                    target_dir = path.parent.parent
                    separated_proteins_dir = os.path.join(target_dir, "separated_proteins")
                    output_csv = os.path.join(target_dir, "cas_info.csv")

                    text = classify_file(tsv_path)
                    classify = classify_text(text)
                    system_coordinates = defaultdict(list)

                    for _, row in df.iterrows():
                        hit_id = row["hit_id"]
                        model_fqn = row["model_fqn"]
                        id = row["hit_pos"]

                        gff_file = os.path.join(separated_proteins_dir, f"*{replicon}*.gff")
                        matched_files = glob.glob(gff_file)
                        if not matched_files or not os.path.exists(matched_files[0]):
                            print(f"  GFF file not found for pattern {gff_file}")
                            continue

                        if classify == "Contains only I-E":
                            coordinates = extract_coordinates(matched_files[0], f"1_{id}")
                            system_coordinates["I-E"].append(coordinates)
                        elif classify == "Contains only I-F":
                            coordinates = extract_coordinates(matched_files[0], f"1_{id}")
                            system_coordinates["I-F"].append(coordinates)
                        elif classify == "Contains both I-E and I-F":
                            if "I-E" in model_fqn:
                                coordinates = extract_coordinates(matched_files[0], f"1_{id}")
                                system_coordinates["I-E"].append(coordinates)
                            elif "I-F" in model_fqn:
                                coordinates = extract_coordinates(matched_files[0], f"1_{id}")
                                system_coordinates["I-F"].append(coordinates)


                    for system, coords in system_coordinates.items():
                        all_starts = [start for start, end in coords]
                        all_ends = [end for start, end in coords]
                        final_start = min(all_starts)
                        final_end = max(all_ends)
                        print(f"Final coordinates for {system}: {final_start} - {final_end}")
                        gene_names = df[df["model_fqn"].str.contains(system, na=False)]["gene_name"].tolist()
                        all_results.append([num_systems, system, final_start, final_end, ', '.join(gene_names)])
                        num_systems += 1
                        print("Results:", all_results, '\n')

            if all_results:
                with open(output_csv, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['loci_ID', 'loci_type', 'loci_start', 'loci_end', 'each_gene'])
                    writer.writerows(all_results)

        print("--------------------------------------------------------------------------------")


def filter_complete_cas_systems(base_folder):
    required_genes = {
        "I-E": ["cas5", "cas7", "cse2", "cas8e", "cas3"],
        "I-F": ["cas7", "cas5", "cas8f", "cas3"]
    }

    for root, dirs, files in os.walk(base_folder):
        if "cas_info.csv" in files:
            cas_info_path = os.path.join(root, "cas_info.csv")
            df = pd.read_csv(cas_info_path)

            filtered_rows = []
            for _, row in df.iterrows():
                loci_type = row.get("loci_type", "")
                each_gene = row.get("each_gene", "")
                each_gene_list = str(each_gene).split(", ") if pd.notna(each_gene) else []

                if loci_type in required_genes:
                    genes_present = all(
                        any(req_gene in gene for gene in each_gene_list)
                        for req_gene in required_genes[loci_type]
                    )
                    if genes_present:
                        filtered_rows.append(row)

            if filtered_rows:
                output_df = pd.DataFrame(filtered_rows)
                output_path = os.path.join(root, "complete_cas_result.csv")
                output_df.to_csv(output_path, index=False)
                print(f"Complete CAS system found: {output_path}")
            else:
                print(f"No complete CAS system in: {cas_info_path}")