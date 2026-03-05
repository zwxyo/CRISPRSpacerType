import os
import json
import csv
import re
from utils import reverse_complement

def extract_spacers(json_file, folder_name):
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    results = []

    for seq in data.get("Sequences", []):
        desc = seq.get("Description", "")
        match = re.match(r"(.+?)(?:,|$)", desc)
        organism_name = match.group(1) if match else desc
        version = seq.get("Version")

        for crispr in seq.get("Crisprs", []):
            if crispr.get("Evidence_Level") in [3, 4]:
            # if crispr.get("Evidence_Level") == 4:
                orientation = crispr.get("Potential_Orientation", "ND")
                crispr_dir = crispr.get("CRISPRDirection", "ND")

                start = crispr.get("Start")
                end = crispr.get("End")
                # sequence_id = crispr.get("Name")
                sequence_id = version
                dr_consensus = crispr.get("DR_Consensus", "")

                if orientation == "ND" and crispr_dir != "ND":
                    orientation = crispr_dir

                if orientation == "-":
                    orientation = "Reversed"
                elif orientation == "+":
                    orientation = "Forward"

                left_flank = ""
                right_flank = ""
                spacers_list = []

            for region in crispr.get("Regions", []):
                region_type = region.get("Type")
                region_seq = region.get("Sequence", "")

                if region_type == "LeftFLANK":
                    left_flank = region_seq
                elif region_type == "RightFLANK":
                    right_flank = region_seq
                elif region_type == "Spacer":
                    spacers_list.append(region_seq)

            joined_spacers = " ".join(spacers_list)

            if orientation == "Reversed":
                typing_repeat = reverse_complement(dr_consensus)
                typing_spacers = [reverse_complement(s) for s in reversed(spacers_list)]
            elif orientation == "Forward":
                typing_repeat = dr_consensus
                typing_spacers = spacers_list
            else:
                typing_repeat = "none"
                typing_spacers = "none"

            if joined_spacers:
                # for region in crispr.get("Regions", []):
                #     if region.get("Type") == "Spacer":
                results.append({
                    "File": folder_name,
                    "Strain Name": organism_name,
                    "Start": start,
                    "End": end,
                    "Sequence Id": sequence_id,
                    "Array Orientation": orientation,
                    "Repeat Sequence": dr_consensus,
                    "Spacer Sequence": joined_spacers,
                    "Typing Repeat": typing_repeat,
                    "Typing Spacers": typing_spacers,
                    "Leader Region": left_flank,
                    "Downstream Region": right_flank
                })
    return results


def batch_extract(root_dir, output_csv):
    all_results = []

    for root, dirs, files in os.walk(root_dir):
        if "result.json" in files:
            json_path = os.path.join(root, "result.json")
            folder_name = os.path.basename(root)
            print(f"Processing {folder_name}...\n")
            spacers = extract_spacers(json_path, folder_name)
            all_results.extend(spacers)

    output_dir = os.path.dirname(output_csv)

    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ["File", "Strain Name", "Start", "End", "Repeat Sequence", "Spacer Sequence", "Typing Repeat", "Typing Spacers", "Array Orientation", "Leader Region", "Downstream Region", "Sequence Id"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print(f"Data has been written to {output_csv}")


def extract_cas_info(json_file):
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    results = []
    loci_counter = 1

    for seq in data.get("Sequences", []):
        replicon = seq.get("Version", seq.get("Id", "Unknown"))

        for cas_locus in seq.get("Cas", []):
            loci_start = cas_locus.get("Start")
            loci_end = cas_locus.get("End")

            raw_type = cas_locus.get("Type", "")
            match = re.search(r"Subtype-([A-Za-z0-9\-]+)", raw_type)
            if match:
                loci_type = match.group(1)
            else:
                loci_type = raw_type

            genes_list = []
            for gene in cas_locus.get("Genes", []):
                sub_type = gene.get("Sub_type")
                if sub_type:
                    genes_list.append(sub_type)

            each_gene = ",".join(genes_list)

            results.append({
                "loci_ID": loci_counter,
                "loci_type": loci_type,
                "loci_start": loci_start,
                "loci_end": loci_end,
                "each_gene": each_gene,
                "replicon": replicon
            })

            loci_counter += 1

    return results


def batch_extract_cas(root_dir, output_base_dir="cas"):
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    processed_count = 0

    for root, dirs, files in os.walk(root_dir):
        if "result.json" in files:
            json_path = os.path.join(root, "result.json")
            folder_name = os.path.basename(root)
            print(f"Processing {folder_name}...\n")

            cas_data = extract_cas_info(json_path)

            target_folder = os.path.join(output_base_dir, folder_name)
            if not os.path.exists(target_folder):
                os.makedirs(target_folder)

            output_csv = os.path.join(target_folder, "cas_info.csv")

            with open(output_csv, 'w', newline='') as csvfile:
                fieldnames = ["loci_ID", "loci_type", "loci_start", "loci_end", "each_gene", "replicon"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

                writer.writeheader()
                writer.writerows(cas_data)

            processed_count += 1

