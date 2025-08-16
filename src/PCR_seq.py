import os
import re
import csv
import ast
import pandas as pd
from utils import *

# regular expression
sequence_id_pattern = r"^>(\S+)"
strain_name_pattern = r">[^ ]+ ([^\n,]+)"
position_pattern = r"CRISPR: (\d+), (\d+)-(\d+),"
strand_pattern = r"Strand: (\w+)"

spacer_repeat_pattern = r'\s+([A-Za-z]+(?:\s[A-Za-z]+)*)\s+s:'
repeat_pattern = r'\s*([ATGCN\s]+)\s+s:'

pattern = r'(CRISPR: \d+.*?Certainty Score: [\d.]+)'
pattern_alt = r'(Alternative CRISPR: \d+.*?Certainty Score: [\d.]+)'

def parse_candidates(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    strain_match = re.search(strain_name_pattern, content)
    sequnence_id_match = re.search(sequence_id_pattern, content)
    matches = re.findall(pattern, content, re.DOTALL)

    if strain_match and sequnence_id_match and matches:
        results = []
        for match in matches:
            position_match = re.search(position_pattern, match)
            strand_match = re.search(strand_pattern, match)

            spacers = re.findall(spacer_repeat_pattern, match)
            repeats = re.findall(repeat_pattern, match)

            repeat_seq_match =repeats[-1]
            spacer_seq_match = spacers[:-1]

            strain_name = strain_match.group(1)
            start = position_match.groups()[1]
            end = position_match.groups()[2]
            strand = strand_match.group(1)
            sequence_id = sequnence_id_match.group(1)


            results.append((strain_name, strand, repeat_seq_match, spacer_seq_match, sequence_id))
        return results
    else:
        strain_name = strain_match.group(1)
        sequence_id = sequnence_id_match.group(1)
        return strain_name, sequence_id

def parse_spacer_dataset(filepath):
    pattern = re.compile(r'^>(.*Bona-fide.*?)(?:_-_)(-?\d+)_-_(-?\d+)_-_(.*)\n([\s\S]+?)(?=>|$)', re.MULTILINE)

    records = []

    with open(filepath, 'r') as file:
        content = file.read()
        matches = pattern.findall(content)

        for match in matches:
            identifier, start, end, repeat, spacers = match
            sequence_id_match = re.match(r'^(\S+)', identifier)
            sequence_id = sequence_id_match.group(1)

            records.append({
                "Sequence_id": sequence_id,
                "Repeat": repeat,
                "Spacer": spacers
            })

    return records


def PCR_process(input_dir, output_csv):

    results = []

    for root, _, files in os.walk(input_dir):

        if root != input_dir:
            print(root)

            parent_folder = os.path.basename(os.path.dirname(root))
            folder_name = os.path.basename(root)

            if parent_folder == os.path.basename(input_dir):
                file = folder_name
            else:
                file = parent_folder

        strain_name = "Unknown"
        strand = "Unknown"
        repeat = "Unknown"
        spacers = "Unknown"
        typing_repeat = "Unknown"
        typing_spacers = "Unknown"
        sequence_id = "Unknown"


        fasta_file = os.path.join(root, "Complete_spacer_dataset.fasta")
        if os.path.exists(fasta_file):
            fasta_records = parse_spacer_dataset(fasta_file)

            for entry in fasta_records:
                sequence_id = entry["Sequence_id"]
                repeat = entry.get("Repeat")
                spacers = entry.get("Spacer")

                results.append([
                    file,
                    strain_name,
                    repeat,
                    spacers,
                    strand,
                    typing_repeat,
                    typing_spacers,
                    sequence_id
                ])

        candidates_file = os.path.join(root, "Bona-Fide_Candidates.txt")

        if os.path.exists(candidates_file):
            candidate_data = parse_candidates(candidates_file)

            if candidate_data:
                if isinstance(candidate_data, list):
                    for (i, row), (entry) in zip(enumerate(results), (candidate_data)):
                        if isinstance(entry, tuple) and len(entry) == 5:
                            strain_name, strand, repeat_seq, spacers_seq, sequence_id = entry

                            for result in results:
                                if result[7] == sequence_id:
                                    # strain_name
                                    if result[1] == "Unknown":
                                        result[1] = strain_name
                                    # strand
                                    if result[4] == "Unknown":
                                        result[4] = strand
                                    # typing_repeat
                                    if result[5] == "Unknown":
                                        result[5] = repeat_seq
                                    # typing_spacers
                                    if result[6] == "Unknown":
                                        result[6] = spacers_seq

                        else:
                            print("error: the tuple format is not as expected", entry)

                else:
                    strain_name, sequence_id = candidate_data

                    strand = "none"
                    repeat = "none"
                    spacers = "none"
                    typing_repeat = "none"
                    typing_spacers = "none"

                    results.append([
                        file,
                        strain_name,
                        repeat,
                        spacers,
                        strand,
                        typing_repeat,
                        typing_spacers,
                        sequence_id
                    ])


    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File", "Strain Name", "Repeat Sequence", "Spacer Sequence", "Array Orientation", "Typing Repeat", "Typing Spacers", "Sequence Id"])
        writer.writerows(results)


    print(f"Data has been written to {output_csv}")

#--------------------------------------------------------------------------------
def extract_species(strain_name):
    match = re.search(r'Cronobacter\s+(\w+)', strain_name)
    if match:
        return match.group(1).lower()
    else:
        return 'unknown'

def extract_crispr_type(seq_id):
    match = re.match(r'(CRISPR\d+)_', seq_id)
    if match:
        return match.group(1).lower()
    else:
        return 'unknown'

def pre_typing_process(input_csv, output_folder):
    data = pd.read_csv(input_csv)
    create_folder(output_folder)

    data['Species Name'] = data['Strain Name'].apply(extract_species)
    data['CRISPR Type'] = data['Sequence Id'].apply(extract_crispr_type)


    for (species, crispr_type), group in data.groupby(['Species Name', 'CRISPR Type']):
        output_file = os.path.join(output_folder, f'{species}_{crispr_type}.csv')
        group.to_csv(output_file, index=False)

