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

leader_pattern = r"Leader region\n([A-Za-z]+)"
downstream_pattern = r"Downstream region\n([A-Za-z]+)"

cer_score = r"Certainty Score:\s*([\d\.]+)"

pattern = r'(CRISPR: \d+.*?Certainty Score: [\d.]+)'
pattern_alt = r'(Alternative CRISPR: \d+.*?Certainty Score: [\d.]+)'


# process Bona-Fide_Candidates.txt
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

            #repeat_seq_match = spacers[-1]
            repeat_seq_match =repeats[-1]
            spacer_seq_match = spacers[:-1]

            strain_name = strain_match.group(1)
            start = position_match.groups()[1]
            end = position_match.groups()[2]
            strand = strand_match.group(1)
            sequence_id = sequnence_id_match.group(1)

            leader_region = re.search(leader_pattern, match)
            # leader_sequence = leader_region.group(1)

            if leader_region:
                leader_sequence = leader_region.group(1)

            else:
                leader_sequence = None
                # print("Leader region sequence: None")

            downstream_region = re.search(downstream_pattern, match)
            if downstream_region:
                downstream_sequence = downstream_region.group(1)

            else:
                downstream_sequence = None
                # print("Downstream region sequence: None")

            results.append((strain_name, start, end, strand, repeat_seq_match, spacer_seq_match, leader_sequence, downstream_sequence, sequence_id))
        return results
    else:
        strain_name = strain_match.group(1)
        sequence_id = sequnence_id_match.group(1)
        return strain_name, sequence_id



# process Complete_spacer_dataset.fasta
def parse_spacer_dataset(filepath):
    #pattern = re.compile(r'^>(.*Bona-fide.*)_(-?\d+)_-_(-?\d+)_.*\n([\s\S]+?)(?=>|$)', re.MULTILINE)
    pattern = re.compile(r'^>(.*Bona-fide.*?)(?:_-_)(-?\d+)_-_(-?\d+)_-_(.*)\n([\s\S]+?)(?=>|$)', re.MULTILINE)

    records = []

    with open(filepath, 'r') as file:
        content = file.read()
        matches = pattern.findall(content)
        for match in matches:
            identifier, start, end, repeat, spacers = match

            records.append({
                "Identifier": identifier,
                "Start": start,
                "End": end,
                "Repeat": repeat,
                "Spacer": spacers
            })
    return records

def Collect_results(input_dir, output_csv):

    results = []

    for root, _, files in os.walk(input_dir):

        # if root == input_dir:
        #     continue

        # # folder_name = os.path.basename(root)
        # # file = folder_name
        # if root != input_dir:
        #     folder_name = os.path.basename(root)
        #     file = folder_name

        if root != input_dir:
            parent_folder = os.path.basename(os.path.dirname(root))
            folder_name = os.path.basename(root)

            if parent_folder == os.path.basename(input_dir):
                file = folder_name
            else:
                file = parent_folder

        strain_name = "Unknown"
        strand = "Unknown"
        start = 0
        end = 0
        repeat = "Unknown"
        spacers = "Unknown"
        typing_repeat = "Unknown"
        typing_spacers = "Unknown"
        leader_region = "Unknown"
        downstream_region = "Unknown"
        sequence_id = "Unknown"

        # analysis Complete_spacer_dataset.fasta
        fasta_file = os.path.join(root, "Complete_spacer_dataset.fasta")
        if os.path.exists(fasta_file):
            fasta_records = parse_spacer_dataset(fasta_file)

            for entry in fasta_records:
                start = entry.get('Start')
                end = entry.get('End')
                repeat = entry.get("Repeat")
                spacers = entry.get("Spacer")

                results.append([
                    file,
                    strain_name,
                    start,
                    end,
                    repeat,
                    spacers,
                    strand,
                    typing_repeat,
                    typing_spacers,
                    leader_region,
                    downstream_region,
                    sequence_id
                ])

        candidates_file = os.path.join(root, "Bona-Fide_Candidates.txt")
        if os.path.exists(candidates_file):
            candidate_data = parse_candidates(candidates_file)

            if candidate_data:
                if isinstance(candidate_data, list):
                    #for entry in candidate_data:
                    for (i, row), (entry) in zip(enumerate(results), (candidate_data)):
                        if isinstance(entry, tuple) and len(entry) == 9:
                            strain_name, start_bona, end_bona, strand, repeat_seq, spacers_seq, leader_seq, downstream_seq, sequence_id = entry

                            repeat_length = len(repeat_seq.replace(" ", ""))

                            if strand == "Reversed":
                                corrected_start = int(end_bona) - (repeat_length - 1)
                                corrected_end = int(start_bona) + (repeat_length - 1)
                            else:
                                corrected_start, corrected_end = int(start_bona), int(end_bona)

                            for result in results:

                                if ((int(result[2]) == corrected_start and int(result[3]) == corrected_end)
                                        or (-4 < int(result[2]) - corrected_start < 4 and -4 < int(result[3]) - corrected_end < 4)):
                                # if ((int(result[2]) == corrected_start and int(result[3]) == corrected_end)):

                                    # strain_name
                                    if result[1] == "Unknown":
                                        result[1] = strain_name
                                    # strand
                                    if result[6] == "Unknown":
                                        result[6] = strand
                                    # typing_repeat
                                    if result[7] == "Unknown":
                                        result[7] = repeat_seq
                                    # typing_spacers
                                    if result[8] == "Unknown":
                                        result[8] = spacers_seq
                                    # leader_region
                                    if result[9] == "Unknown":
                                        result[9] = leader_seq
                                    # downstream_region
                                    if result[10] == "Unknown":
                                        result[10] = downstream_seq
                                    if result[11] == "Unknown":
                                        result[11] = sequence_id

                                    # print("result:", result[2], result[3])
                            # print("\n")
                            # bona candidate < 2
                            if len(candidate_data) < 2:
                                alt_file = os.path.join(root, "Alternative_Candidates.txt")
                                alt_data = alternative_candidates(alt_file)
                                if alt_data:
                                    if isinstance(alt_data, list):
                                        for (i, row), (entry_alt) in zip(enumerate(results), (alt_data)):
                                            if isinstance(entry_alt, tuple) and len(
                                                    entry_alt) == 5:  # 确保每个元素是一个包含5个元素的元组
                                                start_alt, end_alt, strand_alt, repeat_alt, spacers_alt = entry_alt  # 解包元组

                                                # 校正位点信息
                                                if strand_alt == "Reversed":
                                                    corrected_start_alt = int(end_alt)
                                                    corrected_end_alt = int(start_alt)
                                                    repat_init = reverse_complement(repeat_alt)
                                                    spacers_init = ' '.join(
                                                        reverse_complement(seq) for seq in spacers_alt)
                                                else:
                                                    corrected_start_alt, corrected_end_alt = int(
                                                        start_alt), int(end_alt)
                                                    repat_init = repeat_alt
                                                    spacers_init = ' '.join(spacers_alt)

                                                if (corrected_end_alt < corrected_start) or (
                                                        corrected_start_alt > corrected_end):
                                                    start = corrected_start_alt
                                                    end = corrected_end_alt
                                                    repeat = repat_init
                                                    spacers = spacers_init
                                                    typing_repeat = repeat_alt
                                                    typing_spacers = spacers_alt
                                                    leader_region = "none"
                                                    downstream_region = "none"

                                                    results.append([
                                                        file,
                                                        strain_name,
                                                        start,
                                                        end,
                                                        repeat,
                                                        spacers,
                                                        strand,
                                                        typing_repeat,
                                                        typing_spacers,
                                                        leader_region,
                                                        downstream_region,
                                                        sequence_id
                                                    ])
                        else:
                            print("error: the tuple format is not as expected", entry)

                    # print("\n")

                else:
                    strain_name, sequence_id = candidate_data

                    strand = "none"
                    start = 0
                    end = 0
                    repeat = "none"
                    spacers = "none"
                    typing_repeat = "none"
                    typing_spacers = "none"
                    leader_region = "none"
                    downstream_region = "none"

                    results.append([
                        file,
                        strain_name,
                        start,
                        end,
                        repeat,
                        spacers,
                        strand,
                        typing_repeat,
                        typing_spacers,
                        leader_region,
                        downstream_region,
                        sequence_id
                    ])


    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File", "Strain Name", "Start", "End", "Repeat Sequence", "Spacer Sequence", "Array Orientation", "Typing Repeat", "Typing Spacers", "Leader Region", "Downstream Region", "Sequence Id"])
        writer.writerows(results)

    print(f"Data has been written to {output_csv}")

#--------------------------------------------------------------------------------
# Alternative_Candidates
def alternative_candidates(file_path):

    with open(file_path, "r") as f:
        content = f.read()

    matches = re.findall(pattern_alt, content, re.DOTALL)

    if matches:
        results = []
        for match in matches:
            position_match = re.search(position_pattern, match)
            strand_match = re.search(strand_pattern, match)

            spacers = re.findall(spacer_repeat_pattern, match)
            repeats = re.findall(repeat_pattern, match)

            repeat_seq_match = repeats[-1]
            spacer_seq_match = spacers[:-1]

            start = position_match.groups()[1]
            end = position_match.groups()[2]
            strand = strand_match.group(1)

            results.append((start, end, strand, repeat_seq_match, spacer_seq_match))
        return results
    else:
        return None

#--------------------------------------------------------------------------------
# initial data filtering
def Data_filtering(input_file, output_file):

    df = pd.read_csv(input_file)

    df = df[df['Repeat Sequence'] != 'none']

    df['Repeat Sequence'] = df['Repeat Sequence'].str.replace(' ', '', regex=False)
    df['Typing Repeat'] = df['Typing Repeat'].str.replace(' ', '', regex=False)

    bad_files = set()

    for index, row in df.iterrows():
        file_id = row['File']

        # process Leader/Downstream Region
        if not ((row['Leader Region'] == 'none') and (row['Downstream Region'] == 'none')):
            if pd.isna(row['Leader Region']) or pd.isna(row['Downstream Region']) or \
                    len(row['Leader Region']) < 60 or len(row['Downstream Region']) < 60:
                bad_files.add(file_id)
                continue

        # process Typing Spacers
        try:
            spacers = ast.literal_eval(row['Typing Spacers'])
            if not isinstance(spacers, list):
                raise ValueError
        except:

            bad_files.add(file_id)
            continue

        for spacer in spacers:
            if 'N' in spacer or len(spacer) < 18 or len(spacer) > 60:
                bad_files.add(file_id)
                break

    df = df[~df['File'].isin(bad_files)]

    for index, row in df.iterrows():
        if row['Array Orientation'] == "Forward":
            df.at[index, 'Leader Region'], df.at[index, 'Downstream Region'] = row['Downstream Region'], row[
                'Leader Region']

    df.to_csv(output_file, index=False)

    print(f"Processed data has been written to {output_file}")
    print('\n')

    file = df['File'].unique()
    print("--------------------------------------------------------------------------------")
    return file

