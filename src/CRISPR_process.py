import os
import re
import csv
import pandas as pd


# regular expression
strain_name_pattern = r">[^ ]+ ([^\n,]+)"
position_pattern = r"CRISPR: (\d+), (\d+)-(\d+),"
strand_pattern = r"Strand: (\w+)"

spacer_repeat_pattern = r'\s+([A-Za-z]+(?:\s[A-Za-z]+)*)\s+s:'
repeat_pattern = r'\s*([ATGCN\s]+)\s+s:'

leader_pattern = r"Leader region\n([A-Za-z]+)"
downstream_pattern = r"Downstream region\n([A-Za-z]+)"

cer_score = r"Certainty Score:\s*([\d\.]+)"

pattern = r'(CRISPR: \d+.*?Certainty Score: [\d.]+)'

# process Bona-Fide_Candidates.txt
def parse_candidates(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    strain_match = re.search(strain_name_pattern, content)
    matches = re.findall(pattern, content, re.DOTALL)

    if strain_match and matches:
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

            results.append((strain_name, start, end, strand, repeat_seq_match, spacer_seq_match, leader_sequence, downstream_sequence))
        return results
    else:
        strain_name = strain_match.group(1)
        return strain_name



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
                    downstream_region
                ])

        candidates_file = os.path.join(root, "Bona-Fide_Candidates.txt")
        if os.path.exists(candidates_file):
            candidate_data = parse_candidates(candidates_file)

            if candidate_data:
                if isinstance(candidate_data, list):
                    #for entry in candidate_data:
                    for (i, row), (entry) in zip(enumerate(results), (candidate_data)):
                        if isinstance(entry, tuple) and len(entry) == 8:
                            strain_name, start_bona, end_bona, strand, repeat_seq, spacers_seq, leader_seq, downstream_seq = entry

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

                                    # print("result:", result[2], result[3])
                            # print("\n")
                        else:
                            print("error: the tuple format is not as expected", entry)

                    # print("\n")

                else:
                    strain_name = candidate_data

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
                        downstream_region
                    ])


    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File", "Strain Name", "Start", "End", "Repeat Sequence", "Spacer Sequence", "Array Orientation", "Typing Repeat", "Typing Spacers", "Leader Region", "Downstream Region"])
        writer.writerows(results)

    print(f"Data has been written to {output_csv}")

#--------------------------------------------------------------------------------
# initial data filtering
def Data_filtering(input_file, output_file):

    df = pd.read_csv(input_file)

    df = df[df['Repeat Sequence'] != 'none']

    df = df[(df['Leader Region'].notna()) & (df['Downstream Region'].notna())]
    df = df[(df['Leader Region'].str.len() >= 60) & (df['Downstream Region'].str.len() >= 60)]

    df = df[~df['Typing Spacers'].str.contains('N', na=False)]

    # df['Typing Spacers'] = df['Typing Spacers'].apply(
    #     lambda spacers: [s for s in eval(spacers) if 'N' not in s] if pd.notna(spacers) else spacers
    # )

    df['Typing Repeat'] = df['Typing Repeat'].str.replace(' ', '', regex=False)

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

