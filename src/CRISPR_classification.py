import os
import re
import pandas as pd
from utils import create_folder
from utils import needleman_wunsch_similarity
from utils import reverse_complement

repeat_1_reversed = "CGGTTTATCCCCGCTCGCGCGGGGAACAC"
repeat_2_reversed = "CGGTTTATCCCCGCTCGCGCGGGGAACAG"
repeat_36_reversed = "TTTCTAAGCTGCCTGTACGGCAGTGAAC"

repeat_1_forward = reverse_complement(repeat_1_reversed)
repeat_2_forward = reverse_complement(repeat_2_reversed)
repeat_36_forward = reverse_complement(repeat_36_reversed)

# process BLAST results
def get_max_query_positions(df):
    if df.empty:
        return (None, None)
    max_row = df.loc[df['bit_score'].idxmax()]
    return max_row['query_start'], max_row['query_end']

def process_blast_file(blast_file, sequence_id):
    if not os.path.exists(blast_file):
        raise FileNotFoundError("blast file not exist")

    blast_data = pd.read_csv(blast_file, sep='\t', header=None,
                             names=["query_id", "subject_id", "identity", "alignment_length",
                                    "mismatches", "gap_opens", "query_start", "query_end",
                                    "subject_start", "subject_end", "evalue", "bit_score", "source"])

    # cueO、fadM、Y2-aiiA
    cueO_data = blast_data[blast_data['source'].str.contains('cueO', na=False)]
    fadM_data = blast_data[blast_data['source'].str.contains('fadM', na=False)]
    Y2_data = blast_data[blast_data['source'].str.contains('Y2-aiiA', na=False)]
    #print("all data",Y2_data)

    Y2_data.loc[:, 'bit_score'] = pd.to_numeric(Y2_data['bit_score'], errors='coerce')
    cueO_data.loc[:, 'bit_score'] = pd.to_numeric(cueO_data['bit_score'], errors='coerce')
    fadM_data.loc[:, 'bit_score'] = pd.to_numeric(fadM_data['bit_score'], errors='coerce')

    # cueO_max_score = cueO_data.loc[cueO_data['bit_score'].idxmax()]
    # fadM_max_score = fadM_data.loc[fadM_data['bit_score'].idxmax()]
    # Y2_max_score = Y2_data.loc[Y2_data['bit_score'].idxmax()]

    # cueO_query_start, cueO_query_end = cueO_max_score['query_start'], cueO_max_score['query_end']
    # fadM_query_start, fadM_query_end = fadM_max_score['query_start'], fadM_max_score['query_end']
    # Y2_query_start, Y2_query_end = Y2_max_score['query_start'], Y2_max_score['query_end']

    cueO_query_start, cueO_query_end = get_max_query_positions(cueO_data)
    fadM_query_start, fadM_query_end = get_max_query_positions(fadM_data)
    Y2_query_start, Y2_query_end = get_max_query_positions(Y2_data)


    return (cueO_query_start, cueO_query_end), (fadM_query_start, fadM_query_end), (Y2_query_start, Y2_query_end)

#--------------------------------------------------------------------------------
# process cas files
def process_cas_file(file, base_folder):
    if not os.path.exists(base_folder):
        raise FileNotFoundError("cas file not exist")

    results = []
    gcf_folder = None
    #gcf_folder = os.path.join(base_folder, refseq_assembly)
    # for folder in os.listdir(base_folder):
    #     print(folder)
    #     if folder.startswith(refseq_assembly):
    #         gcf_folder = os.path.join(base_folder, folder)
    #         break

    for root, dirs, files in os.walk(base_folder):
        for folder in dirs:
            # print(f"Checking folder: {folder}")
            if folder.startswith(file):
                gcf_folder = os.path.join(root, folder)
                break
        if gcf_folder:
            break


    cas_loci_path = os.path.join(gcf_folder, "cas_info.csv")
    if os.path.exists(cas_loci_path):
        cas_loci_df = pd.read_csv(cas_loci_path)

        for _, row in cas_loci_df.iterrows():

            cas_type = row['loci_type']
            cas_start = row['loci_start']
            cas_end = row['loci_end']
            cas_replicon = row['replicon']

            cas_type = cas_type if pd.notna(cas_type) else 'None'
            cas_start = cas_start if pd.notna(cas_start) else 'None'
            cas_end = cas_end if pd.notna(cas_end) else 'None'
            cas_replicon = cas_replicon if pd.notna(cas_replicon) else 'None'


            results.append((cas_type, cas_start, cas_end, cas_replicon))

        return results

#--------------------------------------------------------------------------------
# CRISPR classification
def parse_range(gene_range):
    return list(map(int, gene_range)) if None not in gene_range else [float('inf'), float('-inf')]

def CRISPR_sort(file, data, blast_base, cas_result, threshold=0.85):
    # cueO_range = list(map(int, cueO))
    # fadM_range = list(map(int, fadM))
    # Y2_range = list(map(int, Y2))

    new_strain_name = None
    results = []

    for _, row in data.iterrows():
        if row['File'] == file:
            CRISPR_strain = row['Strain Name']
            CRISPR_start = int(row['Start'])
            CRISPR_end = int(row['End'])

            CRISPR_repeat = row['Repeat Sequence']
            CRISPR_orientation = row['Array Orientation']
            CRISPR_spacer = row['Typing Spacers']

            CRISPR_sequence_id = row['Sequence Id']

            blast_file = f'{blast_base}/blast_results_{file}.txt'
            cueO, fadM, Y2 = process_blast_file(blast_file, CRISPR_sequence_id)
            cueO_range = parse_range(cueO)
            fadM_range = parse_range(fadM)
            Y2_range = parse_range(Y2)

            # species_name = CRISPR_strain.split()[1]
            match = re.search(r'Cronobacter\s+(\S+)', CRISPR_strain)
            species_name = match.group(1) if match else None

            print(CRISPR_strain, CRISPR_start, CRISPR_end, CRISPR_orientation, CRISPR_sequence_id)
            print(f"Y2-aiiA query_start: {Y2_range[0]}, Y2-aiiA query_end: {Y2_range[1]}")
            print(f"fadM query_start: {fadM_range[0]}, fadM query_end: {fadM_range[1]}")
            print(f"cueO query_start: {cueO_range[0]}, cueO query_end: {cueO_range[1]}")

            found_match = False

            if CRISPR_orientation == "Reversed":

                if not cas_result:
                    print("No cas")

                    if 0 < CRISPR_start - cueO_range[1] < 600:
                        new_strain_name = f"{CRISPR_strain}_crispr3"
                        print("new_strain_name：", new_strain_name)
                    elif ((20000 >= CRISPR_start - Y2_range[1] > -10) and (0 < fadM_range[0] - CRISPR_end < 600)) or (((20000 >= CRISPR_start - Y2_range[1] > -10) or (0 < fadM_range[0] - CRISPR_end < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_reversed) >= threshold):
                        new_strain_name = f"{CRISPR_strain}_crispr2"
                        print("new_strain_name：", new_strain_name)
                    else:
                        new_strain_name = f"{CRISPR_strain}_none"
                        print("new_strain_name：", new_strain_name)

                    results.append({
                        "File": file,
                        "Strain Name": new_strain_name,
                        "Start": CRISPR_start,
                        "End": CRISPR_end,
                        "Array Orientation": CRISPR_orientation,
                        "Typing Spacers": CRISPR_spacer,
                        "Sequence Id": CRISPR_sequence_id,
                        "Species Name": species_name
                    })

                    print('\n')

                elif len(cas_result) > 1:
                    for cas_type, cas_start, cas_end, cas_replicon in cas_result:
                        if CRISPR_sequence_id == cas_replicon:
                            print(f"Cas Type: {cas_type}, Start: {cas_start}, End: {cas_end}, replicon: {cas_replicon}")
                            cas_start = int(cas_start)
                            cas_end = int(cas_end)

                            if (("I-E" in cas_type) and (0 < cas_start - CRISPR_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr1"
                                print("new_strain_name：", new_strain_name)

                            elif (((20000 >= CRISPR_start - Y2_range[1] > -10) and (0 < fadM_range[0] - CRISPR_end < 600)) or (((20000 >= CRISPR_start - Y2_range[1] > -10) or (0 < fadM_range[0] - CRISPR_end < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_reversed) >= threshold) or (("I-E" in cas_type) and (20000 >= CRISPR_start - Y2_range[1] > -10) and (0 < fadM_range[0] - CRISPR_end < 600) and (0 < CRISPR_start - cas_end <= 20000))):
                                new_strain_name = f"{CRISPR_strain}_crispr2"
                                print("new_strain_name：", new_strain_name)

                            elif ((0 < CRISPR_start - cueO_range[1] < 600) or (("I-F" in cas_type) and (0 < cas_start - CRISPR_end < 600) and (0 < CRISPR_start - cueO_range[1] < 600))):
                                new_strain_name = f"{CRISPR_strain}_crispr3"
                                print("new_strain_name：", new_strain_name)

                            elif (("I-F" in cas_type) and ( 0 < CRISPR_start - cas_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr6"
                                print("new_strain_name：", new_strain_name)

                            else:
                                new_strain_name = f"{CRISPR_strain}_none"
                            print("result:", new_strain_name)
                            print('\n')


                            results.append({
                                "File": file,
                                "Strain Name": new_strain_name,
                                "Start": CRISPR_start,
                                "End": CRISPR_end,
                                "Array Orientation": CRISPR_orientation,
                                "Typing Spacers": CRISPR_spacer,
                                "Sequence Id": CRISPR_sequence_id,
                                "Species Name": species_name
                            })
                            found_match = True

                    if not found_match:
                        if 0 < CRISPR_start - cueO_range[1] < 600:
                            new_strain_name = f"{CRISPR_strain}_crispr3"
                            print("new_strain_name：", new_strain_name)
                        elif ((20000 >= CRISPR_start - Y2_range[1] > -10) and (0 < fadM_range[0] - CRISPR_end < 600)) or (((20000 >= CRISPR_start - Y2_range[1] > -10) or (0 < fadM_range[0] - CRISPR_end < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_reversed) >= threshold):
                            new_strain_name = f"{CRISPR_strain}_crispr2"
                            print("new_strain_name：", new_strain_name)
                        else:
                            new_strain_name = f"{CRISPR_strain}_none"
                            print("new_strain_name：", new_strain_name)

                        results.append({
                            "File": file,
                            "Strain Name": new_strain_name,
                            "Start": CRISPR_start,
                            "End": CRISPR_end,
                            "Array Orientation": CRISPR_orientation,
                            "Typing Spacers": CRISPR_spacer,
                            "Sequence Id": CRISPR_sequence_id,
                            "Species Name": species_name
                        })
                        print('\n')


                else:
                    for cas_type, cas_start, cas_end, cas_replicon in cas_result:
                        if CRISPR_sequence_id == cas_replicon:
                            print(f"Cas Type: {cas_type}, Start: {cas_start}, End: {cas_end}, replicon: {cas_replicon}")
                            cas_start = int(cas_start)
                            cas_end = int(cas_end)

                            if (("I-E" in cas_type) and (0 < cas_start - CRISPR_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr1"
                                print("new_strain_name：", new_strain_name)

                            elif (((20000 >= CRISPR_start - Y2_range[1] > -10) and (0 < fadM_range[0] - CRISPR_end < 600)) or (((20000 >= CRISPR_start - Y2_range[1] > -10) or (0 < fadM_range[0] - CRISPR_end < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_reversed) >= threshold) or (("I-E" in cas_type) and (20000 >= CRISPR_start >= Y2_range[1]) and (0 < fadM_range[0] - CRISPR_end < 600) and (0 < CRISPR_start - cas_end <= 20000))):
                                new_strain_name = f"{CRISPR_strain}_crispr2"
                                print("new_strain_name：", new_strain_name)

                            elif ((0 < CRISPR_start - cueO_range[1] < 600) or (("I-F" in cas_type) and (0 < cas_start - CRISPR_end < 600) and (0 < CRISPR_start - cueO_range[1] < 600))):
                                new_strain_name = f"{CRISPR_strain}_crispr3"
                                print("new_strain_name：", new_strain_name)

                            elif (("I-F" in cas_type) and ( 0 < CRISPR_start - cas_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr6"
                                print("new_strain_name：", new_strain_name)

                            else:
                                new_strain_name = f"{CRISPR_strain}_none"
                            print("result:", new_strain_name)
                            print('\n')

                            results.append({
                                "File": file,
                                "Strain Name": new_strain_name,
                                "Start": CRISPR_start,
                                "End": CRISPR_end,
                                "Array Orientation": CRISPR_orientation,
                                "Typing Spacers": CRISPR_spacer,
                                "Sequence Id": CRISPR_sequence_id,
                                "Species Name": species_name
                            })
                            found_match = True

                    if not found_match:
                        if 0 < CRISPR_start - cueO_range[1] < 600:
                            new_strain_name = f"{CRISPR_strain}_crispr3"
                            print("new_strain_name：", new_strain_name)
                        elif ((20000 >= CRISPR_start - Y2_range[1] > -10) and (0 < fadM_range[0] - CRISPR_end < 600)) or (((20000 >= CRISPR_start - Y2_range[1] > -10) or (0 < fadM_range[0] - CRISPR_end < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_reversed) >= threshold):
                            new_strain_name = f"{CRISPR_strain}_crispr2"
                            print("new_strain_name：", new_strain_name)
                        else:
                            new_strain_name = f"{CRISPR_strain}_none"
                            print("new_strain_name：", new_strain_name)

                        results.append({
                            "File": file,
                            "Strain Name": new_strain_name,
                            "Start": CRISPR_start,
                            "End": CRISPR_end,
                            "Array Orientation": CRISPR_orientation,
                            "Typing Spacers": CRISPR_spacer,
                            "Sequence Id": CRISPR_sequence_id,
                            "Species Name": species_name
                        })
                        print('\n')


            else:
                if not cas_result:
                    print("No cas")

                    if 0 < cueO_range[0] - CRISPR_end < 600:
                        new_strain_name = f"{CRISPR_strain}_crispr3"
                        print("new_strain_name：", new_strain_name)
                    elif ((20000 >= Y2_range[0] - CRISPR_end > -10) and (0 < CRISPR_start - fadM_range[1] < 600)) or (((20000 >= Y2_range[0] - CRISPR_end > -10) or (0 < CRISPR_start - fadM_range[1] < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_forward) >= threshold):
                        new_strain_name = f"{CRISPR_strain}_crispr2"
                        print("new_strain_name：", new_strain_name)
                    else:
                        new_strain_name = f"{CRISPR_strain}_none"
                        print("new_strain_name：", new_strain_name)

                    results.append({
                        "File": file,
                        "Strain Name": new_strain_name,
                        "Start": CRISPR_start,
                        "End": CRISPR_end,
                        "Array Orientation": CRISPR_orientation,
                        "Typing Spacers": CRISPR_spacer,
                        "Sequence Id": CRISPR_sequence_id,
                        "Species Name": species_name
                    })
                    print('\n')

                elif len(cas_result) > 1:
                    for cas_type, cas_start, cas_end, cas_replicon in cas_result:
                        if CRISPR_sequence_id == cas_replicon:
                            print(f"Cas Type: {cas_type}, Start: {cas_start}, End: {cas_end}, replicon: {cas_replicon}")
                            cas_start = int(cas_start)
                            cas_end = int(cas_end)

                            #if (("I-E" in cas_type) and ( 0 <= CRISPR_start - cas_end < 500) and (Y2_range[1] < cas_start) and (Y2_range[0] > fadM_range[1])):
                            if (("I-E" in cas_type) and (0 < CRISPR_start - cas_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr1"
                                print("new_strain_name：", new_strain_name)


                            elif (((20000 >= Y2_range[0] - CRISPR_end > -10) and (0 < CRISPR_start - fadM_range[1] < 600)) or (((20000 >= Y2_range[0] - CRISPR_end > -10) or (0 < CRISPR_start - fadM_range[1] < 600)) and needleman_wunsch_similarity(CRISPR_repeat,repeat_2_forward) >= threshold) or (("I-E" in cas_type) and (20000 >= Y2_range[0] - CRISPR_end > -10) and (0 < CRISPR_start - fadM_range[1] < 600) and (0 < cas_start - CRISPR_end <= 20000))):
                                new_strain_name = f"{CRISPR_strain}_crispr2"
                                print("new_strain_name：", new_strain_name)


                            elif ((0 < cueO_range[0] - CRISPR_end < 600) or (("I-F" in cas_type) and (0 < CRISPR_start - cas_end < 600) and (0 < cueO_range[0] - CRISPR_end < 600))):
                                new_strain_name = f"{CRISPR_strain}_crispr3"
                                print("new_strain_name：", new_strain_name)

                            elif (("I-F" in cas_type) and (0 < cas_start - CRISPR_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr6"
                                print("new_strain_name：", new_strain_name)

                            else:
                                new_strain_name = f"{CRISPR_strain}_none"
                            print("result:", new_strain_name)
                            print('\n')

                            results.append({
                                "File": file,
                                "Strain Name": new_strain_name,
                                "Start": CRISPR_start,
                                "End": CRISPR_end,
                                "Array Orientation": CRISPR_orientation,
                                "Typing Spacers": CRISPR_spacer,
                                "Sequence Id": CRISPR_sequence_id,
                                "Species Name": species_name
                            })
                            found_match = True

                    if not found_match:
                        if 0 < cueO_range[0] - CRISPR_end < 600:
                            new_strain_name = f"{CRISPR_strain}_crispr3"
                            print("new_strain_name：", new_strain_name)
                        elif ((20000 >= Y2_range[0] - CRISPR_end > -10) and (0 < CRISPR_start - fadM_range[1] < 600)) or (((20000 >= Y2_range[0] - CRISPR_end > -10) or (0 < CRISPR_start - fadM_range[1] < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_forward) >= threshold):
                            new_strain_name = f"{CRISPR_strain}_crispr2"
                            print("new_strain_name：", new_strain_name)
                        else:
                            new_strain_name = f"{CRISPR_strain}_none"
                            print("new_strain_name：", new_strain_name)

                        results.append({
                            "File": file,
                            "Strain Name": new_strain_name,
                            "Start": CRISPR_start,
                            "End": CRISPR_end,
                            "Array Orientation": CRISPR_orientation,
                            "Typing Spacers": CRISPR_spacer,
                            "Sequence Id": CRISPR_sequence_id,
                            "Species Name": species_name
                        })
                        print('\n')

                else:
                    for cas_type, cas_start, cas_end, cas_replicon in cas_result:
                        if CRISPR_sequence_id == cas_replicon:
                            print(f"Cas Type: {cas_type}, Start: {cas_start}, End: {cas_end}, replicon: {cas_replicon}")
                            cas_start = int(cas_start)
                            cas_end = int(cas_end)

                            #if (("I-E" in cas_type) and ( 0 <= CRISPR_start - cas_end < 500) and (Y2_range[1] < cas_start) and (Y2_range[0] > fadM_range[1])):
                            if (("I-E" in cas_type) and (0 < CRISPR_start - cas_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr1"
                                print("new_strain_name：", new_strain_name)

                            elif (((20000 >= Y2_range[0] - CRISPR_end > -10) and (0 < CRISPR_start - fadM_range[1] < 600)) or (((20000 >= Y2_range[0] - CRISPR_end > -10) or (0 < CRISPR_start - fadM_range[1] < 600)) and needleman_wunsch_similarity(CRISPR_repeat,repeat_2_forward) >= threshold) or (("I-E" in cas_type) and (20000 >= Y2_range[0] >= CRISPR_end) and (0 < CRISPR_start - fadM_range[1] < 600) and (0 < cas_start - CRISPR_end <= 20000))):
                                new_strain_name = f"{CRISPR_strain}_crispr2"
                                print("new_strain_name：", new_strain_name)


                            elif ((0 < cueO_range[0] - CRISPR_end < 600) or (("I-F" in cas_type) and (0 < CRISPR_start - cas_end < 600) and (0 < cueO_range[0] - CRISPR_end < 600))):
                                new_strain_name = f"{CRISPR_strain}_crispr3"
                                print("new_strain_name：", new_strain_name)

                            elif (("I-F" in cas_type) and (0 < cas_start - CRISPR_end < 600)):
                                new_strain_name = f"{CRISPR_strain}_crispr6"
                                print("new_strain_name：", new_strain_name)

                            else:
                                new_strain_name = f"{CRISPR_strain}_none"
                            print("result:", new_strain_name)
                            print('\n')

                            results.append({
                                "File": file,
                                "Strain Name": new_strain_name,
                                "Start": CRISPR_start,
                                "End": CRISPR_end,
                                "Array Orientation": CRISPR_orientation,
                                "Typing Spacers": CRISPR_spacer,
                                "Sequence Id": CRISPR_sequence_id,
                                "Species Name": species_name
                            })
                            found_match = True

                    if not found_match:
                        if 0 < cueO_range[0] - CRISPR_end < 600:
                            new_strain_name = f"{CRISPR_strain}_crispr3"
                            print("new_strain_name：", new_strain_name)
                        elif (20000 >= Y2_range[0] - CRISPR_end > -10) and (0 < CRISPR_start - fadM_range[1] < 600) or (((20000 >= Y2_range[0] - CRISPR_end > -10) or (0 < CRISPR_start - fadM_range[1] < 600)) and needleman_wunsch_similarity(CRISPR_repeat, repeat_2_forward) >= threshold):
                            new_strain_name = f"{CRISPR_strain}_crispr2"
                            print("new_strain_name：", new_strain_name)
                        else:
                            new_strain_name = f"{CRISPR_strain}_none"
                            print("new_strain_name：", new_strain_name)

                        results.append({
                            "File": file,
                            "Strain Name": new_strain_name,
                            "Start": CRISPR_start,
                            "End": CRISPR_end,
                            "Array Orientation": CRISPR_orientation,
                            "Typing Spacers": CRISPR_spacer,
                            "Sequence Id": CRISPR_sequence_id,
                            "Species Name": species_name
                        })
                        print('\n')
    return results


def CRISPR_classification(input_csv, base_folder, output_csv, blast_base):

    data = pd.read_csv(input_csv)

    create_folder(output_csv)

    processed_refseqs = set()

    for idx, row in data.iterrows():
        file = row['File']
        # CRISPR_strain = row['Strain Name']

        if file in processed_refseqs:
            continue

        # blast_file = f'{blast_base}/blast_results_{file}.txt'

        # process blast files
        # cueO_range, fadM_range, Y2_range = process_blast_file(blast_file)

        # process cas files
        cas_result = process_cas_file(file, base_folder)

        # CRISPR classification
        # results = CRISPR_sort(file, data, cas_result, cueO_range, fadM_range, Y2_range)
        results = CRISPR_sort(file, data, blast_base, cas_result)

        if results:
            results_df = pd.DataFrame(results)

            for crispr_type in ['_crispr1', '_crispr2', '_crispr3', '_crispr6']:
                filtered_df = results_df[results_df['Strain Name'].str.contains(crispr_type)]
                file_name = os.path.join(output_csv, f'CRISPR_results_{crispr_type[1:]}.csv')
                filtered_df.to_csv(file_name, mode='a', header=not os.path.exists(file_name), index=False)

        else:
            print("No results found.")

        processed_refseqs.add(file)
        print("--------------------------------------------------------------------------------")