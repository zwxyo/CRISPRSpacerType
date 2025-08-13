import os
import sys
import pandas as pd
from utils import create_folder
import csv
import ast
import glob
from utils import needleman_wunsch_similarity
from collections import defaultdict

def merge_rows(row1, row2):

    start = min(row1['Start'], row2['Start'])
    end = max(row1['End'], row2['End'])

    spacers1 = ast.literal_eval(row1['Typing Spacers'])
    spacers2 = ast.literal_eval(row2['Typing Spacers'])

    max_k = 0

    max_len = min(len(spacers1), len(spacers2))
    for k in range(1, max_len+1):
        if spacers1[-k:] == spacers2[:k]:
            max_k = k

    merged_spacers = spacers1 + spacers2[max_k:]

    row1['Start'] = start
    row1['End'] = end
    #row1['Typing Spacers'] = str(merged_spacers)
    row1['Typing Spacers'] = merged_spacers
    return row1

def classification_file_process(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith('.csv'):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path)
            # print(f"Processed: {filename}, original rows: {len(df)}")
            df = df.drop_duplicates()
            # print(f"final rows: {len(df)}")

            df['Start'] = df['Start'].astype(int)
            df['End'] = df['End'].astype(int)

            result_rows = []

            for _, group in df.groupby('File'):
                group = group.reset_index(drop=True)
                used = [False] * len(group)

                for i in range(len(group)):
                    if used[i]:
                        continue
                    row_i = group.loc[i]
                    merged_row = row_i.copy()

                    for j in range(i + 1, len(group)):
                        if used[j]:
                            continue
                        row_j = group.loc[j]

                        same_seqid = row_i['Sequence Id'] == row_j['Sequence Id']
                        start_i, end_i = row_i['Start'], row_i['End']
                        start_j, end_j = row_j['Start'], row_j['End']
                        overlap = (
                            same_seqid and (
                                abs(start_j - end_i) < 150
                                or
                                abs(start_i - end_j) < 150
                            )
                        )
                        if overlap:
                            merged_row = merge_rows(merged_row, row_j)
                            used[j] = True

                    used[i] = True
                    result_rows.append(merged_row)

            result_rows = pd.DataFrame(result_rows) if result_rows else pd.DataFrame()

            final_rows = []
            if not result_rows.empty:

                result_rows['Spacer Count'] = result_rows['Typing Spacers'].apply(
                    lambda x: len(ast.literal_eval(x)) if isinstance(x, str) else len(x)
                )

                for file_name, group in result_rows.groupby('File'):
                    max_row = group.loc[group['Spacer Count'].idxmax()]
                    final_rows.append(max_row.to_dict())

            df_final = pd.DataFrame(final_rows).drop(columns=['Spacer Count'], errors='ignore')
            if not df_final.empty:
                df_final.to_csv(file_path, index=False)
                #print(f"Saved: {filename}, final rows: {len(df_final)}")
            else:
                print(f"No valid rows for {filename}")

# Classify CRISPR classification files by species
def split_csv_by_species(folder_path, new_folder_name):

    new_folder_path = os.path.join(folder_path, new_folder_name)

    os.makedirs(new_folder_path, exist_ok=True)

    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    for csv_file in csv_files:

        file_path = os.path.join(folder_path, csv_file)
        df = pd.read_csv(file_path)

        for species_name, group in df.groupby('Species Name'):

            file_number = csv_file.split('_')[-1].split('.')[0]
            new_file_name = f"{species_name}_{file_number}.csv"
            new_file_path = os.path.join(new_folder_path, new_file_name)

            group.to_csv(new_file_path, index=False)

            print(f"Saved file: {new_file_path}")

#--------------------------------------------------------------------------------
# spacer number
def load_spacers(spacer_file):

    spacer_dict = {}
    with open(spacer_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            spacer_dict[row['Spacer']] = row['Key']
    return spacer_dict


def process_spacer_csv(input_folder, spacer_folder, output_folder):

    for file_name in os.listdir(input_folder):
        if file_name.endswith('.csv'):
            file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name)

            spacer_file_path = os.path.join(spacer_folder, file_name)

            if os.path.exists(spacer_file_path):

                spacer_dict = load_spacers(spacer_file_path)

                with open(file_path, 'r') as f_in, open(output_file_path, 'w', newline='') as f_out:
                    reader = csv.DictReader(f_in)
                    fieldnames = reader.fieldnames + ['Spacers Order']
                    writer = csv.DictWriter(f_out, fieldnames=fieldnames)
                    writer.writeheader()

                    for row in reader:
                        spacers = ast.literal_eval(row['Typing Spacers'])
                        numbered_spacers = []

                        for spacer in spacers:

                            spacer_clean = spacer.replace(" ", "")
                            best_match = 'Not Found'
                            best_similarity = 0

                            for known_spacer, key in spacer_dict.items():
                                similarity = needleman_wunsch_similarity(spacer_clean, known_spacer)
                                if similarity > best_similarity:
                                    best_similarity = similarity
                                    best_match = key if similarity >= 0.84 else 'Not Found'

                            numbered_spacers.append(best_match)

                        row['Spacers Order'] = ', '.join(numbered_spacers)
                        writer.writerow(row)
            else:
                print(f"Warning: No spacer mapping file found for {file_name}")

def spacer_number(input_folder, spacer_folder, output_folder):

    create_folder(output_folder)
    process_spacer_csv(input_folder, spacer_folder, output_folder)

#--------------------------------------------------------------------------------
# serial number
def load_serials(serial_file):
    serial_dict = {}
    for file in glob.glob(os.path.join(serial_file, "*.csv")):
        with open(file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                serial_dict[row['Spacers Order']] = row['serial number']
    return serial_dict


def serial_number(input_folder, serial_folder, output_folder):
    create_folder(output_folder)
    serial_dict = load_serials(serial_folder)

    for file in glob.glob(os.path.join(input_folder, "*.csv")):
        output_file = os.path.join(output_folder, os.path.basename(file))
        with open(file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
            reader = csv.DictReader(f_in)
            fieldnames = reader.fieldnames + ['Serial Number']
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()

            for row in reader:
                row['Serial Number'] = serial_dict.get(row['Spacers Order'], 'Not Found')
                writer.writerow(row)
    print("Serial number processing completed!")

#--------------------------------------------------------------------------------
# ct number
def parse_strain_name(strain_name):
    parts = strain_name.rsplit('_', 1)
    if len(parts) == 2 and parts[1][6:].isdigit():
        return parts[0], int(parts[1][6:])
    return strain_name, None

def load_spacer_mapping(database_dir, species, crispr_type):
    spacer_file = os.path.join(f'{database_dir}/spacer', f'{species}_{crispr_type}.csv')
    if not os.path.exists(spacer_file):
        return None
    spacer_dict = {}
    with open(spacer_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            spacer_dict[row['Spacer']] = row['Key']
    return spacer_dict

def find_serial_number(database_dir, species, crispr_type, numbers_str):
    serial_file = os.path.join(f'{database_dir}', 'serial', f'{species}_crispr{crispr_type}.csv')
    if not os.path.exists(serial_file):
        return None
    with open(serial_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['Spacers Order'] == numbers_str:
                return row[f'serial number']
    return None

def load_existing_ct_numbers(ct_folder):
    print("ct_folder", ct_folder)
    ct_data = defaultdict(dict)
    ct_dir = ct_folder
    if not os.path.exists(ct_dir):
        return ct_data
    for file_path in glob.glob(os.path.join(ct_dir, '*.csv')):
        species = os.path.basename(file_path).replace('.csv', '')
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = (row['CRISPR1-Serial number'], row['CRISPR2-Serial number'],
                       row['CRISPR3-Serial number'], row['CRISPR6-Serial number'])
                ct_data[species][key] = row['CT_number']
    return ct_data


def process_crispr_file(database_dir, file_path, strain_data):
    base_name = os.path.basename(file_path)
    crispr_type = int(base_name.split('_')[-1].split('.')[0][6:])

    with open(file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row['File']:
                continue
            file_id = row['File']

            species = row['Species Name']
            crispr_type_all = "crispr" + str(crispr_type)

            spacer_map = load_spacer_mapping(database_dir, species, crispr_type_all)
            if not spacer_map:
                continue

            try:
                spacers = ast.literal_eval(row['Typing Spacers'].replace(" ", ""))
                numbers = [spacer_map[s] for s in spacers if s in spacer_map]
            except:
                continue
            if len(numbers) != len(spacers):
                continue

            numbers_str = ', '.join(numbers)
            serial = find_serial_number(database_dir, species, crispr_type, numbers_str)

            strain_data[file_id]['Species'] = species
            if crispr_type in [1, 2, 3, 6]:
                strain_data[file_id][f'CRISPR{crispr_type}'] = serial


def process_all_crispr_results(database_dir, folder_path):

    strain_data = defaultdict(lambda: {
        'Species': None,
        'CRISPR1': None, 'CRISPR2': None,
        'CRISPR3': None, 'CRISPR6': None
    })
    for file_path in glob.glob(f'{folder_path}/CRISPR_results_*.csv'):
        process_crispr_file(database_dir, file_path, strain_data)
    return strain_data

#---------------------------------------------------------------
def load_serial_data(serial_dir):
    serial_data = defaultdict(lambda: {
        'Species': None,
        'CRISPR1': None,
        'CRISPR2': None,
        'CRISPR3': None,
        'CRISPR6': None
    })
    for serial_file in glob.glob(os.path.join(serial_dir, '*.csv')):
        filename = os.path.basename(serial_file)
        base_name = os.path.splitext(filename)[0]
        if '_crispr' not in base_name:
            continue
        try:
            species_part, crispr_part = base_name.rsplit('_crispr', 1)
            crispr_type = int(crispr_part)
            if crispr_type not in [1, 2, 3, 6]:
                continue
        except ValueError:
            continue

        with open(serial_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                file_id = row.get('File')
                if not file_id:
                    continue
                serial_number = row.get('Serial Number', '')
                serial_data[file_id]['Species'] = species_part
                key = f'CRISPR{crispr_type}'
                serial_data[file_id][key] = serial_number
    return serial_data

#---------------------------------------------------------------
def assign_ct_numbers(strain_data, existing_ct_data):
    for strain, data in strain_data.items():
        if not data['Species']:
            continue
        species = data['Species']
        crispr_combination = (
            data['CRISPR1'] or '',
            data['CRISPR2'] or '',
            data['CRISPR3'] or '',
            data['CRISPR6'] or ''
        )
        if species in existing_ct_data and crispr_combination in existing_ct_data[species]:
            strain_data[strain]['CT_number'] = existing_ct_data[species][crispr_combination]
        else:
            strain_data[strain]['CT_number'] = 'Unknown'

# def output_ct_results(strain_data, output_folder):
#     os.makedirs(output_folder, exist_ok=True)
#     species_groups = defaultdict(list)
#     for strain, data in strain_data.items():
#         if not data['Species']:
#             continue
#         species_groups[data['Species']].append({
#             'Strain': strain,
#             'CT_number': data['CT_number'],
#             'CRISPR1': data['CRISPR1'] or '',
#             'CRISPR2': data['CRISPR2'] or '',
#             'CRISPR3': data['CRISPR3'] or '',
#             'CRISPR6': data['CRISPR6'] or ''
#         })
#     for species, records in species_groups.items():
#         output_file = os.path.join(output_folder, f'{species}.csv')
#         fields = ['Strain', 'CT_number', 'CRISPR1-Serial number', 'CRISPR2-Serial number',
#                   'CRISPR3-Serial number', 'CRISPR6-Serial number']
#         with open(output_file, 'w', newline='') as f:
#             writer = csv.DictWriter(f, fieldnames=fields)
#             writer.writeheader()
#             for record in records:
#                 writer.writerow({
#                     'Strain': record['Strain'],
#                     'CT_number': record['CT_number'],
#                     'CRISPR1-Serial number': record['CRISPR1'],
#                     'CRISPR2-Serial number': record['CRISPR2'],
#                     'CRISPR3-Serial number': record['CRISPR3'],
#                     'CRISPR6-Serial number': record['CRISPR6']
#                 })
#     print("Processing completed! Results saved in 'user_output/spacer/ct' directory.")

def output_ct_results(strain_data, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    species_groups = defaultdict(list)

    for file_id, data in strain_data.items():
        if not data['Species']:
            continue
        species_groups[data['Species']].append({
            'File': file_id,
            'CT_number': data['CT_number'],
            'CRISPR1': data['CRISPR1'] or '',
            'CRISPR2': data['CRISPR2'] or '',
            'CRISPR3': data['CRISPR3'] or '',
            'CRISPR6': data['CRISPR6'] or ''
        })

    for species, records in species_groups.items():
        output_file = os.path.join(output_folder, f'{species}.csv')
        fields = ['File', 'CT_number', 'CRISPR1-Serial number', 'CRISPR2-Serial number',
                  'CRISPR3-Serial number', 'CRISPR6-Serial number']
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for record in records:
                writer.writerow({
                    'File': record['File'],
                    'CT_number': record['CT_number'],
                    'CRISPR1-Serial number': record['CRISPR1'],
                    'CRISPR2-Serial number': record['CRISPR2'],
                    'CRISPR3-Serial number': record['CRISPR3'],
                    'CRISPR6-Serial number': record['CRISPR6']
                })
    print("Processing completed! Results saved in 'user_output/spacer/ct' directory.")


def ct_number(database_dir, serial_folder, folder_path, ct_folder, output_folder):
    existing_ct_data = load_existing_ct_numbers(ct_folder)
    # strain_data = process_all_crispr_results(database_dir, folder_path)
    strain_data = load_serial_data(serial_folder)
    assign_ct_numbers(strain_data, existing_ct_data)
    output_ct_results(strain_data, output_folder)

#--------------------------------------------------------------------------------
def summary(folder_path, output_spacer, output_serial, output_ct):

    script_dir = os.path.dirname(os.path.abspath(__file__))

    parent_dir = os.path.dirname(script_dir)

    database_dir = os.path.join(parent_dir, "database/Ct_db_Cronobacter")

    classification_file_process(folder_path)

    if not os.path.exists(folder_path):
        print(f"Error: The folder '{folder_path}' does not exist.")
        sys.exit(1)

    elif not os.listdir(folder_path):
        print(f"Error: The folder '{folder_path}' is empty.")
        sys.exit(1)

    # Classification by species
    # folder_path = 'user_output/CRISPR_sort'
    new_folder_name = 'Cronobacter_Genus'
    split_csv_by_species(folder_path, new_folder_name)
    print('\n')

    # spacer
    print("Numbering spacers is underway")
    input_spacer = f'{folder_path}/Cronobacter_Genus'
    # output_spacer = 'user_output/spacer/spacer_order'
    spacer_folder = f'{database_dir}/spacer'
    spacer_number(input_spacer, spacer_folder, output_spacer)
    print('\n')

    # serial
    print("Numbering serials is underway")
    input_serial = output_spacer
    # output_serial = 'user_output/spacer/serial'
    serial_folder = f'{database_dir}/serial'
    serial_number(input_serial, serial_folder, output_serial)
    print('\n')

    # ct
    print("Numbering ct is underway")
    # output_ct = 'user_output/spacer/ct'
    ct_folder = f'{database_dir}/ct'
    ct_number(database_dir, output_serial, folder_path, ct_folder, output_ct)
    print('\n')
    print("--------------------------------------------------------------------------------")

def summary_pcr(folder_path, output_spacer, output_serial, output_ct):

    script_dir = os.path.dirname(os.path.abspath(__file__))

    parent_dir = os.path.dirname(script_dir)

    database_dir = os.path.join(parent_dir, "database/Ct_db_Cronobacter")

    # classification_file_process(folder_path)

    if not os.path.exists(folder_path):
        print(f"Error: The folder '{folder_path}' does not exist.")
        sys.exit(1)

    elif not os.listdir(folder_path):
        print(f"Error: The folder '{folder_path}' is empty.")
        sys.exit(1)

    # Classification by species
    # new_folder_name = 'Cronobacter_Genus'
    # split_csv_by_species(folder_path, new_folder_name)
    # print('\n')

    # spacer
    print("Numbering spacers is underway")
    input_spacer = f'{folder_path}/Cronobacter_Genus'
    # output_spacer = 'user_output/spacer/spacer_order'
    spacer_folder = f'{database_dir}/spacer'
    spacer_number(input_spacer, spacer_folder, output_spacer)
    print('\n')

    # serial
    print("Numbering serials is underway")
    input_serial = output_spacer
    # output_serial = 'user_output/spacer/serial'
    serial_folder = f'{database_dir}/serial'
    serial_number(input_serial, serial_folder, output_serial)
    print('\n')

    # ct
    print("Numbering ct is underway")
    # output_ct = 'user_output/spacer/ct'
    ct_folder = f'{database_dir}/ct'
    ct_number(database_dir, output_serial, folder_path, ct_folder, output_ct)
    print('\n')
    print("--------------------------------------------------------------------------------")
