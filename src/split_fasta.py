import os

def split_fna(input_fna, output_folder="cas"):

    input_basename = os.path.splitext(os.path.basename(input_fna))[0]

    base_folder = os.path.join(output_folder, input_basename)

    if os.path.exists(base_folder) and any(fname.endswith('.fna') for fname in os.listdir(base_folder)):
        print(f"Skip {input_fna} because {base_folder} already exists and contains a '.fna' file")
        return

    os.makedirs(base_folder, exist_ok=True)

    with open(input_fna, 'r') as infile:
        sequence = ''
        header = None
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:

                    output_path = os.path.join(base_folder, f"{header}.fna")
                    with open(output_path, 'w') as out_file:
                        out_file.write(f">{header}\n{sequence}\n")
                header = line[1:]
                sequence = ''
            else:
                sequence += line

        if header:
            output_path = os.path.join(base_folder, f"{header}.fna")
            with open(output_path, 'w') as out_file:
                out_file.write(f">{header}\n{sequence}\n")

    print(f"File split to directory: {base_folder}")

def find_and_split_fna(root_folder):
    # script_dir = os.path.dirname(os.path.abspath(__file__))
    output_folder = os.path.join(root_folder, "ct_output/cas")
    os.makedirs(output_folder, exist_ok=True)

    # for dirpath, _, filenames in os.walk(root_folder):
    #     for filename in filenames:
    #         if filename.endswith('.fna'):
    #             input_fna = os.path.join(dirpath, filename)
    #             print(f"Processing file: {input_fna}")
    #             split_fna(input_fna, output_folder)

    for file_name in os.listdir(root_folder):
        file_path = os.path.join(root_folder, file_name)

        if os.path.isfile(file_path) and file_name.endswith('.fna'):
            split_fna(file_path, output_folder)

