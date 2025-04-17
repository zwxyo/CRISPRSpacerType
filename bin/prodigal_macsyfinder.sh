#!/bin/bash


#input_folder="external/user_input"
input_folder="$1"


script_dir=$(cd "$(dirname "$0")" && pwd)
# python module path
python_module_dir="$script_dir/../src"
#split_fasta_script="$python_module_dir/split_fasta.py"
#cas_loci_script="$python_module_dir/Cas_loci.py"

# cas path
root_folder="$input_folder/ct_output/cas"

for dir in "$root_folder"/*/; do
    # check if there are cas_info.csv files
    if [ ! -f "$dir/cas_info.csv" ]; then
        echo "Delete folder: $dir"
        rm -rf "$dir"
    fi
done


#python3 split_fasta.py "$input_folder"
# python3 "import sys; sys.path.append('$python_module_dir'); sys.argv = ['split_fasta.py', '$input_folder']; exec(open('$split_fasta_script').read())"
python3 -c "import sys; sys.path.append('$python_module_dir'); from split_fasta import find_and_split_fna; find_and_split_fna('$input_folder')"


find "$root_folder" -type f -name "*.fna" | while read -r file; do

    base_name=$(basename "$file" .fna)

    folder_path=$(dirname "$file")

    output_folder="$folder_path/separated_proteins"
    mkdir -p "$output_folder"

    output_file="$output_folder/$base_name.faa"
    output_file1="$output_folder/$base_name.gff"


    #if [[ -f "$output_file" && -f "$output_file1" ]]; then
        #echo "Skip $file, $output_file and $output_file1 already exist."
        #continue
    #fi

    # Prodigal
    echo "processing $file..."
    prodigal -i "$file" -c -m -f gff -a "$output_file" -o "$output_file1" -p meta
done


find "$root_folder" -type d -name "separated_proteins" | while read -r separated_folder; do

    for file in "$separated_folder"/*.faa; do

        base_name=$(basename "$file" .faa)

        output_folder=$(dirname "$separated_folder")/macsyfinder
        mkdir -p "$output_folder"

        macsyfinder_output="$output_folder/$base_name"

        #if [[ -d "$macsyfinder_output" && "$(ls -A "$macsyfinder_output" 2>/dev/null)" ]]; then
            #echo "Skip $file, MacsyFinder result $macsyfinder_output already exists"
            #continue
        #fi

        # macsyfinder
        echo "processing $file..."
        macsyfinder --db-type gembase --replicon-topology circular --models CASFinder all --accessory-weight 1 --exchangeable-weight 1 --coverage-profile 0.4 --redundancy-penalty 1 --sequence-db "$file" --out-dir "$macsyfinder_output"
    done
done

#python3 Cas_loci.py "."
# python3 "import sys; sys.path.append('$python_module_dir'); sys.argv = ['Cas_loci.py', '$root_folder']; exec(open('$cas_loci_script').read())"
python3 -c "import sys; sys.path.append('$python_module_dir'); from Cas_loci import process_best_solution; process_best_solution('$input_folder/ct_output')"
python3 -c "import sys; sys.path.append('$python_module_dir'); from Cas_loci import filter_complete_cas_systems; filter_complete_cas_systems('$input_folder/ct_output/cas')"

