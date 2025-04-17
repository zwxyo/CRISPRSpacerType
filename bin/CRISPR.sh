#!/bin/bash


# Input directory and output directory
#INPUT_DIR=${1:-"user_input"}
#OUTPUT_DIR=${2:-"user_output"}
INPUT_DIR="$1"
OUTPUT_DIR="$2"


script_dir=$(cd "$(dirname "$0")" && pwd)
crisprcasstack_dir="$script_dir/../external"
ccs_script="$crisprcasstack_dir/CRISPRCasStack.py"

colorize() {
    echo -e "\033[$1m$2\033[0m"
}

show_progress() {
    local pid=$1
    # local spin=('⠁' '⠂' '⠄' '⡀' '⢀' '⠠' '⠐' '⠈' '⠁' '⠂' '⠄' '⡀' '⢀' '⠠' '⠐' '⠈')
    local delay=0.1
    # local i=0
    local progress_started=false

    while kill -0 "$pid" 2>/dev/null; do
        if [[ "$progress_started" = false ]]; then
            progress_started=true
            # printf "\r%s %s" "$(colorize 36 'processing:')" "$(colorize 32 "${spin[$i]}")"
            # printf "\r%s %s[%s%s] %s" "$(colorize $color "${spin[$i]}")" "$(colorize 36 'processing:')"
            printf "\r%s" "$(colorize 36 'processing:')"
        fi
        # ((i = (i + 1) % ${#spin[@]}))

        sleep "$delay"
    done
    printf "\r%s %s\n" "$(colorize 32 '✔')" "$(colorize 36 'processed')"
}


# FASTA format validation
validate_fna() {
    local file="$1"
    local valid=true
    local in_header=0
    local has_sequence=0
    local line_number=0
    local current_header_line=0
    local file_has_header=0

    while IFS= read -r line; do
        ((line_number++))
        # header
        if [[ "$line" =~ ^\> ]]; then
            file_has_header=1

            if ((in_header == 1)) && ((has_sequence == 0)); then
                echo "error: there is no sequence on line ${current_header_line}, a new header appears on line ${line_number}" >&2
                valid=false
            fi

            in_header=1
            current_header_line=$line_number
            has_sequence=0
        else
            # after header
            if ((in_header == 1)); then

                clean_line=$(echo "$line" | tr -d -c "ACGTN")
                if [[ -z "$clean_line" ]]; then
                    echo "error: line ${line_number} (blank line or invalid character) appears after the header (line ${current_header_line})" >&2
                    valid=false
                else
                    has_sequence=1
                    if [[ ! "$clean_line" =~ ^[ACGTN]+$ ]]; then
                        echo "error: line ${line_number} contains invalid character:'$line'" >&2
                        valid=false
                    fi
                fi
            else
                if [[ -z "$line" ]]; then
                    echo "error: line ${line_number} is blank" >&2
                    valid=false
                else
                    echo "error: line ${line_number} is not after any header: '$line'" >&2
                    valid=false
                fi
            fi
        fi
    done < "$file"


    if ((in_header == 1 && has_sequence == 0)); then
        echo "error: the last header (line ${current_header_line}) is missing a sequence" >&2
        valid=false
    fi


    if ((file_has_header == 0)); then
        echo "error: the file does not contain any header lines that start with '>'" >&2
        valid=false
    fi

    [ "$valid" = true ]
}


if ! find "$INPUT_DIR" -type f -iname "*.fna" | grep -q .; then
    echo "No '.fna' files were found, and the script terminated"
    exit 1
fi


rm -rf "$crisprcasstack_dir/result"/*


# find "$INPUT_DIR" -type f \( -iname "*.fna" \) | while read -r file; do
find "$INPUT_DIR" -maxdepth 1 -type f \( -iname "*.fna" \) | while read -r file; do

    # If file exists
    if [[ -f "$file" ]]; then
        # Get the filename and its parent folder name
        filename=$(basename "$file")
        name="${filename%.*}"

        relative_path="${file#$INPUT_DIR/}"
        subdir_parent=$(dirname "$relative_path")


        if [[ "$subdir_parent" == "." ]]; then
            subdir_parent=""
        fi

        if [[ -z "$subdir_parent" ]]; then
            # subdir="$OUTPUT_DIR/$name"
            subdir="$crisprcasstack_dir/result/$name"
        else
            # subdir="$OUTPUT_DIR/$subdir_parent/$name"
            subdir="$crisprcasstack_dir/result/$name"
        fi

        mkdir -p "$subdir"

        # Check if output for this file already exists
        check_file="$subdir/check.txt"
        error_log="$subdir/validation_errors.log"


        if [[ ! -f "$check_file" ]]; then
            echo -e "$(colorize 34 "Documents are being verified:") $file"

            if ! validate_fna "$file" 2> "$error_log"; then
                echo -e "$(colorize 31 "File format error, skip processing. See $error_log for details")"
                continue
            fi

            rm -rf "$subdir"/*

            echo "$file..."

            # python3 "$crisprcasstack_dir/CRISPRCasStack.py" -i "$file" -o "$subdir" -m g
            # python3 -c "import sys; sys.path.append('$crisprcasstack_dir'); sys.argv = ['CRISPRCasStack.py', '-i', '$file', '-o', '$subdir', '-m', 'g']; exec(open('$ccs_script').read())"
            # python3 -c "import sys, os; os.chdir('$crisprcasstack_dir'); sys.path.append('$crisprcasstack_dir'); sys.argv = ['CRISPRCasStack.py', '-i', '$file', '-o', '$subdir', '-m', 'g']; exec(open('$ccs_script').read())"

            run_analysis() {
                python3 -c "import sys, os; os.chdir('$crisprcasstack_dir'); sys.path.append('$crisprcasstack_dir'); sys.argv = ['CRISPRCasStack.py', '-i', '$file', '-o', '$subdir', '-m', 'g']; exec(open('$ccs_script').read())"
            }

            if run_analysis & then
                analysis_pid=$!
                show_progress "$analysis_pid"
                wait "$analysis_pid"
                if [[ $? -eq 0 ]]; then
                    echo "Processing completed for $file" > "$check_file"
                    mkdir -p "$OUTPUT_DIR"

                    rsync -av --delete "$crisprcasstack_dir/result"/* "$OUTPUT_DIR"/
                else
                    echo -e "$(colorize 31 'CRISPR script execution failed')"
                fi
            fi

        else

             echo -e "$(colorize 33 "Skip: $file, a result already exists")"
        fi
    fi
done
