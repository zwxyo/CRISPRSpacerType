#!/bin/bash


script_dir=$(cd "$(dirname "$0")" && pwd)
# python module path
python_module_dir="$script_dir/../src"

# the user's current working directory
user_current_dir=$(pwd)

# Default parameters
CSI="$user_current_dir"
result_summary="$CSI/ct_output"
CSO="$result_summary/CRISPR"
#CTI="$CSO"
#CTO="$current_dir/ct_utput/output_process.csv"
BT=""
BI="$CSI"
DB=""
CAS=false
CT=false
MLST=false
ONLY_CT=false
PCR_seq=false
# THREADS=1


set -m
trap 'echo -e "\nInterrupt signal captured, all processes are being terminated..."; kill -- -$$; exit' SIGINT


# Help Information
usage() {
    echo -e "\n\033[1;36m============================================================\033[0m"
    echo -e "\033[1;32m            Welcome to the CRISPRSpacerType Pipeline!          \033[0m"
    echo -e "\033[1;36m============================================================\033[0m"
    echo -e "\033[1;33mUsage:\033[0m $0 [options]"
    echo ""
    echo -e "\033[1;34mOptions:\033[0m"
    echo -e "  \033[1;32m-i, --csi <path>\033[0m    Input files are used to identify CRISPR (default: user current directory)"
    echo -e "  \033[1;32m-o, --cso <path>\033[0m    Output files of the results of CRISPR identification (default: user current directory/ct_output/CRISPR)"
    echo -e "  \033[1;32m--cas\033[0m           Identify cas proteins"
    # echo -e "  \033[1;32m--cti <path>\033[0m    CRISPRType input file (default: user_output)"
    # echo -e "  \033[1;32m--cto <path>\033[0m    CRISPRType output file (default: current directory/ct_output/output_process.csv)"
    echo -e "  \033[1;32m--bt <m|s>\033[0m      BLAST type: 'm' for custom database, 's' for standard"
    echo -e "  \033[1;32m--bi <path>\033[0m     BLAST query sequence (default: same as --csi)"
    echo -e "  \033[1;32m--db <path>\033[0m     BLAST database"
    echo -e "  \033[1;32m--ct\033[0m            CRISPR Typing Number"
    echo -e "  \033[1;32m--mlst\033[0m          MLST result"
    echo -e "  \033[1;32m--only\033[0m          Only perform CRISPR Typing"
    echo -e "  \033[1;32m--pcr\033[0m           Input pcr seq"
    echo -e "  \033[1;32m-h, --help\033[0m      Show this help message and exit"
    echo -e "\033[1;36m-------------------------------------------------------------\033[0m"
    echo -e "\033[1;34mFor more details, see:\033[0m \033[1;36mhttps://github.com/zwxyo/CRISPRSpacerType\033[0m"
    echo -e "\n"
    exit 0
}


while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i|--csi) CSI="$2"; shift 2;;
    -o|--cso) CSO="$2"; shift 2;;
    --cas) CAS=true; shift;;
#    --cti) CTI="$2"; shift 2;;
#    --cto) CTO="$2"; shift 2;;
    --bt) BT="$2"; shift 2;;
    --bi) BI="$2"; shift 2;;
    --db) DB="$2"; shift 2;;
    --ct) CT=true; CAS=true; shift;;
    --mlst) MLST=true; shift;;
    --only) ONLY_CT=true; shift;;
    --pcr) PCR_seq=true; shift;;
#    --threads)
#      if [[ "$2" == "auto" ]]; then
#        THREADS=$(nproc)
#      else
#        THREADS="$2"
#      fi
#      shift 2;;
    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

#echo "Using $THREADS threads"

#================================================================================
if [[ "$PCR_seq" == true && "$CT" == true ]]; then

  bash "$script_dir/CRISPR.sh" "$CSI" "$CSO"

  if [[ $? -eq 0 ]]; then
    echo "Running PCR sequence processing..."
    CTO="$(dirname "$CSO")"

    python3 -c "import sys; sys.path.append('$python_module_dir'); from PCR_seq import PCR_process; PCR_process('$CSO', '${CTO}/result_pcr.csv')"
    python3 -c "import sys; sys.path.append('$python_module_dir'); from PCR_seq import pre_typing_process; pre_typing_process('${CTO}/result_pcr.csv', '${result_summary}/CRISPR_sort/Cronobacter_Genus')"

    python3 -c "import sys; sys.path.append('$python_module_dir'); from Spacer_ct_numbering import summary_pcr; summary_pcr('${result_summary}/CRISPR_sort', '${result_summary}/spacer/spacer_order', '${result_summary}/spacer/serial', '${result_summary}/spacer/ct')"

  else
    echo -e "\033[1;31mError:\033[0m CRISPR recognition fails, skip the next step..."
    exit 1
  fi

  exit 0

#================================================================================
# CT
elif [[ "$PCR_seq" == false && "$CT" == true && "$ONLY_CT" == false ]]; then

  # CAS identification
  if [[ "$CAS" == true ]]; then
      bash "$script_dir/prodigal_macsyfinder.sh" "$CSI"
  fi

  # CRISPR identification
  # exec bash "$script_dir/CRISPR_final.sh" "$CSI" "$CSO"
  bash "$script_dir/CRISPR.sh" "$CSI" "$CSO"

  #if [[ "$CAS" == true ]]; then
  #    bash "$script_dir/prodigal_macsyfinder.sh" "$CSI"
  #fi


  if [[ $? -eq 0 ]]; then
      # CRISPRType process data
      CTO="$(dirname "$CSO")"
      python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Collect_results; Collect_results('$CSO', '${CTO}/result.csv')"
      python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Data_filtering; Data_filtering('${CTO}/result.csv', '${CTO}/result_process.csv')"
  else
      echo -e "\033[1;31mError:\033[0m CRISPR recognition fails, skip the next step..."
      exit 1
  fi


  #CTO="$(dirname "$CSO")"
  #python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Collect_results; Collect_results('$CSO', '$CTO/result.csv')"
  #python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Data_filtering; Data_filtering('${CTO}/result.csv', '${CTO}/result_process.csv')"

  #--------------------------------------------------------------------------------
  # BLAST
  # check parameter BI
  if [[ "$BT" != "" && "$BT" != "m" && "$BT" != "s" ]]; then
      echo -e "\033[1;31mError:\033[0m --bt must be 'm' or 's'."
      exit 1
  elif [[ "$BT" == "s" && -n "$BI" && "$BI" != "$CSI" ]]; then
      echo -e "\033[1;31mError:\033[0m When using the standard BLAST type (--bt s), the --bi option must either be omitted or set to the same path as --csi."
      exit 1
  elif [[ "$BT" == "m" && -z "$DB" ]]; then
      echo -e "\033[1;31mError:\033[0m If you opt to utilize a custom BLAST database, please ensure the presence of your custom database FNA file."
      exit 1
  elif [[ "$BT" == "m" && ! -f "$DB" ]]; then
      echo -e "\033[1;31mError:\033[0m Custom BLAST DB file $DB not found!"
      exit 1
  elif [[ "$BT" == "s" && -n "$DB" ]]; then
      echo -e "\033[1;31mError:\033[0m Customization of the BLAST database path is not required."
      exit 1
  fi

  if compgen -G "$BI/*.fna" > /dev/null; then
      if [[ "$BT" == "m" && -n "$DB" ]]; then
        # python3 BLAAST/BLAST_M.py "$BI"
        # python3 -c "import sys; sys.path.append('$python_module_dir'); sys.argv = ['BLAST_M.py', '$BI', '$DB', '$result_summary/blast_result']; exec(open('$BLAST_M_script').read())"
        python3 -c "import sys; sys.path.append('$python_module_dir'); from BLAST_M import BLAST_M; BLAST_M('$BI', '$DB', '${result_summary}/blast_result')"

      elif [[ "$BT" == "s" ]]; then
        # python3 BLAST/BLAST_S.py "$BI"
        # python3 -c "import sys; sys.path.append('$python_module_dir'); sys.argv = ['BLAST_S.py', '$BI', '$result_summary/blast_result']; exec(open('$BLAST_S_script').read())"
        python3 -c "import sys; sys.path.append('$python_module_dir'); from BLAST_S import BLAST_S; BLAST_S('$BI', '${result_summary}/blast_result')"

      fi
  else
      echo -e "\033[1;31mError:\033[0m the '.fna' file was not found in directory $BI, skipping the BLAST step..."
      exit 1
  fi

  #--------------------------------------------------------------------------------
  # CRISPR classification
  # python3 CRISPR_classification.py
  # python3 -c "import sys; sys.path.append('$python_module_dir'); sys.argv = ['CRISPR_classification.py', '${CTO}/result_process.csv', '$result_summary/cas', '$result_summary/CRISPR_sort', '$result_summary/blast_result']; exec(open('$CRISPR_classification').read())"
  python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_classification import CRISPR_classification; CRISPR_classification('${CTO}/result_process.csv', '${result_summary}/cas', '${result_summary}/CRISPR_sort', '${result_summary}/blast_result')"

  # python3 -c "import sys; sys.path.append('$python_module_dir'); from Spacer_ct_numbering import summary; summary()"
  # python3 -c "import sys; sys.path.append('$python_module_dir'); sys.argv = ['Spacer_ct_numbering.py', '$result_summary/CRISPR_sort', '$result_summary/spacer/spacer_order', '$result_summary/spacer/serial', '$result_summary/spacer/ct']; exec(open('$CT_number').read())"
  python3 -c "import sys; sys.path.append('$python_module_dir'); from Spacer_ct_numbering import summary; summary('${result_summary}/CRISPR_sort', '${result_summary}/spacer/spacer_order', '${result_summary}/spacer/serial', '${result_summary}/spacer/ct')"

  #--------------------------------------------------------------------------------
  # mlst
  output_mlst="$CTO/mlst"

  if [[ "$MLST" == true ]]; then
      if ! command -v mlst &> /dev/null; then
        echo -e "\033[1;31mError:\033[0m 'mlst' command not found. Please install mlst first."
        exit 1
      else
        mkdir -p "$output_mlst"
        # mlst --scheme cronobacter --legacy --csv *.fna > mlst_results.csv
        mlst --scheme cronobacter --legacy --csv ${CSI}/*.fna > ${output_mlst}/mlst_results.csv
      fi
  fi

#================================================================================
elif [[ "$CT" == true && "$ONLY_CT" == true ]]; then
  python3 -c "import sys; sys.path.append('$python_module_dir'); from Spacer_ct_numbering import summary; summary('${result_summary}/CRISPR_sort', '${result_summary}/spacer/spacer_order', '${result_summary}/spacer/serial', '${result_summary}/spacer/ct')"

  #--------------------------------------------------------------------------------
  if [[ "$MLST" == true ]]; then
      if ! command -v mlst &> /dev/null; then
        echo -e "\033[1;31mError:\033[0m 'mlst' command not found. Please install mlst first."
        exit 1
      else
        mkdir -p "$output_mlst"
        # mlst --scheme cronobacter --legacy --csv *.fna > mlst_results.csv
        mlst --scheme cronobacter --legacy --csv ${CSI}/*.fna > ${output_mlst}/mlst_results.csv
      fi
  fi

#================================================================================
elif [[ "$CT" == false && "$ONLY_CT" == true ]]; then
  echo -e "\033[1;31mError:\033[0m --only must be used together with --ct."
  exit 1

#================================================================================
else
  # CAS identification
  if [[ "$CAS" == true ]]; then
      bash "$script_dir/prodigal_macsyfinder.sh" "$CSI"
  fi

  # CRISPR identification
  # exec bash "$script_dir/CRISPR_final.sh" "$CSI" "$CSO"
  bash "$script_dir/CRISPR.sh" "$CSI" "$CSO"

  #if [[ "$CAS" == true ]]; then
  #    bash "$script_dir/prodigal_macsyfinder.sh" "$CSI"
  #fi


  if [[ $? -eq 0 ]]; then
      # CRISPRType process data
      CTO="$(dirname "$CSO")"
      python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Collect_results; Collect_results('$CSO', '${CTO}/result.csv')"
      python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Data_filtering; Data_filtering('${CTO}/result.csv', '${CTO}/result_process.csv')"
  else
      echo -e "\033[1;31mError:\033[0m CRISPR recognition fails, skip the next step..."
      exit 1
  fi


  #CTO="$(dirname "$CSO")"
  #python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Collect_results; Collect_results('$CSO', '$CTO/result.csv')"
  #python3 -c "import sys; sys.path.append('$python_module_dir'); from CRISPR_process import Data_filtering; Data_filtering('${CTO}/result.csv', '${CTO}/result_process.csv')"

  #--------------------------------------------------------------------------------
  # BLAST
  # check parameter BI
  if [[ "$BT" != "" && "$BT" != "m" && "$BT" != "s" ]]; then
      echo -e "\033[1;31mError:\033[0m --bt must be 'm' or 's'."
      exit 1
  elif [[ "$BT" == "s" && -n "$BI" && "$BI" != "$CSI" ]]; then
      echo -e "\033[1;31mError:\033[0m When using the standard BLAST type (--bt s), the --bi option must either be omitted or set to the same path as --csi."
      exit 1
  elif [[ "$BT" == "m" && -z "$DB" ]]; then
      echo -e "\033[1;31mError:\033[0m If you opt to utilize a custom BLAST database, please ensure the presence of your custom database FNA file."
      exit 1
  elif [[ "$BT" == "m" && ! -f "$DB" ]]; then
      echo -e "\033[1;31mError:\033[0m Custom BLAST DB file $DB not found!"
      exit 1
  elif [[ "$BT" == "s" && -n "$DB" ]]; then
      echo -e "\033[1;31mError:\033[0m Customization of the BLAST database path is not required."
      exit 1
  fi

  if compgen -G "$BI/*.fna" > /dev/null; then
      if [[ "$BT" == "m" && -n "$DB" ]]; then
        # python3 BLAAST/BLAST_M.py "$BI"
        # python3 -c "import sys; sys.path.append('$python_module_dir'); sys.argv = ['BLAST_M.py', '$BI', '$DB', '$result_summary/blast_result']; exec(open('$BLAST_M_script').read())"
        python3 -c "import sys; sys.path.append('$python_module_dir'); from BLAST_M import BLAST_M; BLAST_M('$BI', '$DB', '${result_summary}/blast_result')"

      elif [[ "$BT" == "s" ]]; then
        # python3 BLAST/BLAST_S.py "$BI"
        # python3 -c "import sys; sys.path.append('$python_module_dir'); sys.argv = ['BLAST_S.py', '$BI', '$result_summary/blast_result']; exec(open('$BLAST_S_script').read())"
        python3 -c "import sys; sys.path.append('$python_module_dir'); from BLAST_S import BLAST_S; BLAST_S('$BI', '${result_summary}/blast_result')"

      fi
  else
      echo -e "\033[1;31mError:\033[0m the '.fna' file was not found in directory $BI, skipping the BLAST step..."
      exit 1
  fi

  #--------------------------------------------------------------------------------
  # mlst
  output_mlst="$CTO/mlst"

  if [[ "$MLST" == true ]]; then
      if ! command -v mlst &> /dev/null; then
        echo -e "\033[1;31mError:\033[0m 'mlst' command not found. Please install mlst first."
        exit 1
      else
        mkdir -p "$output_mlst"
        # mlst --scheme cronobacter --legacy --csv *.fna > mlst_results.csv
        mlst --scheme cronobacter --legacy --csv ${CSI}/*.fna > ${output_mlst}/mlst_results.csv
      fi
  fi
fi
#================================================================================
echo -e "\n\033[1;32mAll steps completed successfully!\033[0m"
echo -e "Output directory: \033[1;36m${result_summary}\033[0m"
