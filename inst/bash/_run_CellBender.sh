#!/bin/bash

# Description: Run CellBender sequentially for multiple samples
# Author: David Rodriguez Morales
# Date created: 21-12-23
# Date Modified: 21-05-25
# Version: 1.0



function usage()
{
	echo
	echo "Run CellBender for multiple samples (adapted for CellRanger)"
	echo
	echo "Usage: $0 [-i name of samples (comma separated)] [-o output path] [--cellRanger-output path to CellRanger outputs]"
	echo
	echo "OPTIONS"
	echo " [-i | --inputNames SampleName1,SampleName2,SampleName3] - Sample Names (comma separated)"
	echo " [-o | --outputPath ] - Path to output folder for CellBender. CellBender output files will be saved under this folder. The prefix SampleName{N}_out will be use"
	echo " [--cellRanger-output] - Path to folder containing CellRanger output of multiple samples (subfolder for each sample). SampleName{N} should be the same as the subfolder containing the results of the sample"
	echo " [--cpu-threads] - Number of CPU threads (Default 15)"
	echo " [--epochs] - Number of training epochs (Default 150)"
	echo " [--lr] - Learning rate (Default 1e-5)"
	echo " [--cuda] - if added, the script will run on GPU"
	echo " [--expected-cells] - Expected number of cells. Can be extracted from cellranger report"
	echo " [--total-droplets] - Total number of cells to use to infer cell probabilty. Cells above are surely empty"
	echo " [--force-umi-prior] - Ignore heuristics from cellbender and set the prior of minimum number of counts for a cell (Default 500)"
	echo " [--estimator-multiple-cpu] - if added, the script will run estimator step using multiple CPU, Not recommended for big datasets!"
	echo " [--log] - If added, a log file will be generated for each sample"
	echo
	echo " EXAMPLE"
	echo " With GPU: bash ./run_cellbender -i sample1,sample2,sample3 -o ./cellbender --cellRanger-output ./cellranger --cuda  --log"
	echo " Without GPU: bash ./run_cellbender -i sample1,sample2,sample3 -o ./cellbender --cellRanger-output ./cellranger  --log"
	exit 1
}

# Normalise path, if user provide path/ transform to path
function normalize_path() {
    echo "$1" | sed 's:/*$::'
}

# Default values
epochs=150
cuda_arg=""
estimator_multi_arg=""
cpu_threads=15
lr_arg=0.00001
log_flag=""
expected_cells=""
total_droplets=""
force_umi_prior=500


# Parameters
while test $# -gt 0
do
	case $1 in
		-h | --help)
		usage
		;;
		-i | --inputNames)
		IFS="," read -r -a SAMPLENAMES <<< "$2"
		shift 2
		;;
		-o | --outputPath)
		outPath=$(normalize_path "$2")
		shift 2
		;;
		--cellRanger-output)
		cellRanger_out=$(normalize_path "$2")
		shift 2
		;;
		--cpu-threads)
		cpu_threads=$2
		shift 2
		;;
		--epochs)
		epochs=$2
		shift 2
		;;
    --lr)
    lr_arg=$2
    shift 2
    ;;
    --expected-cells)
    expected_cells=$2
    shift 2
    ;;
    --total-droplets)
    total_droplets=$2
    shift 2
    ;;
    --force-umi-prior)
    force_umi_prior=$2
    shift 2
    ;;
    --cuda)
    cuda_arg="--cuda"
    shift
    ;;
    --estimator-multiple-cpu)
    estimator_multi_arg="--estimator-multiple-cpu"
    shift
    ;;
    --log)
    log_flag="true"
    shift
    ;;
		*) echo "ERROR: unknown argument $1"; exit 1
	esac
done

# Check if mandatory parameters are set
if [[ -z "${SAMPLENAMES[*]}" || -z "$outPath" || -z "$cellRanger_out" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Print raw commands to the log file
command_file="${outPath}/commands_CellBender.txt"
echo "# Raw commands" >> "$command_file"
join_names=$(IFS=, ; echo "${SAMPLENAMES[*]}")
echo "$0 -i $join_names -o ${outPath} --cellRanger-output ${cellRanger_out} --cpu-threads ${cpu_threads} --epochs ${epochs} --learning-rate ${lr_arg} --expected-cells ${expected_cells} --total-droplets ${total_droplets} --force-umi-prior ${force_umi_prior}  ${cuda_arg} ${estimator_multi_arg}  ${log_flag:+--log}"  >> "$command_file"

i=0
for SAMPLENAME in ${SAMPLENAMES[*]}; do
  # Create a log file for each sample
  log_file="${outPath}/${SAMPLENAME}_cellbender.log"

	# Run CellBender
	echo >> "$command_file"
	echo "Run CellBender for ${SAMPLENAME}" >> "$command_file"

  start=$(date +%s)

  if [[ -n "$estimator_multi_arg" ]]; then
    echo 'WARNING!!!! --estimator-multiple-cpu is set, this is not recommended for big datasets'
  else
    echo ''
  fi

    input_file="$(find "${cellRanger_out}/${SAMPLENAME}/outs/"*raw_feature_bc_matrix.h5 | head -n 1)"

	cellbender_cmd="cellbender remove-background \
	--input=${input_file} \
	--output=${outPath}/${SAMPLENAME}_out.h5 \
	--cpu-threads=${cpu_threads} \
	--epochs=${epochs} \
	--learning-rate=${lr_arg} \
	--force-cell-umi-prior=${force_umi_prior} \
	${estimator_multi_arg} \
	${cuda_arg}"

  # Conditionally add --expected-cells if value is set
  if [[ -n "$expected_cells" ]]; then
    cellbender_cmd+=" --expected-cells=${expected_cells}"
  fi

  # Conditionally add --total-droplets-included if value is set
  if [[ -n "$total_droplets" ]]; then
    cellbender_cmd+=" --total-droplets-included=${total_droplets}"
  fi

	echo "${cellbender_cmd}" >> "$command_file"

  # Generate a log file depending on the flag
  if [[ -n "$log_flag" ]]; then
    eval "$cellbender_cmd" 2>&1 | tee -a "${log_file}"
        # cellbender remove-background --input="${cellRanger_out}/${SAMPLENAME}/outs/raw_feature_bc_matrix.h5" --output="${outPath}/${SAMPLENAME}_out.h5" --cpu-threads="${cpu_threads}"  --epochs="${epochs}" --learning-rate="${lr_arg}" --expected-cells="${expected_cells}"  --total-droplets-included="${total_droplets}"  --force-empty-umi-prior="${force_umi_prior}" ${cuda_arg} ${estimator_multi_arg} 2>&1 | tee -a "${log_file}"
  else
    eval "$cellbender_cmd"
        # cellbender remove-background --input="${cellRanger_out}/${SAMPLENAME}/outs/raw_feature_bc_matrix.h5" --output="${outPath}/${SAMPLENAME}_out.h5" --cpu-threads="${cpu_threads}"  --epochs="${epochs}" --learning-rate="${lr_arg}" --expected-cells="${expected_cells}"  --total-droplets-included="${total_droplets}"   --force-empty-umi-prior="${force_umi_prior}" ${cuda_arg} ${estimator_multi_arg}
  fi

  # If it fails, it continues to the next sample
  if [ $? -ne 0 ]; then
        echo "Error: CellBender command failed for ${SAMPLENAME}"
        continue
  fi

  echo "Duration: $((($(date +%s)-$start)/60)) min" >> "$command_file"

	let i++

done
