#!/bin/bash
job_id=$RANDOM$RANDOM
echo "Job ID: $job_id"

#dirs and paths
input_dir="10_pdb_subset"
output_base_dir="fatcat_results"
output_dir="${output_base_dir}/alignments"
matrix_dir="${output_base_dir}/matrices"
log_dir="${output_base_dir}/logs"

mkdir -p "$output_dir"
mkdir -p "$matrix_dir"
mkdir -p "$log_dir"

log_file="${log_dir}/fatcat_timelog_ID$job_id.txt"
fatcat_score_matrix_file="${matrix_dir}/fatcat_score_matrix_ID$job_id.tsv"
pvalue_matrix_file="${matrix_dir}/pvalue_matrix_ID$job_id.tsv"

echo "Pairwise Alignment FATCAT time log" > "$log_file"
echo "Job ID: $job_id" >> "$log_file"
echo "Format: protein1_protein2  START_TIME  END_TIME  DURATION_SEC" >> "$log_file"

#get pdbs
files=("$input_dir"/*.pdb)
num_files=${#files[@]}
num_alignments=0
total_time=0

if [ $num_files -eq 0 ]; then
    echo "no pdbs in $input_dir :("
    exit 1
fi

echo "$num_files pdb files found"

#matrix
echo -n "Protein" > "$fatcat_score_matrix_file"
echo -n "Protein" > "$pvalue_matrix_file"
for ((i=0; i<num_files; i++)); do
    f=$(basename "${files[i]}" .pdb)
    echo -ne "\t$f" >> "$fatcat_score_matrix_file"
    echo -ne "\t$f" >> "$pvalue_matrix_file"
done
echo >> "$fatcat_score_matrix_file"
echo >> "$pvalue_matrix_file"

declare -A fatcat_score_matrix
declare -A pvalue_matrix

#loop over each pair and print progress
total_pairs=$((num_files * (num_files - 1) / 2))
current_pair=0

for ((i=0; i < num_files; i++)); do
    f1=$(basename "${files[i]}" .pdb)
    echo -n "$f1" >> "$fatcat_score_matrix_file"
    echo -n "$f1" >> "$pvalue_matrix_file"

    for ((j=0; j < num_files; j++)); do
        if [ $i -eq $j ]; then
            # Diagonal (self-comparison)
            fatcat_score_matrix["$f1,$f1"]="1000.00"
            pvalue_matrix["$f1,$f1"]="0.00"
            echo -ne "\t1000.00" >> "$fatcat_score_matrix_file"
            echo -ne "\t0.00" >> "$pvalue_matrix_file"
            continue
        elif [ $j -lt $i ]; then
            # Lower triangle - use values from upper triangle
            f2=$(basename "${files[j]}" .pdb)
            echo -ne "\t${fatcat_score_matrix["$f2,$f1"]}" >> "$fatcat_score_matrix_file"
            echo -ne "\t${pvalue_matrix["$f2,$f1"]}" >> "$pvalue_matrix_file"
            continue
        fi

        f2=$(basename "${files[j]}" .pdb)
        prefix="${f1}_${f2}"
        ((current_pair++))
        
        echo "Starting alignment $current_pair/$total_pairs: $prefix"
        start_time=$(date +%s)
        start_time_human=$(date '+%Y-%m-%d %H:%M:%S')

        #run fatcat
        FATCAT -p1 "${files[i]}" -p2 "${files[j]}" -o "${output_dir}/${prefix}" -m -ac -t > "${output_dir}/${prefix}.aln" 2>&1

        end_time=$(date +%s)
        end_time_human=$(date '+%Y-%m-%d %H:%M:%S')
        duration=$((end_time - start_time))
        
        echo "Finished alignment: $prefix in ${duration}s"

        echo -e "${prefix}\t${start_time_human}\t${end_time_human}\t${duration}" >> "$log_file"

        # Extract Score and P-value from output file
        #score=$(grep "Score" "${output_dir}/${prefix}.aln" | awk '{print $(NF)}')
        #pvalue=$(grep "P-value" "${output_dir}/${prefix}.aln" | awk '{print $NF}')
	#using pearl to get percies score
	score=$(grep -oP "Score \K[\d.]+" "${output_dir}/${prefix}.aln")
	pvalue=$(grep -oP "P-value \K[\d.e+-]+" "${output_dir}/${prefix}.aln")

        #save matrices
        fatcat_score_matrix["$f1,$f2"]=$score
        pvalue_matrix["$f1,$f2"]=$pvalue

        echo -ne "\t$score" >> "$fatcat_score_matrix_file"
        echo -ne "\t$pvalue" >> "$pvalue_matrix_file"

        ((num_alignments++))
        total_time=$((total_time + duration))
    done
    echo >> "$fatcat_score_matrix_file"
    echo >> "$pvalue_matrix_file"
done

#print where is everything
echo -e "\n###FATCAT Pairwise Comparison Summary###" | tee -a "$log_file"
echo "Job ID: $job_id" | tee -a "$log_file"
echo "Input directory: $input_dir" | tee -a "$log_file"
echo "Number of PDB files: $num_files" | tee -a "$log_file"
echo "Total alignments performed: $num_alignments" | tee -a "$log_file"
echo "Total computation time: $total_time seconds" | tee -a "$log_file"

if [[ $num_alignments -gt 0 ]]; then
  avg_time=$(echo "scale=2; $total_time / $num_alignments" | bc)
  echo "Average time per alignment: $avg_time seconds" | tee -a "$log_file"
fi

echo -e "\nResults saved in:" | tee -a "$log_file"
echo "1. Alignment files: $output_dir" | tee -a "$log_file"
echo "2. Score matrix: $fatcat_score_matrix_file" | tee -a "$log_file"
echo "3. P-value matrix: $pvalue_matrix_file" | tee -a "$log_file"
echo "4. Log file: $log_file" | tee -a "$log_file"
