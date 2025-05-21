#!/bin/bash
job_id=$RANDOM$RANDOM
echo "Job ID: $job_id"

input_dir="pdb_subset_ID$job_id"
output_dir="fatcat_10_results_ID$job_id"
log_file="fatcat_timelog_ID$job_id.txt"

mkdir -p "$output_dir"
echo "Pairwise Alignment FATCAT time log" > "$log_file"
echo "Job ID: $job_id" > "$log_file"
echo "Format: protein1_protein2  START_TIME  END_TIME  DURATION_SEC" >> "$log_file"

files=("$input_dir"/*.pdb)
num_alignments=0
total_time=0

for ((i=0; i < ${#files[@]}; i++)); do
  for ((j=i+1; j < ${#files[@]}; j++)); do
    f1=$(basename "${files[i]}" .pdb)
    f2=$(basename "${files[j]}" .pdb)
    prefix="${f1}_${f2}"

    echo "Starting alignment: $prefix"
    start_time=$(date +%s)
    start_time_human=$(date '+%Y-%m-%d %H:%M:%S')

    FATCAT -p1 "${files[i]}" -p2 "${files[j]}" -o "${output_dir}/${prefix}" -m -ac -t

    end_time=$(date +%s)
    end_time_human=$(date '+%Y-%m-%d %H:%M:%S')

    duration=$((end_time - start_time))
    echo "Finished alignment: $prefix in ${duration}s"

    echo -e "${prefix}\t${start_time_human}\t${end_time_human}\t${duration}" >> "$log_file"

    ((num_alignments++))
    total_time=$((total_time + duration))
  done
done

echo -e "\nTotal alignments: $num_alignments" | tee -a "$log_file"
echo "Total time (seconds): $total_time" | tee -a "$log_file"

if [[ $num_alignments -gt 0 ]]; then
  avg_time=$(echo "scale=2; $total_time / $num_alignments" | bc)
else
  avg_time=0
fi

echo "Average time per alignment (seconds): $avg_time" | tee -a "$log_file"

