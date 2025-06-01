#!/bin/bash
job_id=$RANDOM$RANDOM
echo "Job ID: $job_id"
#generates results.csv that is saved in current dir, it wil be overwritten if run again (temporary solution for now)
#

#dirs and paths
input_dir="${1:-10_pdb_subset}" #run ./script.sh input_dir, if input_dir not provided it will use 10_pdb_subset
output_base_dir="fatcat_results"
output_dir="${output_base_dir}/alignments"
matrix_dir="${output_base_dir}/matrices"
csv_dir="${output_base_dir}/csv"
log_dir="${output_base_dir}/logs"
tmp_dir="${output_base_dir}/tmp"

mkdir -p "$output_dir"
mkdir -p "$matrix_dir"
mkdir -p "$csv_dir"
mkdir -p "$log_dir"
mkdir -p "$tmp_dir"

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
            fatcat_score_matrix["$f1,$f1"]="Nan"
            pvalue_matrix["$f1,$f1"]="Nan"
            echo -ne "\tNan" >> "$fatcat_score_matrix_file"
            echo -ne "\tNan" >> "$pvalue_matrix_file"
            continue
        elif [ $j -lt $i ]; then
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

        #run fatcat and extract scores
        FATCAT -p1 "${files[i]}" -p2 "${files[j]}" -o "${output_dir}/${prefix}" -m -ac -t > "${output_dir}/${prefix}.aln" 2>&1
        score=$(grep -oP "Score \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
        pvalue=$(grep -oP "P-value \K[\d.e+-]+" "${output_dir}/${prefix}.aln" || echo "Nan")

	#other metrics for csv
	twists=$(grep -oP "Twists \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	ini_len=$(grep -oP "ini-len \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	init_rmsd=$(grep -oP "ini-rmsd \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	opt_equ=$(grep -oP "opt-equ \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	opt_rmsd=$(grep -oP "opt-rmsd \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	chain_rmsd=$(grep -oP "chain-rmsd \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	align_len=$(grep -oP "align-len \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	gaps=$(grep -oP "gaps \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	afp_num=$(grep -oP "Afp-num \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	identity=$(grep -oP "Identity \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")
	similarity=$(grep -oP "Similarity \K[\d.]+" "${output_dir}/${prefix}.aln" || echo "Nan")

        #get num of atoms and residues in each prot
        atoms1=$(grep '^ATOM' "${files[i]}" | wc -l)
        residues1=$(grep '^ATOM' "${files[i]}" | grep ' CA ' | wc -l)
        atoms2=$(grep '^ATOM' "${files[j]}" | wc -l)
        residues2=$(grep '^ATOM' "${files[j]}" | grep ' CA ' | wc -l)
        #{
        #echo "Protein 1: $f1, Atoms: $atoms1, Residues: $residues1"
        #echo "Protein 2: $f2, Atoms: $atoms2, Residues: $residues2"
        #echo "fatcat score: $score"
        #} | tee -a results.txt
#csv
	csv_file="${csv_dir}/results$job_id.csv"
	if [ ! -f "$csv_file" ]; then
		#echo "Protein1,Atoms1,Residues1,Protein2,Atoms2,Residues2,FATCAT_score" > "$csv_file"
echo "4. TM-score matrix: $tmscore_matrix_file" | tee -a "$log_file"
		echo "Protein1,Atoms1,Residues1,Protein2,Atoms2,Residues2,FATCAT_score,P-value,twists,ini_len,init_rmsd,opt_equ,opt_rmsd,chain_rmsd,align_len,gaps,afp_num,identity,similarity" > "$csv_file"
	fi
	#echo "$f1,$atoms1,$residues1,$f2,$atoms2,$residues2,$score" >> "$csv_file"
	echo "$f1,$atoms1,$residues1,$f2,$atoms2,$residues2,$score,$pvalue,$twists,$ini_len,$init_rmsd,$opt_equ,$opt_rmsd,$chain_rmsd,$align_len,$gaps,$afp_num,$identity,$similarity" >> "$csv_file"

        end_time=$(date +%s)
        end_time_human=$(date '+%Y-%m-%d %H:%M:%S')
        duration=$((end_time - start_time))
        echo "Finished alignment: $prefix in ${duration}s"

        echo -e "${prefix}\t${start_time_human}\t${end_time_human}\t${duration}" >> "$log_file"

        #save matrices
        fatcat_score_matrix["$f1,$f2"]=$score
        pvalue_matrix["$f1,$f2"]=$pvalue

        echo -ne "\t$score" >> "$fatcat_score_matrix_file"
        echo -ne "\t$pvalue" >> "$pvalue_matrix_file"

        ((num_alignments++))
        total_time=$((total_time + duration))

        # Clean up foldseek temporary files
        rm -rf "$tmp_output_dir"
    done
    echo >> "$fatcat_score_matrix_file"
    echo >> "$pvalue_matrix_file"
done

#summary
echo -e "\nFATCAT Pairwise Comparison Summary" | tee -a "$log_file"
echo "job ID: $job_id" | tee -a "$log_file"
echo "input dir: $input_dir" | tee -a "$log_file"
echo "num  of PDB files: $num_files" | tee -a "$log_file"
echo "num of alignments: $num_alignments" | tee -a "$log_file"
echo "computation time: $total_time seconds" | tee -a "$log_file"

if [[ $num_alignments -gt 0 ]]; then
  avg_time=$(echo "scale=2; $total_time / $num_alignments" | bc)
  echo "avrg time per alignment: $avg_time s" | tee -a "$log_file"
fi

echo -e "\nresults saved in:" | tee -a "$log_file"
echo "1. alignment files: $output_dir" | tee -a "$log_file"
echo "2. score matrix: $fatcat_score_matrix_file" | tee -a "$log_file"
echo "3. p-value matrix: $pvalue_matrix_file" | tee -a "$log_file"
echo "4. csv file: $csv_file" | tee -a "$log_file"
echo "5. log file: $log_file" | tee -a "$log_file"

rm -rf "$tmp_dir"

