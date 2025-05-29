#!/bin/bash

mkdir -p alphafold_structures

while read -r uniprot_id; do
  echo "downloading af structure for $uniprot_id"

  url="https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-model_v4.pdb"
  output="alphafold_structures/${uniprot_id}.pdb"

  if wget -q --spider "$url"; then
    wget -q "$url" -O "$output"
    echo " saved as $output"
  else
    echo "$uniprot_id download failed"
  fi

done < uniprot_ids.txt

