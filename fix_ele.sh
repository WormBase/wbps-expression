comm -1 -3 elegans_studies_to_keep.tsv elegans_all_studies.tsv  | while read -r study; do
 if [[ -f "studies/caenorhabditis_elegans/${study}/${study}.skipped_runs.tsv" ]] && [[ -f "studies/caenorhabditis_elegans/${study}/${study}.design.tsv" ]]; then
 	cut -f 1 studies/caenorhabditis_elegans/${study}/${study}.design.tsv | tail -n +2 >> studies/caenorhabditis_elegans/${study}/${study}.skipped_runs.tsv
 else 
 	cut -f 1 studies/caenorhabditis_elegans/${study}/${study}.design.tsv > studies/caenorhabditis_elegans/${study}/${study}.skipped_runs.tsv
 fi
 rm studies/caenorhabditis_elegans/${study}/${study}.design.tsv
done
