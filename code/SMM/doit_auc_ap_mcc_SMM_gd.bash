#!/bin/bash

# Paths
RDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/code"
RESDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/results/SMM"

# Output file
output_file="$RESDIR/performance_metrics.csv"
echo "method,allele,AUC,AP,MCC,encoding,hid_neu,epi" > "$output_file"

# Navigate to results directory
cd "$RESDIR"/gradient_decent

# Loop over alleles
for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701 
do

cd $a.res

# Here you can type the lambdas to test
# 0 0.02 0.1 0.001
for l in 0.1   # Just choose the best one (since it's fast we do all of them)
do

#mkdir -p l.$l

cd l.$l

eval_file="concat_performace.eval"

# If output doesn't exist, run the evaluation
if [ ! -f "$eval_file" ]; then
    result=$(python "$RDIR/AP_AUC_MCC.py" -pf concat.eval)
    echo "SMM,$a,$result,GD,$l" >> ../../"$eval_file"
fi

# Append result to performance_metrics.csv
if [ -f "$eval_file" ]; then
    cat "$eval_file" >> "$RESDIR/performance_metrics.csv"
fi

cd ..  
done
cd ..    
done

