#!/bin/bash

# Paths
RDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/code"
RESDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/results/ANN"

# Navigate to results directory
cd "$RESDIR"

# Loop over alleles
for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701 
do

cd $a.res

for sc in blosum sparse
do

cd $sc 

for nh in 5
do

cd hidden.$nh

for epi in 0.05
do

cd epi.$epi

eval_file="concat_performace.eval"

# If output doesn't exist, run the evaluation
if [ ! -f "$eval_file" ]; then
    result=$(python "$RDIR/AP_AUC_MCC.py" -pf concat.eval)
    echo "ANN,$a,$result,$sc,$nh,$epi" >> ../../../../"$eval_file"
fi


cd ..  
done
cd ..    
done
cd ..      
done
cd ..        
done
