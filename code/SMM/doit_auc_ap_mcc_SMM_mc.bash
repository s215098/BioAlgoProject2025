# #!/bin/bash

# # Paths
# RDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/code"
# RESDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/results/SMM"

# # Navigate to results directory
# cd "$RESDIR/monte_carlo"

# # Loop over alleles
# for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701 
# do


# for l in 0.1   # Just choose the best one
# do

# # 10 20 50 100
# for t_steps in 20 100    # Just choose the best one
# do

# #mkdir -p $l"_"$t_steps
# #cd $l"_"$t_steps

# #mkdir -p $a.res/l.$l"_"t.$t_steps
# cd $a.res/l.$l"_"t.$t_steps
# pwd

# eval_file="concat_performace.eval"

# # If output doesn't exist, run the evaluation
# if [ ! -f "$eval_file" ]; then
#     result=$(python "$RDIR/AP_AUC_MCC.py" -pf concat.eval)
#     echo "SMM,$a,$result,MC,$l,$t_steps" >> ../../"$eval_file"
# fi

# # Append result to performance_metrics.csv
# if [ -f "$eval_file" ]; then
#     cat "$eval_file" >> "$RESDIR/performance_metrics.csv"
# fi

# cd ..  
# done
# cd ..
# done
  
# done

#!/bin/bash

# Paths
RDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/code"
RESDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/results/SMM"

# Navigate to results directory
cd "$RESDIR/monte_carlo" || exit 1

# Loop over alleles
for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701; do
  for l in 0.1; do
    for t_steps in 20 100; do

      # Navigate into the correct directory
      dir="$a.res/l.$l""_t.$t_steps"
      if [ ! -d "$dir" ]; then
        echo "Skipping missing directory: $dir"
        continue
      fi

      cd "$dir" || continue
      pwd

      eval_file="concat_performance.eval"

      if [ -f concat.eval ]; then
        if [ ! -f "$eval_file" ]; then
          result=$(python "$RDIR/AP_AUC_MCC.py" -pf concat.eval 2>/dev/null)
          if [[ "$result" == *","*","* ]]; then
            echo "SMM,$a,$result,MC,$l,$t_steps" >> ../../"$eval_file"
          else
            echo "Evaluation failed in $dir"
          fi
        else
          echo "Already evaluated: $eval_file"
        fi
      else
        echo "Missing concat.eval in $dir"
      fi

      cd - > /dev/null
    done
  done
done
