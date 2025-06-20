#! /bin/bash -f

## Define path to your code directory
#RDIR="/Users/kristinetoftjohansen/Desktop/Algo/BioAlgoProject2025/code"
RDIR="/Users/mathildedue/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/master_bioinformatics/1.semester/22125_algorithms_in_bioinformatics/BioAlgoProject2025/code"



## Define path you where you have placed the HLA data sets
#DDIR="/Users/kristinetoftjohansen/Desktop/Algo/BioAlgoProject2025/data/PSSM"
DDIR="/Users/mathildedue/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/master_bioinformatics/1.semester/22125_algorithms_in_bioinformatics/BioAlgoProject2025/data/PSSM"


# Here you can type your allele names
# A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
#for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
for a in A0202
do

# Create and move into a results directory for each allele.
rm -rf $a.res_time # remove the previous created result files
mkdir -p $a.res_time #-p is just mkdir options that ensures it is done "safely"
cd $a.res_time

# Here you can type the to test
for beta in 10 #defining hyperparameters (beta for in PSSM)
do

# Create and move into a subdirectory for each hyperparameter value.
mkdir -p b.$beta
cd b.$beta

# Loop over the 5 cross validation configurations
for n in 1 2 3 4
do

# Do training
if [ ! -f mat.$n ] #if mat.$n file does not exist.
then
	# Train the model with fold n using training (f00n) and evaluation (c00n) sets.
    # Save learned matrix in mat.$n (ignore comments with grep).
	python $RDIR/PSSM/pep2mat.py -b $beta -f $DDIR/$a/f00$n.csv | tee >(grep -i "#Time" >> all_times.log) | grep -v "#" > mat.$n
fi

# Do evaluation
if [ ! -f c00$n.pred ] #if c00$n.pred file does not exist.
then
	# Use the learned matrix to score peptides from the test set.
    # Save predictions to c00n.pred (ignore lines starting with "PCC:").
	# python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
	python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n.csv | tee >(grep "PCC:" > c00$n.pcc) | grep -v "PCC:" > c00$n.pred
fi

done

# Do concatenated evaluation over all 5 folds
# Print allele, lambda, followed by:
# - calculated correlation between predicted and actual values (by xycorr)
# - MSE (Mean Squared Error) over all predictions.
echo $a $beta `cat c00{1..4}.pred | grep -v "#" | gawk '{print $2,$3}' | $RDIR/xycorr` \
	   `cat c00{1..4}.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print "MSE:", e/n}' ` >> ../summary.txt

cd ..

done

cd ..

done