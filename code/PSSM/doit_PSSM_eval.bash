#! /bin/bash -f

## Define path to your code directory
RDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/code"

## Define path you where you have placed the HLA data sets
DDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/data/PSSM"

RESDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/PSSM"

cd ../../
cd "$RESDIR"

find "$RESDIR" -name "e00*.eval" -delete
find "$RESDIR" -name "e00*.pcc" -delete
find "$RESDIR" -name "summary_eval.txt" -delete

# Here you can type your allele names
# A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
do

# Create and move into a results directory for each allele.
# rm -rf $a.res # remove the previous created result files
# mkdir -p $a.res #-p is just mkdir options that ensures it is done "safely"

cd $a.res

# Here you can type the to test
for beta in 5.0 10 50 #defining hyperparameters (beta for in PSSM)
do

# Create and move into a subdirectory for each hyperparameter value.
#mkdir -p b.$beta
cd b.$beta

# Loop over the 5 cross validation configurations
for n in 1 2 3 4
do

# # Do training
# if [ ! -f mat.$n ] #if mat.$n file does not exist.
# then
# 	# Train the model with fold n using training (f00n) and evaluation (c00n) sets.
#     # Save learned matrix in mat.$n (ignore comments with grep).
# 	python $RDIR/PSSM/pep2mat.py -b $beta -f $DDIR/$a/f00$n.csv | grep -v "#" > mat.$n
# fi

# Do evaluation
if [ ! -f e00$n.eval ] #if c00$n.pred file does not exist.
then
	# Use the learned matrix to score peptides from the test set.
    # Save predictions to c00n.pred (ignore lines starting with "PCC:").
	# python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
	python "$RDIR/PSSM/pep2score.py" -mat mat.$n -f "$DDIR/$a/e000.csv" | tee >(grep "PCC:" > e00$n.pcc) | grep -v "PCC:" > e00$n.eval

	# Add binary columns: 1 if value > 0.426, else 0
	awk '{bin2=($2>0.426)?1:0; bin3=($3>0.426)?1:0; print $0, bin2, bin3}' e00$n.eval > e00$n.eval.tmp && mv e00$n.eval.tmp e00$n.eval
fi

done

cat e00{1..4}.eval | grep -v "#" > concat.eval

# Do concatenated evaluation over all 5 folds
# Print allele, lambda, followed by:
# - calculated correlation between predicted and actual values (by xycorr)
# - MSE (Mean Squared Error) over all predictions.
echo $a $beta `cat e00{1..4}.eval | grep -v "#" | gawk '{print $2,$3}' | "$RDIR/xycorr"` \
	   `cat e00{1..4}.eval | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print "MSE:", e/n}' ` >> ../summary_eval.txt

cd ..

done

cd ..

done