#! /bin/bash -f

## Define path to your code directory
RDIR="../"

## Define path you where you have placed the HLA data sets
DDIR="../../data/PSSM/"

# Here you can type your allele names
for a in A0201 A3002 
do

# Create and move into a results directory for each allele.
mkdir -p $a.res
cd $a.res

# Here you can type the to test
for l in 0.1 0.02 #defining hyperparameters. 
do

# Create and move into a subdirectory for each hyperparameter value.
mkdir -p l.$l
cd l.$l

# Loop over the 5 cross validation configurations
for n in 0 1 2 3 4 
do

# Do training
if [ ! -f mat.$n ] #if mat.$n file does not exist.
then
	# Train the model with fold n using training (f00n) and evaluation (c00n) sets.
    # Save learned matrix in mat.$n (ignore comments with grep).
	python $RDIR/smm_gradient_descent.py -l $l -t $DDIR/$a/f00$n -e $DDIR/$a/c00$n | grep -v "#" > mat.$n
	
fi

# Do evaluation
if [ ! -f c00$n.pred ] #if c00$n.pred file does not exist.
then
	# Use the learned matrix to score peptides from the test set.
    # Save predictions to c00n.pred (ignore lines starting with "PCC:").
	python $RDIR/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
fi

done

# Do concatenated evaluation over all 5 folds
# Print allele, lambda, followed by:
# - calculated correlation between predicted and actual values (by xycorr)
# - MSE (Mean Squared Error) over all predictions.
echo $a $l `cat c00{0..4}.pred | grep -v "#" | gawk '{print $2,$3}' | xycorr` \
	   `cat c00{0..4}.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' `

cd ..

done

cd ..

done
