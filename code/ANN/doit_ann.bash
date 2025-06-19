#! /bin/bash -f

## Define path to your code directory
RDIR="ANN_forward_clean.py"

## Define path you where you have placed the HLA data sets
DDIR="../../data/AllFiles"

# Here you can type your allele names
for a in A0201 B0801 B5701 
do

mkdir -p $a.res

cd $a.res

# Here you can type the lambdas to test
for l in 0 0.02 0.05 0.1
do

mkdir -p l.$l

cd l.$l

# Loop over the 5 cross validation configurations
for n in 0 1 2 3 4 
do

# Do training
if [ ! -f mat.$n ] 
then
	python $RDIR/SMM/smm_gradient_descent.py -l $l -t $DDIR/$a/f00$n -e $DDIR/$a/c00$n | grep -v "#" > mat.$n
fi

# Do evaluation
if [ ! -f c00$n.pred ] 
then
	python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
fi

done

# Do concatinated evaluation
echo $a $l `cat c00{0..4}.pred | grep -v "#" | gawk '{print $2,$3}' | ../../xycorr` \
	   `cat c00{0..4}.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' `

cd ..

done

cd ..

done