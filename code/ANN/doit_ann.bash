#! /bin/bash -f

## Define path to your code directory
RDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/code"

## Define path you where you have placed the HLA data sets
DDIR="/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/data/AllFiles"

# A0201 A0202 A1101 A3001 B0702 B1501 B5401
# Here you can type your allele names
for a in A0101
do

rm -rf $a.res
mkdir -p $a.res

cd $a.res

# Here you can type the encoding to test
for sc in blosum sparse
do

mkdir -p $sc

cd $sc

# Loop over the 4 cross validation configurations
for n in 0 1 2 3
do

# Do training
if [ ! -f $a-$sc-$n.syn ] 
then
	python "$RDIR/ANN/ANN_train.py" -sc $sc -t "$DDIR/$a/f00$n" -e "$DDIR/$a/c00$n" -syn $a-$sc-$n.syn -stop | tee train.log.$n
fi

# Do evaluation
if [ ! -f c00$n.pred ] 
then
	python "$RDIR/ANN/ANN_forward.py" -sc $sc -e "$DDIR/$a/c00$n" -syn $a-$sc-$n.syn  | tee >(grep "PCC:" > c00$n.pcc) | grep -v "PCC:" > c00$n.pred
fi

done

# Do concatinated evaluation
echo $a $sc `cat c00{0..3}.pred | grep -v "#" | gawk '{print $2,$3}' | ../../../xycorr` \
	   `cat c00{0..3}.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print "MSE:", e/n}' ` >> ../summary.txt

cd ..

done

cd ..

done