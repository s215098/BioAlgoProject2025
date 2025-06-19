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

mkdir -p scheme.$sc

cd scheme.$sc

# Loop over the 4 cross validation configurations
for n in 0 1 2 3
do

# Do training
if [ ! -f $a-$sc-$n.syn ] 
then
	python "$RDIR/ANN/ANN_train.py" -sc $sc -t "$DDIR/$a/f00$n" -e "$DDIR/$a/c00$n" -syn $a-$sc-$n.syn -stop | tee train.log.$n
fi

done

# Do concatinated evaluation
# Concat and summarize PCCs for this scheme
# echo "Summary for $a scheme=$sc" >> ../../summary.txt
# cat c00{0..3}.out >> ../../summary.txt
# echo "" >> ../../summary.txt

cd ..

done

cd ..

done