#! /bin/bash -f

## Define path to your code directory
#RDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/code"
RDIR="/Users/mathildedue/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/master_bioinformatics/1.semester/22125_algorithms_in_bioinformatics/BioAlgoProject2025/code"


## Define path you where you have placed the HLA data sets
#DDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/data/all_files_new"

DDIR="/Users/mathildedue/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/master_bioinformatics/1.semester/22125_algorithms_in_bioinformatics/BioAlgoProject2025/data/all_files_new"

rm -rf gradient_decent
mkdir -p gradient_decent
cd gradient_decent

# A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
# Here you can type your allele names
#for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
for a in A0202 
do

mkdir -p $a.res_time

cd $a.res_time

# Here you can type the lambdas to test
for l in 0.1
do

mkdir -p l.$l

cd l.$l

# Loop over the 4 cross validation configurations
for n in 1 2 3 4
do

# Do training
if [ ! -f mat.$n ] 
then
	python "$RDIR/SMM/smm_gradient_descent.py" -l $l -t "$DDIR/$a/f00$n".csv | tee >(grep -i "#Time" >> all_times.log) | grep -v "#" > mat.$n
fi

# Do test
if [ ! -f c00$n.pred ]
then
	#python "$RDIR/PSSM/pep2score.py" -mat mat.$n -f "$DDIR/$a/c00$n" | grep -v "PCC:" > c00$n.pred
	python "$RDIR/PSSM/pep2score.py" -mat mat.$n -f "$DDIR/$a/c00$n".csv | tee >(grep "PCC:" > c00$n.pcc) | grep -v "PCC:" > c00$n.pred
fi

done

# Do concatinated test
echo $a $l `cat c00{1..4}.pred | grep -v "#" | gawk '{print $2,$3}' | $RDIR/xycorr` \
	   `cat c00{1..4}.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print "MSE:", e/n}' ` >> ../summary.txt

cd ..

done

cd ..

done
