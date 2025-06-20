#! /bin/bash -f

## Define path to your code directory
RDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/code"

## Define path you where you have placed the HLA data sets
DDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/data/all_files_new"

rm -rf monte_carlo
mkdir -p monte_carlo
cd monte_carlo

# A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
# Here you can type your allele names
for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
do

mkdir -p $a.res

#cd $a.res

# Here you can type the lambdas to test
for l in 0 0.02 0.1 0.001
do

for t_steps in 10 20 50 100
do

#mkdir -p $l"_"$t_steps
#cd $l"_"$t_steps

mkdir -p $a.res/l.$l"_"t.$t_steps
cd $a.res/l.$l"_"t.$t_steps

# Loop over the 4 cross validation configurations
for n in 1 2 3 4
do

# Do training
if [ ! -f mat.$n ] 
then
	python "$RDIR/SMM/smm_monte_carlo.py" -l $l -nT $t_steps -t "$DDIR/$a/f00$n".csv | grep -v "#" > mat.$n
fi

# Do test
if [ ! -f c00$n.pred ] 
then
	python "$RDIR/PSSM/pep2score.py" -mat mat.$n -f  "$DDIR/$a/c00$n".csv | tee >(grep "PCC:" > c00$n.pcc) | grep -v "PCC:" > c00$n.pred
fi

done

# Do concatinated test
echo $a $l $t_steps `cat c00{1..4}.pred | grep -v "#" | gawk '{print $2,$3}' | ../../../../xycorr` \
	   `cat c00{1..4}.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print "MSE:", e/n}' ` >> ../summary.txt

cd ..

cd ..
done

done
done