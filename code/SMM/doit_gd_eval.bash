#! /bin/bash -f

## Define path to your code directory
RDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/code"

## Define path you where you have placed the HLA data sets
DDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/data/all_files_new"

RESDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/SMM"

#rm -rf gradient_decent
#mkdir -p gradient_decent

cd ../../
cd "$RESDIR/gradient_decent"

find "$RESDIR/gradient_decent" -name "e00*.eval" -delete
find "$RESDIR/gradient_decent" -name "e00*.pcc" -delete
find "$RESDIR/gradient_decent" -name "summary_eval.txt" -delete

# A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
# Here you can type your allele names
for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
do

#mkdir -p $a.res

cd $a.res

# Here you can type the lambdas to test
# 0 0.02 0.1 0.001
for l in 0 0.02 0.1 0.001    # Just choose the best one (since it's fast we do all of them)
do

#mkdir -p l.$l

cd l.$l

# Loop over the 4 cross validation configurations
for n in 1 2 3 4
do

# Do evaluation
if [ ! -f e00$n.eval ]
then
	python "$RDIR/PSSM/pep2score.py" -mat mat.$n -f "$DDIR/$a/e000.csv" | tee >(grep "PCC:" > e00$n.pcc) | grep -v "PCC:" > e00$n.eval
	
	# Add binary columns: 1 if value > 0.426, else 0
	awk '{bin2=($2>0.426)?1:0; bin3=($3>0.426)?1:0; print $0, bin2, bin3}' e00$n.eval > e00$n.eval.tmp && mv e00$n.eval.tmp e00$n.eval
fi

done

cat e00{1..4}.eval | grep -v "#" > concat.eval

# Do concatinated evaluation
echo $a $l `cat e00{1..4}.eval | grep -v "#" | gawk '{print $2,$3}' | "$RDIR/xycorr"` \
	   `cat e00{1..4}.eval | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print "MSE:", e/n}' ` >> ../summary_eval.txt

cd ..

done

cd ..

done
