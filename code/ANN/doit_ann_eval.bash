#! /bin/bash -f

## Define path to your code directory
#RDIR="/Users/mathildedue/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/master_bioinformatics/1.semester/22125_algorithms_in_bioinformatics/BioAlgoProject2025/code"

## Define path you where you have placed the HLA data sets
#DDIR="/Users/mathildedue/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/master_bioinformatics/1.semester/22125_algorithms_in_bioinformatics/BioAlgoProject2025/data/all_files_new"

#RESDIR="/Users/mathildedue/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/master_bioinformatics/1.semester/22125_algorithms_in_bioinformatics/BioAlgoProject2025/results/ANN"

RDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/code"

## Define path you where you have placed the HLA data sets
DDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/data/all_files_new"

RESDIR="/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/ANN"

cd ../../
cd "$RESDIR"

find "$RESDIR" -name "e00*.eval" -delete
find "$RESDIR" -name "e00*.pcc" -delete
find "$RESDIR" -name "summary_eval.txt" -delete

# A0201 A0202 A1101 A3001 B0702 B1501 B5401
# Here you can type your allele names
for a in A0201 A0202 A1101 A3001 B0702 B1501 B5401 B5701
do

cd $a.res

# Here you can type the encoding to test
for sc in blosum
do

cd $sc

# for i in 50 100 200
# do
# mkdir -p epocs.$i
# cd epocs.$i

for nh in 5
do 

cd hidden.$nh

for epi in 0.05 0.1
do

cd epi.$epi

# Loop over the 4 cross validation configurations
for n in 1 2 3 4
do

# # Do training
# if [ ! -f $a-$sc-$nh-$epi-$n.syn ] 
# then
# 	python "$RDIR/ANN/ANN_train.py" -sc $sc -nh $nh -epi $epi -t "$DDIR/$a/f00$n.csv" -e "$DDIR/$a/c00$n.csv" -syn $a-$sc-$nh-$epi-$n.syn -stop | tee train.log.$n
# fi

# Do evaluation
if [ ! -f e00$n.eval ] 
then
	python "$RDIR/ANN/ANN_forward.py" -sc $sc -e "$DDIR/$a/e000.csv" -syn $a-$sc-$nh-$epi-$n.syn  | tee >(grep "PCC:" > e00$n.pcc) | grep -v "PCC:" > e00$n.eval
	
	# Add binary columns: 1 if value > 0.426, else 0
	awk '{bin2=($2>0.426)?1:0; bin3=($3>0.426)?1:0; print $0, bin2, bin3}' e00$n.eval > e00$n.eval.tmp && mv e00$n.eval.tmp e00$n.eval

fi

done

# Save concatenated predictions
cat e00{1..4}.eval | grep -v "#" > concat.eval

# Do concatinated evaluation
echo $a $sc Hidden neurons: $nh Epsilon: $epi `cat e00{1..4}.eval | grep -v "#" | gawk '{print $2,$3}' | "$RDIR/xycorr"` \
	   `cat e00{1..4}.eval | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print "MSE:", e/n}' ` >> ../../../summary_eval.txt


cd ..

done

cd ..

done

cd ..

done

cd ..

done

# cd ..

# done
