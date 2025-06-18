##########################################################################################
####### Constructing Weight Matrix including pseudo counts and sequence weighting #######
##########################################################################################

### Python Imports
import numpy as np
import math
import copy
from pprint import pprint
from argparse import ArgumentParser


### ARG PARSE
parser = ArgumentParser(description="A Pep2mat program")
parser.add_argument("-b", action="store", dest="beta", type=float, default=0, help="Beta value for the algorithm")
parser.add_argument("-w", action="store_true", dest="sequence_weighting", help="Using sequence weighting")
parser.add_argument("-f", action="store", dest="peptides_file", type=str, default = "PSSM/A0101/A0101_bind.dat", help="File with peptides")
args = parser.parse_args()
beta = args.beta
peptides_file = args.peptides_file
sequence_weighting = args.sequence_weighting


### Data directories
# data directory
data_dir = "data/"

# directory for constructed PSSM
output_dir = "data/Outputs/PSSM/Constructed_PSSM"

### Define options for run
# sequence_weighting = True
# sequence_weighting = False
# define weight on pseudo count
# beta = 50
# beta = 0


### Loading alphabet
alphabet_file = data_dir + "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)


### Load Background Frequencies
bg_file = data_dir + "Matrices/bg.freq.fmt"
_bg = np.loadtxt(bg_file, dtype=float)
bg = {} #creates a dictionary for the background frequencies

for i in range(0, len(alphabet)):
    bg[alphabet[i]] = _bg[i]


### Load Blosum62 Matrix
blosum62_file = data_dir + "Matrices/blosum62.freq_rownorm"
_blosum62 = np.loadtxt(blosum62_file, dtype=float).T
blosum62 = {} #creates a dictionary for the blosum62 matrix

for i, letter_1 in enumerate(alphabet):
    
    blosum62[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):
        
        blosum62[letter_1][letter_2] = _blosum62[i, j]


### Load Peptides
# peptides_file = data_dir + "PSSM/A0201.single_lig"
# peptides_file = data_dir + "PSSM/A0201.small_lig"
# peptides_file = data_dir + "PSSM/A0201.large_lig"

peptides_file_path = data_dir + peptides_file
peptides = np.loadtxt(peptides_file_path, dtype=str).tolist()[] # load peptides from file

#obs, added extra [0] slicing as we have lists in li
#saving the length of the peptides.
if len(peptides[0][0]) == 1:
    peptide_length = len(peptides)
    peptides = [peptides]
else:
    peptide_length = len(peptides[0][0])

# check if all peptides are of the same length
for i in range(0, len(peptides)):
    if len(peptides[i][0]) != peptide_length:
        print("Error, peptides differ in length!")


### Initialize Matrix

#Function for initializing matrix.
def initialize_matrix(peptide_length, alphabet):

    init_matrix = [0]*peptide_length # makes a list of zeros with length equal to peptide_length

    for i in range(0, peptide_length):

        row = {}

        for letter in alphabet: 
            row[letter] = 0.0 #adds values of 0 for each possible substitution.

        #fancy way:  row = dict( zip( alphabet, [0.0]*len(alphabet) ) )

        init_matrix[i] = row
        
    return init_matrix


### Amino Acid Count Matrix (c)
c_matrix = initialize_matrix(peptide_length, alphabet)

for position in range(0, peptide_length):
    print(position)

    #peptides[]
        
    for peptide in peptides:
        c_matrix[position][peptide[0][position]] += 1
        # print(peptide[0][position])

### Sequence Weighting
# w = 1 / r * s
# where 
# r = number of different amino acids in column
# s = number of occurrence of amino acid in column
weights = {}

for peptide in peptides:

    # apply sequence weighting
    if sequence_weighting:
    
        w = 0.0
        neff = 0.0
        
        for position in range(0, peptide_length):

            r = 0

            for letter in alphabet:        

                if c_matrix[position][letter] != 0: #if the count of the amino acid is not zero
                    
                    r += 1 #counting how many different amino acids there are in the column

            s = c_matrix[position][peptide[position]] #number of occurrences of the amino acid in the column

            w += 1.0/(r * s) #calculating the weight for the peptide

            neff += r #adding the number of different amino acids in the column to neff
                
        neff = neff / peptide_length # average number of different amino acids per position
  
    # do not apply sequence weighting
    else:
        
        w = 1  
        
        neff = len(peptides)

    weights[peptide] = w

pprint( "W:")
pprint( weights )
pprint( "Nseq:")
pprint( neff )


### Observed Frequencies Matrix (f)

f_matrix = initialize_matrix(peptide_length, alphabet) #initialize zero matrix for frequencies

for position in range(0, peptide_length):
  
    n = 0;
  
    for peptide in peptides:
    
        f_matrix[position][peptide[position]] += weights[peptide] 
    
        n += weights[peptide] #saving the sum of weights (corresponds to n) to normalize the frequencies later
        
    for letter in alphabet: 
        
        f_matrix[position][letter] = f_matrix[position][letter]/n
      
pprint( f_matrix[0] )


### Pseudo Frequencies Matrix (g)

# Remember g(b) = sum f(a)* q(b|a), and blosum[a,b] = q(a|b)

g_matrix = initialize_matrix(peptide_length, alphabet)

for position in range(0, peptide_length):

    for letter_1 in alphabet:
        for letter_2 in alphabet:

          g_matrix[position][letter_1] += f_matrix[position][letter_2] * blosum62[letter_1][letter_2]

pprint(g_matrix[0]) #just printing the first row of the g_matrix


### Combined Frequencies Matrix (p)

p_matrix = initialize_matrix(peptide_length, alphabet)

alpha = neff - 1 # alpha is the number of different amino acids minus one

for position in range(0, peptide_length):

    for a in alphabet:
        p_matrix[position][a] = (alpha * f_matrix[position][a] + beta * g_matrix[position][a]) / (alpha + beta)

pprint(p_matrix[0]) #just printing the first position of the p_matrix


### Log Odds Weight Matrix (w)

w_matrix = initialize_matrix(peptide_length, alphabet)

for position in range(0, peptide_length):
    
    for letter in alphabet:
        if p_matrix[position][letter] > 0: #check if the probability is greater than zero because log(0) is not handled.
            w_matrix[position][letter] = 2 * math.log(p_matrix[position][letter]/bg[letter])/math.log(2)
        else:
            w_matrix[position][letter] = -999.9

pprint(w_matrix[0])

### Write Matrix to PSI-BLAST format
def to_psi_blast(matrix):

    header = ["", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    print ('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*header)) 

    letter_order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    for i, row in enumerate(matrix):

        scores = []

        scores.append(str(i+1) + " A")

        for letter in letter_order:

            score = row[letter]

            scores.append(round(score, 4))

        print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*scores)) 

### convert w_matrix to PSI-BLAST format and print to file

def to_psi_blast_file(matrix, file_name):
    
    with open(file_name, 'w') as file:

        header = ["", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        file.write ('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n'.format(*header)) 

        letter_order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        for i, row in enumerate(matrix):

            scores = []

            scores.append(str(i+1) + " A")

            for letter in letter_order:

                score = row[letter]

                scores.append(round(score, 4))

            file.write('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n'.format(*scores)) 


### convert  w_matrix to PSI-BLAST format
to_psi_blast(w_matrix)

### convert w_matrix to PSI-BLAST format and print to file

# Write out PSSM in Psi-Blast format to file
to_psi_blast_file(w_matrix, file_name=output_dir)



### Evaluation Udkommenteret fordi pep2score klarer det?

# #evaluation_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_2/A0201.eval"
# evaluation_file = data_dir + "PSSM/A0201.eval"
# #evaluation_file = evaluation_upload.values()

# evaluation = np.loadtxt(evaluation_file, dtype=str).reshape(-1,2)
# evaluation_peptides = evaluation[:, 0]
# evaluation_targets = evaluation[:, 1].astype(float)

# evaluation_peptides, evaluation_targets

# %%
# def score_peptide(peptide, matrix):
#     acum = 0
#     for i in range(0, len(peptide)):
#         acum += w_matrix[i][peptide[i]]
#     return acum

# %%
# evaluation_predictions = []
# for evaluation_peptide in evaluation_peptides:
#     evaluation_predictions.append(score_peptide(evaluation_peptide, w_matrix))

# %%
# from scipy.stats import pearsonr
# import matplotlib.pyplot as plt

# pcc = pearsonr(evaluation_targets, evaluation_predictions)
# print("PCC: ", pcc[0])

# plt.scatter(evaluation_targets, evaluation_predictions);

# %%



