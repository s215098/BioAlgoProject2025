# %% [markdown]
# # SMM with Gradient Descent

# %% [markdown]
# ## Python Imports

# %%
import numpy as np
import random
import copy
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from argparse import ArgumentParser

# %% [markdown]
# ## Arguments from command line

# %%
parser = ArgumentParser(description="SMM gradient descent")
parser.add_argument("-l", action="store", dest="LAMB", default=0.01, type=float, help="Lambda")
parser.add_argument("-t", action="store", dest="TRAINING_FILE", type=str, help="File with training data")
parser.add_argument("-epi", action="store", dest="EPSILON", default=0.05, type=float, help="Epsilon")
parser.add_argument("-s", action="store", dest="SEED", default=1, type=int, help="Seed for random numbers")
parser.add_argument("-i", action="store", dest="EPOCHS", default=100, type=int, help="Number of epochs to train")
parser.add_argument("-w", action="store", dest="w_bound", default=0.1, type=int, help="Weight bound")

args = parser.parse_args()
lamb = args.LAMB
training_file = args.TRAINING_FILE
epsilon = args.EPSILON
seed = args.SEED
epochs = args.EPOCHS
w_bound = args.w_bound
# %% [markdown]
# ## Data Imports

# %% [markdown]
# ## DEFINE THE PATH TO YOUR COURSE DIRECTORY

# RUN: python code/SMM/smm_gradient_descent.py -t data/SMM/A0201_training
# %%
data_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/data/"

# %% [markdown]
# ### Training Data

# %%
#training_file = data_dir + "SMM/A0201_training"
#training_file = data_dir + "SMM/A2403_training"

training = np.loadtxt(training_file, dtype=str)

# %% [markdown]
# ### Evaluation Data


# %% [markdown]
# ### Alphabet

# %%
alphabet_file = data_dir + "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)
#print(f"Alphabet: \n {alphabet}")

# %% [markdown]
# ### Sparse Encoding Scheme

# %%
sparse_file = data_dir + "Matrices/sparse"
_sparse = np.loadtxt(sparse_file, dtype=float)
sparse = {}

for i, letter_1 in enumerate(alphabet):
    
    sparse[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):
        
        sparse[letter_1][letter_2] = _sparse[i, j]
#print(f"Sparse encoding scheme: \n {sparse}")

# %% [markdown]
# ## Peptide Encoding

# %%
# Puts 1 if peptide_letter is the same as alphabet_letter, otherwise 0
# The peptides are encoded as numarical vectors
def encode(peptides, encoding_scheme, alphabet):
    
    encoded_peptides = []

    for peptide in peptides:

        encoded_peptide = []

        for peptide_letter in peptide:

            for alphabet_letter in alphabet:

                encoded_peptide.append(encoding_scheme[peptide_letter][alphabet_letter])

        encoded_peptides.append(encoded_peptide)
        
    return np.array(encoded_peptides)

# %% [markdown]
# ## Error Function

# %%
def cumulative_error(peptides, y, lamb, weights):

    error = 0
    
    for i in range(0, len(peptides)):
        
        # get peptide
        peptide = peptides[i]

        # get target prediction value
        y_target = y[i]
        
        # get prediction
        y_pred = np.dot(peptide, weights)
            
        # calculate error
        error += 1.0/2 * (y_pred - y_target)**2
        
    gerror = error + lamb*np.dot(weights, weights)
    error /= len(peptides)
        
    return gerror, error

# %% [markdown]
# ## Predict value for a peptide list

# %%
def predict(peptides, weights):

    pred = []
    
    for i in range(0, len(peptides)):
        
        # get peptide
        peptide = peptides[i]
        
        # get prediction
        y_pred = np.dot(peptide, weights)
        
        pred.append(y_pred)
        
    return pred

# %% [markdown]
# ## Calculate MSE between two vectors

# %%
def cal_mse(vec1, vec2):
    
    mse = 0
    
    for i in range(0, len(vec1)):
        mse += (vec1[i] - vec2[i])**2
        
    mse /= len(vec1)
    
    return( mse)

# %% [markdown]
# ## Gradient Descent

# %%
def gradient_descent(y_pred, y_target, peptide, weights, lamb_N, epsilon):
    
    # do is dE/dO
    do = (y_pred - y_target)
        
    for i in range(0, len(weights)):
        
        de_dw_i = do * peptide[i] + 2*lamb_N * weights[i]

        weights[i] -= epsilon * de_dw_i

# %% [markdown]
# ## Main Loop
# 
# 

# %%
# Random seed 
np.random.seed(seed)

# peptides
peptides = training[:, 0]
# Create numerical vectors for peptides
peptides = encode(peptides, sparse, alphabet)
N = len(peptides)

# target values
y = np.array(training[:, 1], dtype=float)

'''
#evaluation peptides
evaluation_peptides = evaluation[:, 0]
evaluation_peptides = encode(evaluation_peptides, sparse, alphabet)

#evaluation targets
evaluation_targets = np.array(evaluation[:, 1], dtype=float)
'''

# weights
input_dim  = len(peptides[0])
output_dim = 1
w_bound = 0.1
weights = np.random.uniform(-w_bound, w_bound, size=input_dim)

# training epochs
#epochs = 100

# regularization lambda
#lamb = 0.01
#lamb = 1
#lamb = 10

# regularization lambda per target value
lamb_N = lamb/N

# learning rate
#epsilon = 0.01
#epsilon = 0.05

# error  plot
gerror_plot = []
mse_plot = []
train_mse_plot = []
eval_mse_plot = []
train_pcc_plot = []
eval_pcc_plot = []

# for each training epoch
for e in range(0, epochs):

    # for each peptide
    for i in range(0, N):

        # random index
        ix = np.random.randint(0, N)
        
        # get peptide       
        peptide = peptides[ix]

        # get target prediction value
        y_target = y[ix]
       
        # get initial prediction
        y_pred = np.dot(peptide, weights)

        # gradient descent 
        gradient_descent(y_pred, y_target, peptide, weights, lamb_N, epsilon)
'''
    # compute error
    gerr, mse = cumulative_error(peptides, y, lamb, weights) 
    gerror_plot.append(gerr)
    mse_plot.append(mse)
    
    # predict on training data
    train_pred = predict( peptides, weights )
    train_mse = cal_mse( y, train_pred )
    train_mse_plot.append(train_mse)
    train_pcc = pearsonr( y, train_pred )
    train_pcc_plot.append( train_pcc[0] )
        
    # predict on evaluation data
    eval_pred = predict(evaluation_peptides, weights )
    eval_mse = cal_mse(evaluation_targets, eval_pred )
    eval_mse_plot.append(eval_mse)
    eval_pcc = pearsonr(evaluation_targets, eval_pred)
    eval_pcc_plot.append( eval_pcc[0] )
    
    #print ("Epoch: ", e, "Gerr:", gerr, train_pcc[0], train_mse, eval_pcc[0], eval_mse)


# %% [markdown]
# ## Error Plot

# %%
fig = plt.figure(figsize=(10, 10), dpi= 80)

x = np.arange(0, len(gerror_plot))

plt.subplot(2, 2, 1)
plt.plot(x, gerror_plot)
plt.ylabel("Global Error", fontsize=10);
plt.xlabel("Iterations", fontsize=10);

plt.subplot(2, 2, 2)
plt.plot(x, mse_plot)
plt.ylabel("MSE", fontsize=10);
plt.xlabel("Iterations", fontsize=10);


x = np.arange(0, len(train_mse_plot))

plt.subplot(2, 2, 3)
plt.plot(x, train_mse_plot, label="Training Set")
plt.plot(x, eval_mse_plot, label="Evaluation Set")
plt.ylabel("Mean Squared Error", fontsize=10);
plt.xlabel("Iterations", fontsize=10);
plt.legend(loc='upper right');


plt.subplot(2, 2, 4)
plt.plot(x, train_pcc_plot, label="Training Set")
plt.plot(x, eval_pcc_plot, label="Evaluation Set")
plt.ylabel("Pearson Correlation", fontsize=10);
plt.xlabel("Iterations", fontsize=10);
plt.legend(loc='upper left');
#plt.show()
'''
# %% [markdown]
# ## Get PSSM Matrix

# %% [markdown]
# ### Vector to Matrix

# %%
# our matrices are vectors of dictionaries
def vector_to_matrix(vector, alphabet):
    
    rows = int(len(vector)/len(alphabet))
    
    matrix = [0] * rows
    
    offset = 0
    
    for i in range(0, rows):
        
        matrix[i] = {}
        
        for j in range(0, 20):
            
            matrix[i][alphabet[j]] = vector[j+offset] 
        
        offset += len(alphabet)

    return matrix

# %% [markdown]
# ### Matrix to Psi-Blast

# %%
def to_psi_blast(matrix):

    # print to user
    
    header = ["", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*header)) 

    letter_order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    for i, row in enumerate(matrix):

        scores = []

        scores.append(str(i+1) + " A")

        for letter in letter_order:

            score = row[letter]

            scores.append(round(score, 4))

        print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*scores)) 


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

# %% [markdown]
# ### Print

# %%
matrix = vector_to_matrix(weights, alphabet)
to_psi_blast(matrix)
#print(matrix)
'''
file_name = "code/SMM/data_SMM/SMM_PSSM_matrix"
to_psi_blast_file(matrix, file_name)
'''


