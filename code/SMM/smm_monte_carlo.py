# %% [markdown]
# # SMM with Monte Carlo

# %% [markdown]
# ## Python Imports

# RUN: python code/SMM/smm_monte_carlo.py -t data/SMM/A0201_training -e data/SMM/A0201_evaluation
# %%
import numpy as np
import random
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from argparse import ArgumentParser

# %% [markdown]
# ## Arguments from command line
# 

# %%
parser = ArgumentParser(description="SMM Monte Carlo")
parser.add_argument("-l", action="store", dest="LAMB", default=1, type=float, help="Lambda")
parser.add_argument("-t", action="store", dest="TRAINING_FILE", type=str, help="File with training data")
parser.add_argument("-i", action="store", dest="ITERS", default=100, type=int, help="Number of epochs to train")
parser.add_argument("-s", action="store", dest="SEED", default=1, type=int, help="Seed for random numbers")
parser.add_argument("-Ts", action="store", dest="T_I", default=0.01, type=float, help="Start temperature")
parser.add_argument("-Te", action="store", dest="T_F", default=0.000001, type=float, help="End temperature")
parser.add_argument("-nT", action="store", dest="T_STEPS", default=10, type=int, help="Number of temperature steps")

args = parser.parse_args()
lamb = args.LAMB
training_file = args.TRAINING_FILE
iters = args.ITERS
seed = args.SEED
T_i = args.T_I
T_f = args.T_F
T_steps = args.T_STEPS


# %% [markdown]
# ## Data Imports

# %% [markdown]
# ## DEFINE THE PATH TO YOUR COURSE DIRECTORY

# %%
data_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/data/"

# %% [markdown]
# ## Define run time parameters

# %%
# temperature vector
#T_i = 0.01
#T_f = 0.000001
#T_steps = 10
T_delta = (T_f - T_i) / T_steps
T = np.linspace(T_i,T_f,T_steps )

# iterations
#iters = 1000

# regularization lambda
#lamb = 1
#lamb = 0.00001
#lamb = 10

# %% [markdown]
# ### Training Data

# %%
#training_file = data_dir + "SMM/A0201_training"
#training_file = data_dir + "SMM/A2403_training"
training = np.loadtxt(training_file, dtype=str)

# %% [markdown]
# ### Evaluation Data

# %%
#evaluation_file = data_dir + "SMM/A0201_evaluation"
#evaluation_file = data_dir + "SMM/A2403_evaluation"
#evaluation = np.loadtxt(evaluation_file, dtype=str)

# %% [markdown]
# ### Alphabet

# %%
alphabet_file = data_dir + "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)

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

# %% [markdown]
# ## Peptide Encoding

# %%
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
# ## Main Loop

# %%
# Random seed
np.random.seed(seed)

# peptides
peptides = training[:, 0]
peptides = encode(peptides, sparse, alphabet)

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

# error plots
gerror_plot = []
mse_plot = []
train_mse_plot = []
eval_mse_plot = []
train_pcc_plot = []
eval_pcc_plot = []

perturbation_value = 0.1
# The scale variable defines the scale of the changes in de. This is defined empirically for each problem investigated
scale = 1.0/100

number_of_tries = 0
number_of_accepted = 0
        
# calculate initial error
gerror_initial, mse = cumulative_error(peptides, y, lamb, weights)
        
# for each temperature
for t in T:
  
    # for each iteration
    for i in range(0, iters):
        
        
        # get 2 random weight indexes
        weight_index_1 = np.random.randint(len(weights))
        weight_index_2 = np.random.randint(len(weights))
    
        # ensure they are different
        # we don't want to perturb the same weight twice
        while weight_index_1 == weight_index_2:
            weight_index_2 = np.random.randint(len(weights))
        
        # store original weight values
        original_weight_1 = weights[weight_index_1]    #XX
        original_weight_2 = weights[weight_index_2]    #XX
       
    
        # apply random perturbation to both weights
        # it is not defined how to do it, but here it's chosen that one weight is increased and the other is decreased
        # you can choose how you want and what makes sense
        perturbation = np.random.uniform(0, perturbation_value)
        weights[weight_index_1] += perturbation    #XX
        weights[weight_index_2] -= perturbation    #XX
        
        
        # calculate new error
        gerror_new, mse = cumulative_error(peptides, y, lamb, weights) 
            
        # compute error difference
        de = (gerror_new - gerror_initial) * scale    #XX 

        # compute acceptance probability
        if ( de < 0): 
            p = 1
        else:
            p = np.exp(-de/t)    #XX
            
        # throw coin
        coin = np.random.uniform(0.0, 1.0, 1)[0] 
        
        # weight change is accepted
        if ( coin < p ):
            gerror_initial = gerror_new
            gerror_plot.append(gerror_new)
            mse_plot.append(mse)
            number_of_accepted += 1
        # weight change is declined, restore previous weights
        else: 
            weights[weight_index_1] = original_weight_1 
            weights[weight_index_2] = original_weight_2

        number_of_tries += 1
        
        # define size of move so that on avarage 50% are accepted
        if number_of_tries == 100:
            
            if float(number_of_accepted)/number_of_tries > 0.5:
                perturbation_value *= 1.1
            else: 
                perturbation_value *= 0.9
            
            number_of_tries = 0
            number_of_accepted = 0
            
        # predict on training data
        train_pred = predict(peptides, weights)
        train_mse = cal_mse(y, train_pred)
        train_pcc = pearsonr(y, train_pred)[0]
        train_mse_plot.append(train_mse)
        train_pcc_plot.append(train_pcc)
        '''
        # predict on evaluation data
        eval_pred = predict(evaluation_peptides, weights)
        eval_mse = cal_mse(evaluation_targets, eval_pred)
        eval_pcc = pearsonr(evaluation_targets, eval_pred)[0]
        eval_mse_plot.append(eval_mse)
        eval_pcc_plot.append(eval_pcc)    
        '''
    #print ("t:", t, gerror_new, perturbation_value, train_mse, train_pcc, eval_mse, eval_pcc)

# %% [markdown]
# ## Error Plot
'''
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
# 
# 

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

# %% [markdown]
# ### Print

# %%
matrix = vector_to_matrix(weights, alphabet)
to_psi_blast(matrix)

# %% [markdown]
# ## Performance Evaluation
'''
# %%
evaluation_peptides = evaluation[:, 0]
evaluation_peptides = np.array(encode(evaluation_peptides, sparse, alphabet))

evaluation_targets = np.array(evaluation[:, 1], dtype=float)

y_pred = []
for i in range(0, len(evaluation_peptides)): 
    y_pred.append(np.dot(evaluation_peptides[i], weights))
    
y_pred = np.array(y_pred)

# %%
pcc = pearsonr(evaluation_targets, np.array(y_pred))
print("PCC: ", pcc[0])
plt.scatter(y_pred, evaluation_targets);

#### PCC ####
# lambda = 0.01 pcc = 0.5106454185156161
# lambda = 1 pcc = 0.5176533901125014
# lambda = 10 pcc = 0.49981612680472554


# %% [markdown]
# 


'''