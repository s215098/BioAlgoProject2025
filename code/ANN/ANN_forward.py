# %% [markdown]
# # Predicting with Neural Networks
# 
# ### Fill the parts with X

# %% [markdown]
# ## Python Imports

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from argparse import ArgumentParser
parser = ArgumentParser(description="Predicting with Forward Neural Networks")
parser.add_argument("-e", action="store", dest="evaluation_file", type=str, help="File with evaluation data")
parser.add_argument("-syn", action="store", dest="synfile_name", type=str, help="Name of synaps file")
parser.add_argument("-sc", action="store", dest="scheme", type=str, default="sparse", help="Set encoding scheme (default: sparse)")
args = parser.parse_args()
evaluation_file = args.evaluation_file
synfile_name = args.synfile_name
scheme = args.scheme

# %% [markdown]
# ## Data Imports

# %% [markdown]
# ### DEFINE THE PATH TO YOUR COURSE DIRECTORY

# %%
data_dir = "/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/data/"

# %% [markdown]
# ### define run time parameters

# %%
# Define if we are using blosum or sparse encoding
# blosum_scheme = False
#blosum_scheme = True

# %% [markdown]
# ### Alphabet

# %%
alphabet_file = data_dir + "Matrices/alphabet"
#alphabet_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)

# %% [markdown]
# ### Blosum50 Encoding Scheme

# %%
blosum_file = data_dir + "Matrices/blosum50"
#blosum_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/blosum50"

_blosum50 = np.loadtxt(blosum_file, dtype=float).reshape((24, -1)).T

blosum50 = {}

for i, letter_1 in enumerate(alphabet):
    
    blosum50[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):
        
        blosum50[letter_1][letter_2] = _blosum50[i, j] / 5.0

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
        
        # add a 1 (bias)
        encoded_peptide.append(1)
        
        # store peptide
        encoded_peptides.append(encoded_peptide)
        
    return np.array(encoded_peptides)

# %% [markdown]
# ## Neural Network Functions

# %% [markdown]
# ### Activation (Sigmoid)

# %%
def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))

# %% [markdown]
# ### Forward Propagation

# %%
def forward(X, w1, w2):
    
    # X contains the output from each layer, i.e the input values in the first layer
    # w1 are weights connecting input to hidden, and w2 weights connecting hidden to output
    # In w[i,j]; i is from and j is to
   
    # get dimension, substracting the bias
    input_layer_dim = w1.shape[0] - 1 
    hidden_layer_dim = w2.shape[0] - 1
    
    ################
    # hidden layer #
    ################
    
    # activity of hidden layer
    # Remember z_j = sum_i w(i,j)*input(i)
    for j in range(hidden_layer_dim):
        z = 0.0
        for i in range(input_layer_dim+1):
            # z += XXX
            z += X[0][i] * w1[i, j] # This is the first layer, so we use X[0]. w1 connects input to hidden, so we use w1[i][j] because we have j neurons in the hidden layer.
        # X[1][j] = XXX
        X[1][j] = sigmoid(z) # This is the second layer, so we use X[1]
    
    ################
    # output layer #
    ################
    
    z = 0
    for i in range(hidden_layer_dim+1):
        # z += XXX
        z += X[1][i] * w2[i, 0] # This is the second layer, so we use X[1]. w2 connects hidden to output, so we use w2[i][0] because we have only one output neuron.
    # X[2][0] = XXX
    X[2][0] = sigmoid(z) # This is the third layer, so we use X[2]


# %% [markdown]
# ## Prediction Data

# %%
# evaluation_file = data_dir + "ANN/A2403_evaluation"
#evaluation_file = data_dir + "ANN/A0201_evaluation"
evaluation_data = np.loadtxt(evaluation_file, dtype=str)

peptides = evaluation_data[:, 0]
if scheme == "blosum":
    x_eval = encode(peptides, blosum50, alphabet)
elif scheme == "sparse":
    x_eval = encode(peptides, sparse, alphabet)

y_eval = np.array(evaluation_data[:, 1], dtype=float)

# %% [markdown]
# ## Function to load previously saved Network

# %%
def load_network(file_name):

    f = open(file_name, "r")

    n_line = 0

    weight_list = []

    for line in f:


        # clean and separate line
        sline = line.strip().split()


        # input layer dimension
        if n_line == 1:
            input_layer_dim = int(sline[0])

        # hidden layer dimension    
        if n_line == 2:
            hidden_layer_dim = int(sline[0])

        # output layer dimension
        if n_line == 3:
            output_layer_dim = int(sline[0])

        # model weights
        if n_line >= 5:
            for i in range(0, len(sline)):
                weight_list.append(float(sline[i]))

        n_line += 1

    # HIDDEN LAYER WEIGHTS
    # w_h[i, j] is the weight that links input's feature "i" to neuron "j" of the hidden layer        
    w_h_load = np.zeros(shape=(input_layer_dim+1, hidden_layer_dim))

    for i in range(0, (input_layer_dim+1)*hidden_layer_dim, hidden_layer_dim):

        for j in range(0, hidden_layer_dim):

            row = i // hidden_layer_dim

            w_h_load[row, j] = weight_list[i+j]

            
    # OUTPUT LAYER WEIGHTS
    # w_o[i, j] is the weight that links hidden layer's neuron "i" to neuron "j" of the output layer
    w_o_load = np.zeros(shape=(hidden_layer_dim+1, output_layer_dim))

    w_h_end = (input_layer_dim+1) * hidden_layer_dim

    for i in range(w_h_end, w_h_end+hidden_layer_dim+1, output_layer_dim):

        for j in range(0, output_layer_dim):

            row = (i - w_h_end) // output_layer_dim
            w_o_load[row, j] = weight_list[i+j]
            
            
    # return weight matrices
    return w_h_load, w_o_load

# %% [markdown]
# ## Main code

# %%
# Load network
# synfile_name = data_dir + "ANN/A2403_sp.syn"
# synfile_name = data_dir + "ANN/A2403_bl.syn"
# synfile_name = data_dir + "ANN/A0201_sp.syn"
# synfile_name = data_dir + "ANN/A0201_bl.syn"
w_h, w_o = load_network(synfile_name)

# X matrix 
input_layer_dim = w_h.shape[0]
hidden_layer_dim = w_o.shape[0]
output_layer_dim = w_o.shape[1]

# Find max network dimensions
X_dim = max(input_layer_dim, hidden_layer_dim, output_layer_dim)
X = np.zeros(shape=(3, X_dim))

# The last column in each X layer is set to 1 to deal with the bias weights
X[0][input_layer_dim-1] = 1.0 
X[1][hidden_layer_dim-1] = 1.0
    
# data for plotting
y_preds_eval = []

# loop
for i in range(0, len(x_eval)):        

    # fetch training point
    x = x_eval[i]
    y = y_eval[i]

    if len(x) == input_layer_dim:
        
        X[0] = x

        # forward propagation
        forward(X, w_h, w_o)
        # y_pred = XXX
        y_pred = X[2][0] # This is the third layer, so we use X[2][0] to get the output prediction
        
        y_preds_eval.append(y_pred)
        
        print(peptides[i], y_pred, y)
    else:
        print("Error. Peptide length", len(x),"does not match network sizs", input_layer_dim, "Skip")

# store training performance
eval_perf = pearsonr(y_eval, np.asarray(y_preds_eval))[0]

# PERFORMANCE REPORT
# fig = plt.figure(figsize=(5, 5), dpi = 70)

# plt.scatter(y_preds_eval, y_eval)
# plt.ylabel("Target Value", fontsize=10);
# plt.xlabel("Prediction Value", fontsize=10);

# print performance
print("# Prediction PCC:", round(eval_perf, 4))

# %% [markdown]
# You have an input with 9 amino acids. How does that relate to the number of neurons in the first layer?
# - 9 * 20 = 180
# - There are 20 possible AAs per sequence and the sequence is 9 AAs long
# 
# How many hidden neurons does the network have?
# - The network has 5 hidden neurons
# 
# And how many output values?
# - There is one output value
# 
# 
# Can you understand the number of synaps weights in the files, i.e.,
# 
# ```cat ../data/A2403_sp.syn | grep -v TEST | grep -v ":" | wc```
# 
# The second column in this command gives the number of weights in the synaps file. Can you make sense of this number (911)?
# - There are 900 weights from the input layer to the hidden layer. 5 from the hidden to the output layer. And then one bias per layer, meaning 5 weights to the hidden layer and one to the output.
# 
# What is the predictive performance of the neural network (in terms of the Pearsons correlation)? And how does this value compare to the value listed in the first line of the synaps file (T_PCC column)?
# 
# head ANN/A2403_sp.syn 
# TESTRUNID EPOCH: 99 L_PCC: 0.818935942557 L_ERR: 0.0130974863122 T_PCC: 0.644656416379 T_ERR: 0.0164147955673
# 
# - The PCC is 0.6447, which is the same as the T_PCC

# %%



