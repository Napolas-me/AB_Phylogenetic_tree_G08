import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

file = "./Species/alighn_40_k.aln-clustalw"
# Open the alignment file as a MultipleSeqAlignment object 

with open(file,"r") as clw: 
    alignment = AlignIO.read(clw,"clustal")
    
# Open and initiate the Distance Calculator using the Identity model 
from Bio.Phylo.TreeConstruction import DistanceCalculator 
calculator = DistanceCalculator('identity')


# Write the Distance Matrix 
distance_matrix = calculator.get_distance(alignment)
print(distance_matrix)

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor()

# Construct the phlyogenetic tree using UPGMA algorithm
UPGMATree = constructor.upgma(distance_matrix)

# Order the clades
UPGMATree.ladderize()

# Function to remove InnerXX labels
def remove_inner_labels(tree):
    for clade in tree.find_clades():
        if 'Inner' in clade.name:
            clade.name = ''

# Remove InnerXX labels
remove_inner_labels(UPGMATree)

# Make a better looking tree using the features of matplotlib 
fig_size_y = 17
fig_size_x = 24

fig = plt.figure(figsize=(fig_size_x, fig_size_y), dpi=150) # create figure & set the size 
plt.rc('font', size=22)              # fontsize of the leaf and node labels 
plt.rc('xtick', labelsize=18)       # fontsize of the tick labels
plt.rc('ytick', labelsize=18)       # fontsize of the tick labels
axes = fig.add_subplot(1, 1, 1)

# drawing the tree
Phylo.draw(UPGMATree, axes=axes)
fig.savefig("./python/images/UPGMATree", dpi=200)

# Construct the phlyogenetic tree using NJ algorithm
NJTree = constructor.nj(distance_matrix)

# Order the clades
NJTree.ladderize()

# Remove InnerXX labels
remove_inner_labels(NJTree)

# Create a new figure for the NJ tree
fig2 = plt.figure(figsize=(fig_size_x, fig_size_y), dpi=150)
axes2 = fig2.add_subplot(1, 1, 1)
Phylo.draw(NJTree, axes=axes2)

# Save the figure
fig2.savefig("./python/images/NJTree", dpi=200)


