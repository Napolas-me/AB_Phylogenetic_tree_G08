from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import os 

# Concatenate all fasta files in a directory

directory = './Species/data/'
output = './Species/concatenated.fasta'

#Create output file
with open(output, 'w') as outfile:
    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            #read the fasta file
            records = list(SeqIO.parse(directory + filename, "fasta"))
            #add ids to the sequences
            for record in records:
                record.id = filename.split('(')[0] 
                #print id and description
                print(record.id, record.description)
            for record in records:
                SeqIO.write(record, outfile, "fasta")
                

        