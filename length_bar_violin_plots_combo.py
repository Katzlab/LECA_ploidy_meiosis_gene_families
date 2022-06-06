"""Produce size tsvs but no plots"""
# Katzlab SURF 2022
# Created by Emma Schumacher 05/17/22

#Import statements...
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from random import sample
from pathlib import Path
import os
import shutil


#used to get size of sample (20 % unless the sample is <100 seq)
def get_sampsize(og_fasta_file):
    #count up
    i = 0
    with open(og_fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                i += 1
    #get size 
    if (i > 100):
        s_size = int(i * 0.2)
    else:
        s_size = i 
    return (s_size)
    

#used to randomly sample 500 OGs from fasta to make stuff run better
def get_random(og_fasta_file, size, multi = None):
    #ogpath = Path(og_fasta_file).parent.absolute()
    #print(ogpath)
        
    if not multi:
        name = og_fasta_file[0:og_fasta_file.index(".fasta")] + "smallerseq.fasta"
    else:
        name = og_fasta_file[0:og_fasta_file.index(".fasta")] + str(multi) + "smallerseq.fasta" 

    #open(name, "w") #str(Path(og_fasta_file).parent) + name
    corrected_file = name
    
    with open(og_fasta_file, "r") as original, open(corrected_file, 'w') as corrected:
        seqs = SeqIO.parse(og_fasta_file, "fasta")
        
        for seq in sample(list(seqs), size):
            SeqIO.write(seq, corrected, 'fasta')
            
    return(corrected_file)

# Provide a FASTA file of sequences, along with an OG number, to make a quick
# table of sequence lengths (note this is per OG).
def get_sizes(fasta_file, og_number): 
    # Snag the sequence names and their lengths from a given FASTA file.
    # Note that you can convert the dictionary into a dataframe, but that is
    # more of a pain...
    seq_sizes = {i.id:len(i.seq) for i in SeqIO.parse(fasta_file,'fasta')}
        
    # Save the data as a spreadsheet with Tab-Separated-Values (TSV)
    with open(f'{og_number}.SeqLength.tsv','w+') as w:
        w.write('OG\tSequence_Name\tLength\tDomain\n')

        for k, v in seq_sizes.items():
            # Add in the taxonomic "domains" when possible! This is mostly for
            # fancy plotting. 
            if k[:2] == 'Ba':
                w.write(f'{og_number}\t{k}\t{v}\tBacteria\n')

            elif k[:2] == 'Za':
                w.write(f'{og_number}\t{k}\t{v}\tArchaea\n')

            else:
                w.write(f'{og_number}\t{k}\t{v}\tEukaryota\n')
        return(og_number + ".SeqLength.tsv")

#used to call size tsvs for multiple files 
def make_multiple_OGs_sizes(fasta_files, og_number, runs):
    all_OG_tsvs = []
    if (runs == 1):
        # Calls function to make a tsv for next function
        get_sizes(fasta_files.pop(), og_number)
        all_OG_tsvs.append(og_number + ".SeqLength.tsv")
    else:
        for i in range(0,runs):
            get_sizes(fasta_files.pop(), (og_number + "run_" + str(i)))
            all_OG_tsvs.append(og_number + "run_" + str(i) + ".SeqLength.tsv")
    return (all_OG_tsvs)

#used to call multiple random samples, great for seeing if distributions are consistent
def make_multiple_randoms(file, runs, size):
    all_samples = []
    if (runs == 1):
        # Calls function to make a tsv for next function
        fn = get_random(file, size)
        all_samples.append(fn)
        
    else:
        for i in range(0, runs):
            fn = get_random(file, size, multi = i)
            all_samples.append(fn)
    return (all_samples) 
 

 
if __name__ == '__main__':
    ## EDIT HERE: Modify this so it is the path your fasta file
    og_fasta_file = "/Users/emmaschumacher/Desktop/Testnewscript/gettsv/combined_mre11taxid.fasta"
    ## EDIT HERE: Modify this so it is whatever label you want on your X axis 
    og_number = "MRE11"
    ## EDIT HERE: However many random samples you want (if none, 0)
    runs = 1


    
    if (runs == 0):
        get_sizes(og_fasta_file, og_number) 
    
    else:
        #gets sample size
        size = get_sampsize(og_fasta_file)
        
        #gets a random sample
        fasta_files = make_multiple_randoms(og_fasta_file, runs, size)
        
        # Calls function to make a tsv for next function
        all_OG_tsvs = make_multiple_OGs_sizes(fasta_files, og_number, runs) 









    
