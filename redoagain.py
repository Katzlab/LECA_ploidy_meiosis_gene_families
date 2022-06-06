from Bio import SeqIO
from random import sample

original_file = "/Users/katzlab_admin/Desktop/MAE/gex1bactseqtranslated.fasta"
corrected_file = "/Users/katzlab_admin/Desktop/MAE/gex1bactseqtranslatedsmallerseq.fasta"

with open(original_file) as original, open(corrected_file, 'w') as corrected:
    seqs = SeqIO.parse(original_file, "fasta")
    for seq in sample(list(seqs),500):
        SeqIO.write(seq, corrected, 'fasta')
