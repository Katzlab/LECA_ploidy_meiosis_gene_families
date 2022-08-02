# AUTHOR: Caitlin Timmons and Emma Schumacher
# Purpose: BLAST against Katzlab hook database to identify best hit OGs for a protein sequence
# Usage: python diamond_automate.py <folder with diamond_meiosis.py script and HookDB file> <folder of input FASTAs>

import os
import sys

#does the translation
def main(wd, folder, fasta):    
	
	# USER input: 
	os.chdir(wd)
	
	# makes folders to store the diamond output
	if not os.path.exists('diamond_output'):
		os.makedirs('diamond_output') 
	
	os.system('touch ' + fasta.split(".")[0]+".Meiosis_OGs.tsv") #make output file
	
	os.system('python diamond_meiosis.py ' + folder+"/"+fasta + " " + fasta.split(".")[0]) 
	# deletes cluster file (we don't use)
	os.system('mv *Meiosis_OGs.tsv full_hits')
	
	os.system('mv *Meiotic_Genes_to_OGs.tsv diamond_output')


def automate(wd, folder):
	tranfiles = []
	path = os.getcwd() + '/' + folder
	for file in os.listdir(path):
		if not file.startswith('.'):
			if file.endswith("fasta"):  
				tranfiles.append(file)
	
	for thing in tranfiles:  
		main(wd, folder, thing)
	


if __name__ == '__main__':
	wd = sys.argv[1]
	input_folder = sys.argv[2]
	automate(wd, input_folder)