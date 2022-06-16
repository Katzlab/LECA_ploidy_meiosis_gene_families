# Usage: bash eggnog.sh <gene family name> <kog/cog number>

import pandas as pd
import os
import sys

def main(id_list):
	ids = pd.read_csv(id_list, sep=",", header=0)
	for ind in ids.index:
		if not str(ids['cog'][ind]) == "nan":
			os.system("curl --compressed http://eggnogapi5.embl.de/nog_data/text/fasta/" + ids['cog'][ind] + " --output " + ids['gene'][ind]+"_"+ids['cog'][ind] + ".fasta")
			os.system("mv " + ids['gene'][ind]+"_"+ids['cog'][ind] + ".fasta" + " original")
		if not str(ids['kog'][ind]) == "nan":
			os.system("curl --compressed http://eggnogapi5.embl.de/nog_data/text/fasta/" + ids['kog'][ind] + " --output " + ids['gene'][ind]+"_"+ids['kog'][ind] + ".fasta")
			os.system("mv " + ids['gene'][ind]+"_"+ids['kog'][ind] + ".fasta" + " original")
		if not str(ids['arcog'][ind]) == "nan":
			os.system("curl --compressed http://eggnogapi5.embl.de/nog_data/text/fasta/" + ids['arcog'][ind] + " --output " + ids['gene'][ind]+"_"+ids['arcog'][ind] + ".fasta")
			os.system("mv " + ids['gene'][ind]+"_"+ids['arcog'][ind] + ".fasta" + " original")
		elif not str(ids['enog'][ind]) == "nan":
			os.system("curl --compressed http://eggnogapi5.embl.de/nog_data/text/fasta/" + ids['enog'][ind] + " --output " + ids['gene'][ind]+"_"+ids['enog'][ind] + ".fasta")
			os.system("mv " + ids['gene'][ind]+"_"+ids['enog'][ind] + ".fasta" + " original")
		else :
			print(ids['gene'][ind] + " is missing an ID.")
			
if __name__ == '__main__':
	id_list = sys.argv[1]
	main(id_list)



	

