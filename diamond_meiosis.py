#!/usr/bin/env python3

import os
import sys

def grab_ogs(fas, taxon):
    dmnd_cmd = f'diamond blastp -k 0 -d HookDB -q {fas} -e 1e-10 -f 6 '\
        f'-o {taxon}.Meiosis_OGs.tsv'
    os.system(dmnd_cmd)
    return f'{taxon}.Meiosis_OGs.tsv'

def parse_table(tsv, taxon):
    ogs = {}
    for i in open(tsv).readlines():
        ogs.setdefault((i.split('\t')[0].split('_')[0]+"_"+i.split('\t')[0].split('_')[1]),[]).append(
            i.split('\t')[1].split('_')[-1])

    with open(f'{taxon}.Meiotic_Genes_to_OGs.tsv','w+') as w:
        w.write('Gene Name\tOG5\n')
        for k, v in ogs.items():
            for o in list(set(v)):
                w.write(f'{k}\t{o}\n')


if __name__ == '__main__':
    if len(sys.argv[1:]) != 2:
        print('Usage:\n\n   python diamond_meiosis.py FASTAwithGenes.fas Taxon\n\n')
        sys.exit(1)
    fas = sys.argv[1]
    taxon = sys.argv[2]
    tsv = grab_ogs(fas, taxon)
    parse_table(tsv, taxon)
