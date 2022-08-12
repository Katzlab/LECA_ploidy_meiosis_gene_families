#!/usr/bin/env python3
# coding=utf-8

#Author: Auden Cote-L'Heureux
#Contact at acotelheureux@smith.edu, audenemil@gmail.com, or on Slack if a Katzlabber.
#Last updated 08/17/21
#Notes: run with python3; pipeline written for LGT work (corresponding R script esp.) and may need tweaking for other applications

import os
import sys
import re
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from Bio.Blast.NCBIWWW import qblast
from CUB import *


def help():

	print('\nThis script is for building databases to compare the composition of two sets of sequences. It creates two files, gc3_master.csv and rscu_master.csv, measuring G+C content at third-position four-fold degenerate sites and the relative synonymous codon usage respectively. To run this script, give it a .fasta file of nucleotide sequences:\n\n\tpython CompositionComparison.py -i <.fasta file>\n')
	exit()


def get_args():

	if('-i' in sys.argv):
		try:
			seqs1 = sys.argv[sys.argv.index('-i') + 1]
		except IndexError:
			help()
	else:
		help()

	return seqs1



#seqs1 = os.listdir('GGjelly')


def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_temp.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 
		
	tree_corrected.write(tree2correct)
	tree_corrected.close()


def get_GC3s(seqs1):

	m = open('gc3_master.csv', 'a+')
	m.write('Seq,Type,GC3-Degen,ExpWrightENc,ObsWrightENc_6Fold,ObsWrightENc_No6Fold,ObsWeightedENc_6Fold,ObsWeightedENc_No6Fold\n')
		
	for type in [seqs1]:
		print('\nCompiling GC3 data for input file ' + type + '\n')
					
		os.system('python3 CUB.py ' + type + ' ' + type.split('.')[0] + '_comp universal')
		
		for l, line in enumerate(open(type.split('.')[0].split('/')[-1] + '_comp/Spreadsheets/' + type.split('.')[0].split('/')[-1] + '_comp.CompTrans.ENc.Raw.tsv')):
			line = [cell.strip() for cell in line.split('\t') if cell.strip() != '']
			if(l != 0):
				m.write(line[0] + ',' + type.split('.')[0].split('/')[-1] + ',' + ','.join(line[6:]) + '\n')
		
		
	m.close()
		
				
def make_plots(seqs1):

	taxa = list(dict.fromkeys([rec.description.strip()[:8] for rec in SeqIO.parse(seqs1, 'fasta')]))
	print(taxa)
	for taxon in taxa:
		os.system('/Library/Frameworks/R.framework/Versions/4.0/Resources/Rscript PlotComps.r ' + taxon)
		

def main():
		
	#seqs1 = get_args()
	
	#get_GC3s(seqs1)
	
	files = ["Am_ab_Mabd_filtered_NTD.fasta",
"Am_ab_Pesp_filtered_NTD.fasta",
"Am_ab_Relo_filtered_NTD.fasta",
"Am_ar_Mabd_filtered_NTD.fasta",
"Am_ar_Pesp_filtered_NTD.fasta",
"Am_ar_Relo_filtered_NTD.fasta",
"Am_db_Gfon_filtered_NTD.fasta",
"Am_db_Paes_filtered_NTD.fasta",
"Am_di_Acas_filtered_NTD.fasta",
"Am_di_Apyr_filtered_NTD.fasta",
"Am_di_Gfon_filtered_NTD.fasta",
"Am_di_Paes_filtered_NTD.fasta",
"Am_ib_Eexu_filtered_NTD.fasta",
"Am_is_Eexu_filtered_NTD.fasta",
"Am_mb_Capo_filtered_NTD.fasta",
"Am_mb_Cfur_filtered_NTD.fasta",
"Am_mb_Eoli_filtered_NTD.fasta",
"Am_mb_Pnoc_filtered_NTD.fasta",
"Am_my_Asub_filtered_NTD.fasta",
"Am_my_Capo_filtered_NTD.fasta",
"Am_my_Cfur_filtered_NTD.fasta",
"Am_my_Eoli_filtered_NTD.fasta",
"Am_my_Pnoc_filtered_NTD.fasta",
"Am_tb_As01_filtered_NTD.fasta",
"Am_tb_He21_filtered_NTD.fasta",
"Am_tb_Hele_filtered_NTD.fasta",
"Am_tb_Hp20_filtered_NTD.fasta",
"Am_tb_Hp92_filtered_NTD.fasta",
"Am_tb_Hpap_filtered_NTD.fasta",
"Am_tb_Nesp_filtered_NTD.fasta",
"Am_tb_Qs01_filtered_NTD.fasta",
"Am_tu_Aint_filtered_NTD.fasta",
"Am_tu_As01_filtered_NTD.fasta",
"Am_tu_He21_filtered_NTD.fasta",
"Am_tu_Hele_filtered_NTD.fasta",
"Am_tu_Hp20_filtered_NTD.fasta",
"Am_tu_Hp92_filtered_NTD.fasta",
"Am_tu_Hpap_filtered_NTD.fasta",
"Am_tu_Nesp_filtered_NTD.fasta",
"Am_tu_Qs01_filtered_NTD.fasta",
"Ba_ac_Cwoe_filtered_NTD.fasta",
"Ba_ac_Maur_filtered_NTD.fasta",
"Ba_ac_Nmul_filtered_NTD.fasta",
"Ba_ac_Rhsp_filtered_NTD.fasta",
"Ba_ac_Srap_filtered_NTD.fasta",
"Ba_ad_Bagg_filtered_NTD.fasta",
"Ba_ad_Cthe_filtered_NTD.fasta",
"Ba_ad_Gtun_filtered_NTD.fasta",
"Ba_ad_Hfoe_filtered_NTD.fasta",
"Ba_ad_Susi_filtered_NTD.fasta",
"Ba_aq_Aaeo_filtered_NTD.fasta",
"Ba_aq_Desp_filtered_NTD.fasta",
"Ba_aq_Pmar_filtered_NTD.fasta",
"Ba_aq_Ttak_filtered_NTD.fasta",
"Ba_Ba_Assp_filtered_NTD.fasta",
"Ba_Ba_Bcae_filtered_NTD.fasta",
"Ba_Ba_Ccha_filtered_NTD.fasta",
"Ba_Ba_Opac_filtered_NTD.fasta",
"Ba_Ba_Sgra_filtered_NTD.fasta",
"Ba_bc_Cchl_filtered_NTD.fasta",
"Ba_bc_Clut_filtered_NTD.fasta",
"Ba_bc_Ppha_filtered_NTD.fasta",
"Ba_bi_Ialb_filtered_NTD.fasta",
"Ba_ca_Cexi_filtered_NTD.fasta",
"Ba_cd_Cavi_filtered_NTD.fasta",
"Ba_cd_Nesp_filtered_NTD.fasta",
"Ba_ch_Caur_filtered_NTD.fasta",
"Ba_ch_Cbac_filtered_NTD.fasta",
"Ba_ch_Psub_filtered_NTD.fasta",
"Ba_ch_Tros_filtered_NTD.fasta",
"Ba_cv_Ccap_filtered_NTD.fasta",
"Ba_cv_Ppar_filtered_NTD.fasta",
"Ba_cv_Vesp_filtered_NTD.fasta",
"Ba_cy_Acyl_filtered_NTD.fasta",
"Ba_cy_Gvio_filtered_NTD.fasta",
"Ba_cy_Npun_filtered_NTD.fasta",
"Ba_cy_Smaj_filtered_NTD.fasta",
"Ba_de_Dswu_filtered_NTD.fasta",
"Ba_de_Mhyd_filtered_NTD.fasta",
"Ba_df_Ddes_filtered_NTD.fasta",
"Ba_df_Fsin_filtered_NTD.fasta",
"Ba_el_Emin_filtered_NTD.fasta",
"Ba_fc_Ccel_filtered_NTD.fasta",
"Ba_fu_Fnec_filtered_NTD.fasta",
"Ba_pa_Shen_filtered_NTD.fasta",
"Ba_pb_Pdur_filtered_NTD.fasta",
"Ba_pd_Hesp_filtered_NTD.fasta",
"Ba_pg_Paer_filtered_NTD.fasta",
"Ba_Pl_Gobs_filtered_NTD.fasta",
"Ba_Pl_Plsp_filtered_NTD.fasta",
"Ba_Pl_Ssal_filtered_NTD.fasta",
"Ba_sp_Bmay_filtered_NTD.fasta",
"Ba_sp_Lint_filtered_NTD.fasta",
"Ba_sp_Slut_filtered_NTD.fasta",
"Ba_sy_Sjon_filtered_NTD.fasta",
"Ba_ts_Tthe_filtered_NTD.fasta",
"EE_ab_Asig_filtered_NTD.fasta",
"EE_ab_Nlon_filtered_NTD.fasta",
"EE_ap_Asig_filtered_NTD.fasta",
"EE_ap_Nlon_filtered_NTD.fasta",
"EE_ap_Ttra_filtered_NTD.fasta",
"EE_cb_Cpar_filtered_NTD.fasta",
"EE_cb_Rlen_filtered_NTD.fasta",
"EE_ce_Acsp_filtered_NTD.fasta",
"EE_ce_Rhet_filtered_NTD.fasta",
"EE_cr_Cpar_filtered_NTD.fasta",
"EE_cr_Gthe_filtered_NTD.fasta",
"EE_cr_Rlen_filtered_NTD.fasta",
"EE_ha_Cpol_filtered_NTD.fasta",
"EE_ha_Ispa_filtered_NTD.fasta",
"EE_ha_Ppar_filtered_NTD.fasta",
"EE_ha_Saps_filtered_NTD.fasta",
"EE_hb_Cpol_filtered_NTD.fasta",
"EE_hb_Ispa_filtered_NTD.fasta",
"EE_hb_Ppar_filtered_NTD.fasta",
"EE_hb_Saps_filtered_NTD.fasta",
"EE_ib_Hkuk_filtered_NTD.fasta",
"EE_ib_Smul_filtered_NTD.fasta",
"EE_is_Hkuk_filtered_NTD.fasta",
"EE_is_Smul_filtered_NTD.fasta",
"EE_is_Tsub_filtered_NTD.fasta",
"Ex_eb_Egym_filtered_NTD.fasta",
"Ex_eu_Bsal_filtered_NTD.fasta",
"Ex_eu_Egym_filtered_NTD.fasta",
"Ex_eu_Pemi_filtered_NTD.fasta",
"Ex_eu_Scul_filtered_NTD.fasta",
"Ex_fo_Ssal_filtered_NTD.fasta",
"Ex_hb_Pcos_filtered_NTD.fasta",
"Ex_hb_Pkir_filtered_NTD.fasta",
"Ex_he_Ngru_filtered_NTD.fasta",
"Ex_he_Pcos_filtered_NTD.fasta",
"Ex_he_Pkir_filtered_NTD.fasta",
"Ex_ja_Agod_filtered_NTD.fasta",
"Ex_ja_Sinc_filtered_NTD.fasta",
"Ex_jb_Sinc_filtered_NTD.fasta",
"Ex_ma_Goke_filtered_NTD.fasta",
"Ex_ma_Mjak_filtered_NTD.fasta",
"Ex_mb_Goke_filtered_NTD.fasta",
"Ex_mb_Mjak_filtered_NTD.fasta",
"Ex_ob_Mono_filtered_NTD.fasta",
"Ex_ox_Mono_filtered_NTD.fasta",
"Ex_pa_Tfoe_filtered_NTD.fasta",
"Op_fu_Bden_filtered_NTD.fasta",
"Op_fu_Ccor_filtered_NTD.fasta",
"Op_fu_Dpri_filtered_NTD.fasta",
"Op_fu_Gpro_filtered_NTD.fasta",
"Op_fu_Lcor_filtered_NTD.fasta",
"Op_fu_Mmel_filtered_NTD.fasta",
"Op_fu_Wmel_filtered_NTD.fasta",
"Op_me_Ctel_filtered_NTD.fasta",
"Op_me_Dpul_filtered_NTD.fasta",
"Op_me_Hrob_filtered_NTD.fasta",
"Op_me_Odio_filtered_NTD.fasta",
"Op_me_Skow_filtered_NTD.fasta",
"Pl_gb_Glwi_filtered_NTD.fasta",
"Pl_gl_Glwi_filtered_NTD.fasta",
"Pl_gr_Atri_filtered_NTD.fasta",
"Pl_gr_Dsal_filtered_NTD.fasta",
"Pl_gr_Egui_filtered_NTD.fasta",
"Pl_gr_Hann_filtered_NTD.fasta",
"Pl_gr_Heli_filtered_NTD.fasta",
"Pl_gr_Knit_filtered_NTD.fasta",
"Pl_gr_Mpol_filtered_NTD.fasta",
"Pl_gr_Oluc_filtered_NTD.fasta",
"Pl_rb_Ccoe_filtered_NTD.fasta",
"Pl_rb_Eden_filtered_NTD.fasta",
"Pl_rb_Toli_filtered_NTD.fasta",
"Pl_rh_Ccoe_filtered_NTD.fasta",
"Pl_rh_Eden_filtered_NTD.fasta",
"Pl_rh_Gsul_filtered_NTD.fasta",
"Pl_rh_Toli_filtered_NTD.fasta",
"Sr_ab_Magi_filtered_NTD.fasta",
"Sr_ap_Emax_filtered_NTD.fasta",
"Sr_ap_Magi_filtered_NTD.fasta",
"Sr_cb_Ba02_filtered_NTD.fasta",
"Sr_cb_Bt03_filtered_NTD.fasta",
"Sr_cb_Casp_filtered_NTD.fasta",
"Sr_cb_Cu03_filtered_NTD.fasta",
"Sr_cb_Dm01_filtered_NTD.fasta",
"Sr_cb_Dn01_filtered_NTD.fasta",
"Sr_cb_Gs02_filtered_NTD.fasta",
"Sr_cb_Ls03_filtered_NTD.fasta",
"Sr_cb_Lx03_filtered_NTD.fasta",
"Sr_cb_Mstr_filtered_NTD.fasta",
"Sr_cb_Nvar_filtered_NTD.fasta",
"Sr_cb_Rx04_filtered_NTD.fasta",
"Sr_cb_Sa01_filtered_NTD.fasta",
"Sr_cb_TR03_filtered_NTD.fasta",
"Sr_ch_Vbra_filtered_NTD.fasta",
"Sr_ci_Ba02_filtered_NTD.fasta",
"Sr_ci_Bt03_filtered_NTD.fasta",
"Sr_ci_Casp_filtered_NTD.fasta",
"Sr_ci_Cu03_filtered_NTD.fasta",
"Sr_ci_Dm01_filtered_NTD.fasta",
"Sr_ci_Dn01_filtered_NTD.fasta",
"Sr_ci_Gs02_filtered_NTD.fasta",
"Sr_ci_Imul_filtered_NTD.fasta",
"Sr_ci_Ls03_filtered_NTD.fasta",
"Sr_ci_Lx03_filtered_NTD.fasta",
"Sr_ci_Mstr_filtered_NTD.fasta",
"Sr_ci_Nvar_filtered_NTD.fasta",
"Sr_ci_Otri_filtered_NTD.fasta",
"Sr_ci_Pper_filtered_NTD.fasta",
"Sr_ci_Ptet_filtered_NTD.fasta",
"Sr_ci_Rx04_filtered_NTD.fasta",
"Sr_ci_Sa01_filtered_NTD.fasta",
"Sr_ci_Slem_filtered_NTD.fasta",
"Sr_ci_TR03_filtered_NTD.fasta",
"Sr_db_Aspi_filtered_NTD.fasta",
"Sr_db_Gaus_filtered_NTD.fasta",
"Sr_db_Gspi_filtered_NTD.fasta",
"Sr_db_Kven_filtered_NTD.fasta",
"Sr_db_Nsci_filtered_NTD.fasta",
"Sr_db_Shan_filtered_NTD.fasta",
"Sr_di_Amas_filtered_NTD.fasta",
"Sr_di_Aspi_filtered_NTD.fasta",
"Sr_di_Gaus_filtered_NTD.fasta",
"Sr_di_Gspi_filtered_NTD.fasta",
"Sr_di_Kven_filtered_NTD.fasta",
"Sr_di_Nsci_filtered_NTD.fasta",
"Sr_di_Shan_filtered_NTD.fasta",
"Sr_di_Smic_filtered_NTD.fasta",
"Sr_pb_Olen_filtered_NTD.fasta",
"Sr_pb_Pche_filtered_NTD.fasta",
"Sr_pe_Olen_filtered_NTD.fasta",
"Sr_pe_Pche_filtered_NTD.fasta",
"Sr_pe_Perk_filtered_NTD.fasta",
"Sr_rb_As02_filtered_NTD.fasta",
"Sr_rb_As05_filtered_NTD.fasta",
"Sr_rb_AspA_filtered_NTD.fasta",
"Sr_rb_Blon_filtered_NTD.fasta",
"Sr_rb_Crep_filtered_NTD.fasta",
"Sr_rb_Emar_filtered_NTD.fasta",
"Sr_rb_HgpA_filtered_NTD.fasta",
"Sr_rb_Lamo_filtered_NTD.fasta",
"Sr_rb_Loce_filtered_NTD.fasta",
"Sr_rb_Mf02_filtered_NTD.fasta",
"Sr_rb_Pglo_filtered_NTD.fasta",
"Sr_rb_PspA_filtered_NTD.fasta",
"Sr_rb_Qj01_filtered_NTD.fasta",
"Sr_rh_Aa01_filtered_NTD.fasta",
"Sr_rh_As02_filtered_NTD.fasta",
"Sr_rh_As05_filtered_NTD.fasta",
"Sr_rh_AspA_filtered_NTD.fasta",
"Sr_rh_Blon_filtered_NTD.fasta",
"Sr_rh_Crep_filtered_NTD.fasta",
"Sr_rh_Emar_filtered_NTD.fasta",
"Sr_rh_HgpA_filtered_NTD.fasta",
"Sr_rh_Lamo_filtered_NTD.fasta",
"Sr_rh_Loce_filtered_NTD.fasta",
"Sr_rh_Mf02_filtered_NTD.fasta",
"Sr_rh_Pglo_filtered_NTD.fasta",
"Sr_rh_Pmic_filtered_NTD.fasta",
"Sr_rh_PspA_filtered_NTD.fasta",
"Sr_rh_Qj01_filtered_NTD.fasta",
"Sr_sb_Atth_filtered_NTD.fasta",
"Sr_sb_Caro_filtered_NTD.fasta",
"Sr_sb_Croe_filtered_NTD.fasta",
"Sr_sb_Cspa_filtered_NTD.fasta",
"Sr_sb_Espi_filtered_NTD.fasta",
"Sr_sb_Fspa_filtered_NTD.fasta",
"Sr_sb_Goce_filtered_NTD.fasta",
"Sr_sb_Hseo_filtered_NTD.fasta",
"Sr_sb_Ppar_filtered_NTD.fasta",
"Sr_sb_Svul_filtered_NTD.fasta",
"Sr_sb_Trot_filtered_NTD.fasta",
"Sr_st_Aana_filtered_NTD.fasta",
"Sr_st_Aast_filtered_NTD.fasta",
"Sr_st_Acof_filtered_NTD.fasta",
"Sr_st_Ainv_filtered_NTD.fasta",
"Sr_st_Alim_filtered_NTD.fasta",
"Sr_st_Apal_filtered_NTD.fasta",
"Sr_st_Atth_filtered_NTD.fasta",
"Sr_st_Bhom_filtered_NTD.fasta",
"Sr_st_Caro_filtered_NTD.fasta",
"Sr_st_Cneo_filtered_NTD.fasta",
"Sr_st_Croe_filtered_NTD.fasta",
"Sr_st_Cspa_filtered_NTD.fasta",
"Sr_st_Esil_filtered_NTD.fasta",
"Sr_st_Espi_filtered_NTD.fasta",
"Sr_st_Fspa_filtered_NTD.fasta",
"Sr_st_Goce_filtered_NTD.fasta",
"Sr_st_Hseo_filtered_NTD.fasta",
"Sr_st_Ngad_filtered_NTD.fasta",
"Sr_st_Pinf_filtered_NTD.fasta",
"Sr_st_Pins_filtered_NTD.fasta",
"Sr_st_Ppar_filtered_NTD.fasta",
"Sr_st_Ptri_filtered_NTD.fasta",
"Sr_st_Sdic_filtered_NTD.fasta",
"Sr_st_Spar_filtered_NTD.fasta",
"Sr_st_Svul_filtered_NTD.fasta",
"Sr_st_Trot_filtered_NTD.fasta",
"Za_as_Heia_filtered_NTD.fasta",
"Za_as_Heib_filtered_NTD.fasta",
"Za_as_Loka_filtered_NTD.fasta",
"Za_as_Loki_filtered_NTD.fasta",
"Za_as_Odin_filtered_NTD.fasta",
"Za_as_Thob_filtered_NTD.fasta",
"Za_as_Thor_filtered_NTD.fasta",
"Za_Ba_Crea_filtered_NTD.fasta",
"Za_cr_Aper_filtered_NTD.fasta",
"Za_cr_Asac_filtered_NTD.fasta",
"Za_cr_Clag_filtered_NTD.fasta",
"Za_cr_Ffon_filtered_NTD.fasta",
"Za_cr_Paer_filtered_NTD.fasta",
"Za_cr_Sisl_filtered_NTD.fasta",
"Za_cr_Tpen_filtered_NTD.fasta",
"Za_cr_Vdis_filtered_NTD.fasta",
"Za_ea_Aven_filtered_NTD.fasta",
"Za_ea_Fpla_filtered_NTD.fasta",
"Za_eb_Mrum_filtered_NTD.fasta",
"Za_eb_Msta_filtered_NTD.fasta",
"Za_ec_Maeo_filtered_NTD.fasta",
"Za_ec_Minf_filtered_NTD.fasta",
"Za_eh_Hcar_filtered_NTD.fasta",
"Za_eh_Hrub_filtered_NTD.fasta",
"Za_eh_Hvol_filtered_NTD.fasta",
"Za_eh_Hxan_filtered_NTD.fasta",
"Za_eh_Nmag_filtered_NTD.fasta",
"Za_em_Mlim_filtered_NTD.fasta",
"Za_em_Mlum_filtered_NTD.fasta",
"Za_em_Mpsy_filtered_NTD.fasta",
"Za_ep_Taci_filtered_NTD.fasta",
"Za_et_Tkod_filtered_NTD.fasta",
"Za_eu_Aspa_filtered_NTD.fasta",
"Za_ey_Mkan_filtered_NTD.fasta"
]
	for item in files:
		seqs1 = item
		get_GC3s(seqs1)
	
#	make_plots(seqs1)
		
	
main()
	