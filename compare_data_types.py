#! /usr/bin/env/ python3

import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import itertools
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

##########################################################################################
## METHODS
##########################################################################################
def load_tabular_as_dict(file):
	codon2fts_dict = defaultdict(list)
	with open(file, 'rt') as f:
		records = f.read().splitlines()
		for record in records:
			fields = record.split("\t")
			codon2fts_dict[fields[0]].append(fields[1:])
	return codon2fts_dict

def add_tRNA_data(annotated_prots_dict, codon2trna_dict):
	codon2trna_dict = {codon: float(codon2trna_dict[codon][0][0]) for codon in codon2trna_dict}
	for protein in annotated_prots_dict:
		for record in annotated_prots_dict[protein]:
			record.append(codon2trna_dict[record[0]])

def get_tRNA_avg_prof(data, feature):
	data_before_init = []
	data_after_init = []

	for protein in data:
		if feature in str(data[protein]):
			before_init = []
			for record in data[protein]:
				if feature not in record[-2]:
					before_init.append(record[-1])
				else:
					break
			init_to_end	= [x[-1] for x in data[protein][len(before_init):]]

			data_before_init.append(before_init)
			data_after_init.append(init_to_end)

	if feature != 'SIGNAL' and feature != 'INIT_MET':
		before_init_filled = np.array(list(itertools.zip_longest(*[x[::-1] for x in data_before_init], fillvalue=np.nan))).astype('float64')
		init_to_end_filled = np.array(list(itertools.zip_longest(*data_after_init, fillvalue=np.nan))).astype('float64')	
		average_init_to_end = np.nanmean(init_to_end_filled, axis=1)
		average_before_init = np.nanmean(before_init_filled, axis=1)[::-1]
		average_profile = np.concatenate((average_before_init, average_init_to_end), axis=None)
		number_of_profiles = before_init_filled.shape[1]
		point_0 = len(average_before_init)
			
	else:
		init_to_end_filled = np.array(list(itertools.zip_longest(*data_after_init, fillvalue=np.nan))).astype('float64')	
		average_init_to_end = np.nanmean(init_to_end_filled, axis=1)
		average_profile = average_init_to_end
		number_of_profiles = init_to_end_filled.shape[1]
		point_0 = 0

	average_profile = average_profile[point_0-51:point_0+50]

	return average_profile

def obtain_average_of_profiles(RSCU_profiles, codon_type_profiles):
	rscu_avg_profile = []
	rpcu_profile = []

	rscu_avg_profile = [
		np.nanmean(np.array([x[0].split(",")
		for x in RSCU_profiles['coli']]).astype('float64'), axis=0)
	]

	rpcu_profile = [
		np.nanmean(np.array([x[0].split(",")
		for x in codon_type_profiles['coli']]).astype('float64'), axis=0)
	]


	return rscu_avg_profile, rpcu_profile

def visualize(trna_avg_prof, rscu_avg_prof, rpcu_prof, feature, smoothing_window, results_path):
		profiles_data = [trna_avg_prof, rscu_avg_prof, rpcu_prof]
		smoothed_profiles = []
		for profile in profiles_data:
			cumsum_vec = np.cumsum(np.insert(profile, 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			smoothed_profiles.append(smoothing_vec)	
			
		smoothing_positions = [x for x in range(smoothing_window//2, len(smoothed_profiles[0])+smoothing_window//2)]	

		host = host_subplot(111, axes_class=AA.Axes)
		plt.subplots_adjust(right=0.75)

		par1 = host.twinx()
		par2 = host.twinx()

		offset = 60
		new_fixed_axis = par2.get_grid_helper().new_fixed_axis
		par2.axis["right"] = new_fixed_axis(loc="right", axes=par2,
        		                                offset=(offset, 0))

		par1.axis["right"].toggle(all=True)
		par2.axis["right"].toggle(all=True)

		host.set_ylim(2.9, 3.5)
		par1.set_ylim(1.125, 1.3)
		par2.set_ylim(0.42, 0.6)

		host.set_xlabel("Aa position", labelpad=10)
		host.set_ylabel("average [tRNA]")
		par1.set_ylabel("average RSCU")
		par2.set_ylabel("RPCU")

		p1, = host.plot(smoothing_positions, smoothed_profiles[0], label="average [tRNA]")
		p2, = par1.plot(smoothing_positions, smoothed_profiles[1], label="average RSCU")
		p3, = par2.plot(smoothing_positions, smoothed_profiles[2], label="RPCU")

		host.axvline(x=51, color='black', ls=':', lw=1)
		plt.xticks([10, 30, 51, 70, 90], ['-40', '-20', '0', '20', '40'])
		plt.suptitle(feature)

		host.axis["left"].label.set_color(p1.get_color())
		par1.axis["right"].label.set_color(p2.get_color())
		par2.axis["right"].label.set_color(p3.get_color())

		plt.savefig(f'{results_path}/data_types_compared.svg')

		plt.close()


######################################################################################
## ARGPARSE
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("codon2fts_file", type=str, help="codon2fts file to be processed")
parser.add_argument("RSCU_profiles", type=str, help="profiles of RSCU values for the analyzed feature to be compared")
parser.add_argument("codon_type_profile", type=str, help="profiles of codon_type values for the analyzed feature to be compared")
parser.add_argument("codon2trna_file", type=str, help="codon_to_[tRNA] file to be processed")
parser.add_argument('output_path', nargs='?', type=str, help="desired path for writting output file", default=os.getcwd())
parser.add_argument('-f', '--feature', type=str, help='feature to be compared')
parser.add_argument('-w', '--window', type=int, help='smoothing window size (must be an even number)')
args = parser.parse_args()

######################################################################################
## MAIN
######################################################################################
results_path = os.path.join(args.output_path, f'results/coli')

annotated_proteins_dict = load_tabular_as_dict(args.codon2fts_file)
RSCU_prof = load_tabular_as_dict(args.RSCU_profiles)
codon_type_prof = load_tabular_as_dict(args.codon_type_profile)
codon2trna_dict = load_tabular_as_dict(args.codon2trna_file)

add_tRNA_data(annotated_proteins_dict, codon2trna_dict)
trna_avg_prof = get_tRNA_avg_prof(annotated_proteins_dict, args.feature)
rscu_avg_prof, rpcu_prof = obtain_average_of_profiles(RSCU_prof, codon_type_prof)
visualize(trna_avg_prof, rscu_avg_prof, rpcu_prof, args.feature, args.window, results_path)