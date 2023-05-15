#! /usr/bin/env python3

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

##########################################################################################
## METHODS
##########################################################################################
def load_tabular_as_dict(file):
	data_dict = defaultdict(list)
	with open(file, 'rt') as f:
		records = f.read().splitlines()
		for record in records:
			fields = record.split("\t")
			data_dict[fields[0]].append(fields[1:])
	return data_dict

def obtain_average_of_profiles(all_profiles):
	organism_avg_profiles = {
		organism: np.nanmean(np.array([x[0].split(",")
		for x in all_profiles[organism]]).astype('float64'), axis=0)
		for organism in all_profiles.keys()
	}
	return organism_avg_profiles

def visualize(profile_data, feature, smoothing_window, data_type, results_path):
		smoothed_dict = {}
		for organism in profile_data:
			cumsum_vec = np.cumsum(np.insert(profile_data[organism], 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			smoothed_dict[organism] = smoothing_vec	
			
		smoothing_positions = [x for x in range(smoothing_window//2, len(smoothed_dict['coli'])+smoothing_window//2)]	
		if feature != 'SIGNAL' and feature != 'INIT_MET':
				plt.plot(smoothing_positions, smoothed_dict['coli'])
				plt.plot(smoothing_positions, smoothed_dict['c.elegans'])
				plt.plot(smoothing_positions, smoothed_dict['drosophila'])
				plt.plot(smoothing_positions, smoothed_dict['human'])
				plt.plot(smoothing_positions, smoothed_dict['yeast'])
				plt.axvline(x=51, color='black', ls=':', lw=1)
				plt.xticks([10, 30, 51, 70, 90], ['-40', '-20', '0', '20', '40'])
		else:
				plt.plot(smoothing_positions, smoothed_dict['coli'])
				plt.plot(smoothing_positions, smoothed_dict['c.elegans'])
				plt.plot(smoothing_positions, smoothed_dict['drosophila'])
				plt.plot(smoothing_positions, smoothed_dict['human'])
				plt.plot(smoothing_positions, smoothed_dict['yeast'])		
				plt.axvline(x=0, color='black', ls=':', lw=1)
				plt.xticks([0, 20, 40, 60, 80], ['0', '20', '40', '60', '80'])

		plt.xlabel("Aa position", labelpad=10)
		plt.suptitle(feature)
		plt.legend(['E. coli', 'C. elegans', 'D. melanogaster', 'H. sapiens', 'S. cerevisiae'], loc='center left', bbox_to_anchor=(1, 0.5))
		plt.tight_layout()
		if data_type == 'RSCU':
			plt.ylabel("Average RSCU Values", labelpad=15)
			plt.savefig(f'{results_path}/graphs/{feature}_comparison_Average_RSCU.svg')
		elif data_type == 'codon_type':
			plt.ylabel("Relative Preferred Codon Usage", labelpad=15)
			plt.savefig(f'{results_path}/graphs/{feature}_comparison_RPCU.svg')

		plt.close()

######################################################################################
## ARGPARSE
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("all_profiles_file", type=str, help="File containing the profiles of each organism around a feature for a data_type")
parser.add_argument('output_path', nargs='?', type=str, help="desired path for writting output file", default=os.getcwd())
parser.add_argument('-d', '--data_type', type=str, help='type of data wanted to be analyzed: RSCU or codon_type')
parser.add_argument('-f', '--feature', type=str, help='feature to be analyzed')
parser.add_argument('-w', '--window', type=int, help='smoothing window size (must be an even number)')
args = parser.parse_args()

######################################################################################
## MAIN
######################################################################################
profiles = load_tabular_as_dict(args.all_profiles_file)
results_path = os.path.join(args.output_path, f'results/{args.data_type}_all_organism_profiles')
organism_avg_profiles = obtain_average_of_profiles(profiles)
visualize(organism_avg_profiles, args.feature, args.window, args.data_type, results_path)