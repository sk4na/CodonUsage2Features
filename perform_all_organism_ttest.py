#! /usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy import stats

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

def obtain_average_of_profiles(all_profiles, all_random_profiles):
	organism_profiles = [profile for profile in all_profiles.values()]
	organism_random_profiles = [profile for profile in all_random_profiles.values()]
	profiles = []
	random_profiles	= []

	for organism in organism_profiles:
		for profile in organism:
			profiles.append(profile[0].split(","))
	profiles = np.array(profiles).astype('float64')

	for organism in organism_random_profiles:
		for profile in organism:
			random_profiles.append(profile[0].split(","))
	random_profiles = np.array(profiles).astype('float64')

	average_profiles = np.nanmean(profiles, axis=0)
	average_random_profiles = np.nanmean(random_profiles, axis=0)

	return average_profiles, average_random_profiles

def visualize(profile_data, feature, smoothing_window, data_type, results_path):	
		if feature != 'SIGNAL' and feature != 'INIT_MET':
				data = profile_data
				plt.plot(data)
				plt.axvline(x=51, color='black', ls=':', lw=1)
				plt.xticks([10, 30, 51, 70, 90], ['-40', '-20', '0', '20', '40'])
		else:
				data = profile_data
				plt.plot(data)		
				plt.axvline(x=0, color='black', ls=':', lw=1)
				plt.xticks([0, 20, 40, 60, 80], ['0', '20', '40', '60', '80'])

		if smoothing_window:
			cumsum_vec = np.cumsum(np.insert(data, 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			smoothing_positions = [x for x in range(smoothing_window//2, len(smoothing_vec)+smoothing_window//2)]
			plt.plot(smoothing_positions, smoothing_vec)

		plt.xlabel("Aa position", labelpad=10)
		plt.suptitle(feature)

		if data_type == 'RSCU':
			plt.ylabel("Average RSCU Values", labelpad=15)
			plt.savefig(f'{results_path}/graphs/{feature}_all_organisms_Average_RSCU.svg')
		elif data_type == 'codon_type':
			plt.ylabel("Relative Preferred Codon Usage", labelpad=15)
			plt.savefig(f'{results_path}/graphs/{feature}_all_organisms_RPCU.svg')

		plt.close()

def get_pair_wise_correlations(profiles):
	profiles_df = pd.DataFrame(profiles)
	profiles_corr = profiles_df.T.corr(method='pearson')
	unique_profiles_corr = profiles_corr.mask(np.tril(np.ones(profiles_corr.shape)).astype(bool))
	profiles_corr_dist = [value for value in unique_profiles_corr.stack().reset_index()[0]]

	return profiles_corr_dist


def perform_ttest(all_profiles, all_random_profiles, feature, data_type, results_path):
	all_organism_profiles = np.array([np.array(vec[0].split(',')).astype(float) for vec in [profile for profile in all_profiles.values()][0]])
	all_organism_random_profiles = np.array([np.array(vec[0].split(',')).astype(float) for vec in [profile for profile in all_random_profiles.values()][0]])

	if data_type == 'RSCU':
		profiles_corr_dist = get_pair_wise_correlations(all_organism_profiles)
		random_profiles_corr_dist = get_pair_wise_correlations(all_organism_random_profiles)

		pvalue = stats.ttest_ind(profiles_corr_dist, random_profiles_corr_dist, alternative='greater')[1]

	elif data_type == 'codon_type':
		window_avg_profiles = []
		window_avg_rand_profiles = []
		smoothing_window = 3

		for profile in all_organism_profiles:
			cumsum_vec = np.cumsum(np.insert(profile, 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			window_avg_profiles.append(smoothing_vec)

		for rand_profile in all_organism_random_profiles:			
			cumsum_vec = np.cumsum(np.insert(rand_profile, 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			window_avg_rand_profiles.append(smoothing_vec)
		
		window_prof_corr_dist = get_pair_wise_correlations(window_avg_profiles)
		window_random_prof_corr_dist = get_pair_wise_correlations(window_avg_rand_profiles)

		pvalue = stats.ttest_ind(window_prof_corr_dist, window_random_prof_corr_dist, alternative='greater')[1]

	if pvalue <= 0.05:
		with open(f'{results_path}/significant_profiles.txt', 'a+') as f:
			f.write(feature + "\t" f'significant with a p_value of: {pvalue}' + "\n")
					 			

######################################################################################
## ARGPARSE
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("all_profiles_file", type=str, help="File containing the profiles of each organism around a feature for a data_type")
parser.add_argument("all_random_profiles_file", type=str, help="File containing the random profiles of each organism for a data_type")
parser.add_argument('output_path', nargs='?', type=str, help="desired path for writting output file", default=os.getcwd())
parser.add_argument('-d', '--data_type', type=str, help='type of data wanted to be analyzed: RSCU or codon_type')
parser.add_argument('-f', '--feature', type=str, help='feature to be analyzed')
parser.add_argument('-w', '--window', type=int, help='smoothing window size (must be an even number)')
args = parser.parse_args()

######################################################################################
## MAIN
######################################################################################
profiles = load_tabular_as_dict(args.all_profiles_file)
random_profiles	= load_tabular_as_dict(args.all_random_profiles_file)
results_path = os.path.join(args.output_path, f'results/{args.data_type}_all_organism_profiles')
average_prof, average_random_prof = obtain_average_of_profiles(profiles, random_profiles)
perform_ttest(profiles, random_profiles, args.feature, args.data_type, results_path)
visualize(average_prof, args.feature, args.window, args.data_type, results_path)
