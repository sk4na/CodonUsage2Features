#! /usr/bin/env python3

import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
from collections import defaultdict
import random
from scipy import stats

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

def create_random_profiles(data, data_type, results_path, organism):
	profile_dict = defaultdict(list)
	random_segments = []
	for protein in data:
		if data_type == 'RSCU':
			profile_dict[protein] = [float(x[4]) for x in codon2fts[protein]]
		elif data_type == 'codon_type':
			profile_dict[protein] = [1.0 if x[3] == 'TRUE' else 0.0 for x in codon2fts[protein]]

	for i in range(0, 2000):
		prot_len = 0
		while prot_len < 101:
			random_protein = random.choice(list(profile_dict))
			prot_len = len(profile_dict[random_protein])

		random_position = random.choice(range(0, len(profile_dict[random_protein])))
		seed_value = np.array([profile_dict[random_protein][random_position]])
		segment_before_randpos_len = len(profile_dict[random_protein][:random_position])
		segment_after_randpos_len = len(profile_dict[random_protein][random_position+1:])

		if segment_before_randpos_len != 0 and segment_after_randpos_len != 0:
			if segment_after_randpos_len >= 50 and segment_before_randpos_len >= 50:
				values_before_seed = np.array(profile_dict[random_protein][segment_before_randpos_len-50:random_position])
				values_after_seed = np.array(profile_dict[random_protein][random_position+1:random_position+51])
				random_segment = np.concatenate((values_before_seed, seed_value, values_after_seed), axis=None)
			elif segment_after_randpos_len >= 50 and segment_before_randpos_len < 50:
				values_before_seed = np.array(profile_dict[random_protein][:random_position])
				values_after_seed = np.array(profile_dict[random_protein][random_position+1:random_position+(101-len(values_before_seed))])
				random_segment = np.concatenate((values_before_seed, seed_value, values_after_seed), axis=None)
			elif segment_after_randpos_len < 50 and segment_before_randpos_len >= 50:
				values_after_seed = np.array(profile_dict[random_protein][random_position+1:])
				values_before_seed = np.array(profile_dict[random_protein][random_position-(100-len(values_after_seed)):random_position])
				random_segment = np.concatenate((values_before_seed, seed_value, values_after_seed), axis=None)
		elif segment_after_randpos_len == 0:
			random_segment = np.array(profile_dict[random_protein][random_position-100:])
		elif segment_before_randpos_len == 0:
			random_segment = np.array(profile_dict[random_protein][:101])

		random_segments.append(random_segment)

	with open(f'{results_path}/{data_type}/random_profiles.txt', 'w') as f:
		for i in range(0, len(random_segments)):
			f.write(organism + "\t" + (',').join(random_segments[i].astype('str')) + "\n")

	return random_segments

def align_profiles_on_ft(data, feature, data_type, results_path):
	profiles = False
	data_before_init = []
	data_after_init = []
	all_profiles = []

	for protein in data:
		if feature in str(data[protein]):
			profiles = True
			before_init = []
			for record in data[protein]:
				if feature not in record[-1]:
					if data_type == 'RSCU':
						before_init.append(float(record[4]))
					elif data_type == 'codon_type':
						before_init.append(record[3] == 'TRUE')
				else:
					break
			if data_type == 'RSCU':		
				init_to_end	= [float(x[4]) for x in data[protein][len(before_init):]]
			elif data_type == 'codon_type':
				init_to_end	= [x[3] == 'TRUE' for x in data[protein][len(before_init):]]	

			data_before_init.append(before_init)
			data_after_init.append(init_to_end)

	if profiles == True:
		if feature != 'SIGNAL' and feature != 'INIT_MET':
			before_init_filled = np.array(list(itertools.zip_longest(*[x[::-1] for x in data_before_init], fillvalue=np.nan))).astype('float64')
			init_to_end_filled = np.array(list(itertools.zip_longest(*data_after_init, fillvalue=np.nan))).astype('float64')	
			average_init_to_end = np.nanmean(init_to_end_filled, axis=1)
			average_before_init = np.nanmean(before_init_filled, axis=1)[::-1]
			average_profile = np.concatenate((average_before_init, average_init_to_end), axis=None)
			number_of_profiles = before_init_filled.shape[1]
			point_0 = len(average_before_init)
			
			with open(f'{results_path}/{data_type}/{feature}_profiles.txt', 'w') as f:
				for i in range(0, number_of_profiles):
					prof = np.concatenate((before_init_filled[:, i][::-1], init_to_end_filled[:, i]), axis=None)[point_0-51:point_0+50]
					all_profiles.append(prof)
					f.write(args.organism + "\t" + (',').join(prof.astype('str')) + "\n")

		else:
			init_to_end_filled = np.array(list(itertools.zip_longest(*data_after_init, fillvalue=np.nan))).astype('float64')	
			average_init_to_end = np.nanmean(init_to_end_filled, axis=1)
			average_profile = average_init_to_end
			number_of_profiles = init_to_end_filled.shape[1]
			point_0 = 0

			with open(f'{results_path}/{data_type}/{feature}_profiles.txt', 'w') as f:
				for i in range(0, number_of_profiles):
					prof = init_to_end_filled[:, i][:101]
					all_profiles.append(prof)
					f.write(args.organism + "\t" + (',').join(prof.astype('str')) + "\n")

	elif profiles == False:
		average_profile, point_0 = [0, 0]

	return average_profile, point_0, profiles, all_profiles


def visualize(profile_data, output_path, feature, point_0, smoothing_window, data_type, graphs_path):	
		if feature != 'SIGNAL' and feature != 'INIT_MET':
				data = profile_data[point_0-51:point_0+50]
				plt.plot(data)
				plt.axvline(x=51, color='black', ls=':', lw=1)
				plt.xticks([10, 30, 51, 70, 90], ['-40', '-20', '0', '20', '40'])
		else:
				data = profile_data[:101]
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
			plt.savefig(f'{graphs_path}/{feature}_Average_RSCU.svg')
		elif data_type == 'codon_type':
			plt.ylabel("Relative Preferred Codon Usage", labelpad=15)
			plt.savefig(f'{graphs_path}/{feature}_RPCU.svg')

		plt.close()

def get_pair_wise_correlations(profiles):
	profiles_df = pd.DataFrame(profiles)
	profiles_corr = profiles_df.T.corr(method='pearson')
	unique_profiles_corr = profiles_corr.mask(np.tril(np.ones(profiles_corr.shape)).astype(bool))
	profiles_corr_dist = [value for value in unique_profiles_corr.stack().reset_index()[0]]

	return profiles_corr_dist


def perform_ttest(profiles, random_profiles, feature, data_type, results_path):
	if data_type == 'RSCU':
		profiles_corr_dist = get_pair_wise_correlations(profiles)
		random_profiles_corr_dist = get_pair_wise_correlations(random_profiles)

		pvalue = stats.ttest_ind(profiles_corr_dist, random_profiles_corr_dist, alternative='greater')[1]
		if pvalue <= 0.05:
			with open(f'{results_path}/RSCU_significant_profiles.txt', 'a+') as f:
				f.write(feature + "\t" f'significant with a p_value of: {pvalue}' + "\n")

	if data_type == 'codon_type':
		window_avg_profiles = []
		window_avg_rand_profiles = []
		smoothing_window = 3
		for profile in profiles:
			cumsum_vec = np.cumsum(np.insert(profile, 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			window_avg_profiles.append(smoothing_vec)

		for rand_profile in random_profiles:			
			cumsum_vec = np.cumsum(np.insert(rand_profile, 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			window_avg_rand_profiles.append(smoothing_vec)

		window_prof_corr_dist = get_pair_wise_correlations(window_avg_profiles)
		window_random_prof_corr_dist = get_pair_wise_correlations(window_avg_rand_profiles)

		pvalue = stats.ttest_ind(window_prof_corr_dist, window_random_prof_corr_dist, alternative='greater')[1]
		if pvalue <= 0.05:
			with open(f'{results_path}/RPCU_significant_profiles.txt', 'a+') as f:
				f.write(feature + "\t" f'significant with a p_value of: {pvalue}' + "\n")
					

######################################################################################
## ARGPARSE
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("codon2fts_file", type=str, help="codon2fts file to be processed")
parser.add_argument('output_path', nargs='?', type=str, help="desired path for writting output file", default=os.getcwd())
parser.add_argument('-f', '--feature_list', type=argparse.FileType('r'), help='file containing lists of features to be analyzed (one column, one feature per row)')
parser.add_argument('-d', '--data_type', type=str, help='type of data wanted to be analyzed: RSCU or codon_type')
parser.add_argument('-w', '--window', type=int, help='smoothing window size (must be an even number)')
parser.add_argument('-o', '--organism', type=str, help='name of the organism to be analyzed')
args = parser.parse_args()

######################################################################################
## MAIN
######################################################################################
organism_name = args.organism
results_path = os.path.join(args.output_path, f'results/{organism_name}')
graphs_path = os.path.join(results_path, 'average_profiles_correlations')
isExist_results = os.path.exists(results_path)
isExist_RSCU_results = os.path.exists(f'{results_path}/RSCU')
isExist_codon_type_results = os.path.exists(f'{results_path}/codon_type')

if isExist_results:
	if not isExist_RSCU_results and not isExist_codon_type_results:
		os.mkdir(f'{results_path}/RSCU')
		os.mkdir(f'{results_path}/codon_type')
	isExist_graph = os.path.exists(graphs_path)
	if not isExist_graph:
		os.mkdir(graphs_path)

codon2fts = load_tabular_as_dict(args.codon2fts_file)
features = args.feature_list.read().splitlines()
rand_profiles = create_random_profiles(codon2fts, args.data_type, results_path, organism_name)

print(f'List of features being analyzed for {organism_name}: {features}')

for feature in features:
		avg_profile, index_0, profiles, all_profiles = align_profiles_on_ft(codon2fts, feature, args.data_type, results_path)
		if profiles == True:
			perform_ttest(all_profiles, rand_profiles, feature, args.data_type, results_path)
			visualize(avg_profile, graphs_path, feature, index_0, args.window, args.data_type, graphs_path)

