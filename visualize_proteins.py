#! /usr/bin/env python3

import os
import argparse
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO
import numpy as np

##########################################################################################
## METHODS
##########################################################################################
class MyCustomTranslator(BiopythonTranslator):

	graphic_record_parameters = {
		"labels_spacing": 4,
	}

	def compute_feature_color(self, feature):
		if feature.type == "CHAIN":
			return "green"
		elif feature.type == "DOMAIN":
		    return "red"
		elif feature.type == "TOPO_DOM":
		    return "gold"
		elif feature.type == "TRANSMEM":
		   	return "gray"	    
		else:
		    return "blue"

	def compute_feature_label(self, feature):
		if feature.type == "BINDING":
			if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
				return f"BINDING: {BiopythonTranslator.compute_feature_label(self, feature)}"
			else:
				return feature.type	
		elif feature.type == "ACT_SITE":
			if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
				return f"ACT_SITE: {BiopythonTranslator.compute_feature_label(self, feature)}"
			else:
				return feature.type
		#elif feature.type == "MUTAGEN":
		#	if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
		#		return f"MUTAGEN: {BiopythonTranslator.compute_feature_label(self, feature)}"
		#	else:
		#		return feature.type
		#elif feature.type == "TOPO_DOM":
		#	if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
		#		return f"{BiopythonTranslator.compute_feature_label(self, feature)}"
		#	else:
		#		return feature.type									
		#elif feature.type == "SITE":
		#	if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
		#		return f"SITE: {BiopythonTranslator.compute_feature_label(self, feature)}"
		#	else:
		#		return feature.type	
		#elif feature.type == "REGION":
		#	if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
		#		return f"REGION: {BiopythonTranslator.compute_feature_label(self, feature)}"
		#	else:
		#		return feature.type	
		#elif feature.type == "MOTIF":
		#	if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
		#		return f"MOTIF: {BiopythonTranslator.compute_feature_label(self, feature)}"
		#	else:
		#		return feature.type
		elif feature.type == "CHAIN":
			if BiopythonTranslator.compute_feature_label(self, feature) != feature.type:
				return f"{BiopythonTranslator.compute_feature_label(self, feature)}"
			else:
				return feature.type																		
		else:
			return feature.type
		#return feature.type

	def compute_filtered_features(self, features):
		return[
			feature for feature in features
			if (feature.type != "MUTAGEN" and feature.type != "VARIANT" and feature.type != "BINDING" and feature.type != "LIPID" and feature.type != "CONFLICT")
		]
		
def load_tabular_file(file):
	fields = []
	with open(file, 'rt') as f:
		records = f.read().splitlines()
		for record in records:
			fields.append(record.split("\t"))
	return fields

def get_protein_list(codon2fts, proteins):
	prot_list = []
	for protein in proteins:
		prot_list.append([x for x in codon2fts if x[0] == protein])

	return prot_list

def visualize(protein_data, output_path, smoothing_window):
	for protein in protein_data:
		protein_name = protein[0][0]
		seq = [x[2] for x in protein if x[4] == 'TRUE']
		pos = [float(x[3]) for x in protein]
		pos_label = [float(x[3])-1 for x in protein if x[4] == 'TRUE']
		RSCU = [float(x[5]) for x in protein]
		prot_features = [x[-1].split(',') if "," in x[-1] else x[-1] for x in protein]

		if smoothing_window:
			cumsum_vec = np.cumsum(np.insert(RSCU, 0, 0)) 
			smoothing_vec = (cumsum_vec[smoothing_window:] - cumsum_vec[:-smoothing_window]) / smoothing_window
			smoothing_positions = [x for x in range(smoothing_window//2, len(smoothing_vec)+smoothing_window//2)]

		
		fig, (ax1, ax2) = plt.subplots(
			2, 1, figsize=(16, 7), sharex=True, gridspec_kw={'height_ratios': [2, 2]})

		record = SeqIO.read(args.protein_records, "swiss")


		graphic_record = MyCustomTranslator().translate_record(record)
		graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=1)
		
		plt.sca(ax2)
		ax2.plot(smoothing_positions, smoothing_vec)
		plt.xticks(pos_label, seq, fontsize = 5)
		ax2.tick_params(axis='x', colors='red', width=3, direction= 'in', labelbottom=False)
		plt.xlabel("AA sequence", labelpad=10)
		plt.ylabel("RSCU Values", labelpad=15)

		fig.suptitle(f'{protein_name} Features and RSCU values')
		fig.savefig(f'{output_path}/{protein_name}_fts_RSCU.svg')

######################################################################################
## ARGPARSE
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("codon2fts_file", type=str, help="codon2fts file to be processed")
parser.add_argument("protein_records", type=str, help="parsed protein records file")
parser.add_argument('output_path', nargs='?', type=str, help="desired path for writting output file", default=os.getcwd())
parser.add_argument('-p', '--protein', type=str, help='protein AC to be specificially returned')
#parser.add_argument('-d', '--data_type', type=str, help='type of profile data to display')
parser.add_argument('-w', '--window', type=int, help='smoothing window size (must be an even number)')
args = parser.parse_args()

######################################################################################
## MAIN
######################################################################################

codon2fts = load_tabular_file(args.codon2fts_file)
proteins = args.protein.split(',')
protein_data = get_protein_list(codon2fts, proteins)
graphs_path = os.path.join(args.output_path, "Protein examples graphs")
isExist = os.path.exists(graphs_path)
if not isExist:
	os.mkdir(graphs_path)

visualize(protein_data, graphs_path, args.window)

for protein in protein_data:
	print(f'-------------------------------------{protein[0][0]}-----------------------------------------')	
	for codon in protein:
		line = '\t'.join(codon)
		print(line)
print('\n')	
print(f"################ SVG graph saved at {graphs_path}/{protein[0][0]}_fts_RSCU.svg ################")