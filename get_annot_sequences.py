#! /usr/bin/env python3

import os
import argparse
import gzip
import re
import pickle

##########################################################################################
## METHODS
##########################################################################################

def load_gz_file(file):
	with gzip.open(file, 'rt') as f:
		records = f.read().splitlines()
	return records

def load_file(file):
	with open(file, 'rt') as f:
		records = f.read().splitlines()
	return records

def load_tabular_file(file):
	fields = []
	with open(file, 'rt') as f:
		records = f.read().splitlines()
		for record in records:
			fields.append(record.split("\t"))
	return fields

def parse_kegg_file(kegg_file):
	uniprot_nt_aa_dict = dict({})
	uniprot_code = []
	used_uniprot_codes = []
	in_aaseq_line = False
	in_ntseq_line = False

	for line in kegg_file:
			fields = line.split(maxsplit = 1)
			if fields[0] == 'UniProt:':
				uniprot_code = fields[1]
				is_repeated = used_uniprot_codes.count(uniprot_code)
				used_uniprot_codes.append(uniprot_code)
				if is_repeated > 0 and len(uniprot_code.split()) > 1:
					uniprot_code = uniprot_code.split()[0]
				uniprot_nt_aa_dict[uniprot_code] = [[],[]]
			elif fields[0] == 'AASEQ':
				in_aaseq_line = True
				next
			elif fields[0] == 'NTSEQ':
				in_ntseq_line = True
				in_aaseq_line = False
				next
			elif fields[0] == '///':
				if uniprot_code:
					uniprot_nt_aa_dict[uniprot_code][0] = ''.join(uniprot_nt_aa_dict[uniprot_code][0]).upper()
					uniprot_nt_aa_dict[uniprot_code][1] = ''.join(uniprot_nt_aa_dict[uniprot_code][1])
				in_ntseq_line = False
				uniprot_code = False
			elif uniprot_code:
				if in_aaseq_line:
					uniprot_nt_aa_dict[uniprot_code][1].append(line.strip())
				elif in_ntseq_line:
					uniprot_nt_aa_dict[uniprot_code][0].append(line.strip())

	return uniprot_nt_aa_dict

def parse_uniprot_file(uniprot_file):
	id_aa_ft_dict = dict({})
	uniprot_code_reps = []
	remove_no_ft_entries = []
	SQ_lines = False
	first_AC_line_checkpoint = 0

	for line in uniprot_file:
		fields = line.split(maxsplit = 1)
		if fields[0] == 'AC':
			first_AC_line_checkpoint += 1 
			if first_AC_line_checkpoint == 1:
				uniprot_code = fields[1].split(';')[0]
				id_aa_ft_dict[uniprot_code] = [[],[]]
		elif re.match('^FT\s{3}\w+\s+\d+.?.?\d*$', line):
			if fields[1].split()[0] != 'CHAIN' and fields[1].split()[0] != 'CONFLICT':
				id_aa_ft_dict[uniprot_code][1].append(fields[1])
		elif fields[0] == 'SQ':
			SQ_lines = True
			next
		elif line == '//':
			SQ_lines = False
			first_AC_line_checkpoint = 0
		elif SQ_lines:
			id_aa_ft_dict[uniprot_code][0].append(line.strip().replace(' ',''))

	for uniprot_code in id_aa_ft_dict:
		if len(id_aa_ft_dict[uniprot_code][1]) == 0:
			remove_no_ft_entries.append(uniprot_code)		
		else:
			id_aa_ft_dict[uniprot_code][0] = ''.join(id_aa_ft_dict[uniprot_code][0])

	for entry in remove_no_ft_entries:
		del id_aa_ft_dict[entry]			

	if args.verbose:
		print('{} entries have been removed due to not having any FT of interest'.format(len(remove_no_ft_entries)))

	return id_aa_ft_dict

def check_kegg2uniprot_aasq(id_nt_aa_dict, id_aa_ft_dict):
	equals = []
	not_equals = []
	for i in id_nt_aa_dict:
		if i in id_aa_ft_dict.keys():
			if id_aa_ft_dict[i][0] == id_nt_aa_dict[i][1]:
					equals.append(i)
			else:
					not_equals.append(i)
	if args.verbose:				
			print('{} do not have the same sequence in both databases, so they will be excluded from the analysis'.format(str(not_equals)))					

	return not_equals

def parse_codon_usage_file(codon_usage_file):
	codon_stats = dict({})
	for record in codon_usage_file:
		if len(record) > 1 and record[0] != 'CODON':
			codon_stats[record[0]] = [record[2], record[6]]
	return codon_stats

def split_in_codons(nt_sequence):
	used_codons = [nt_sequence[i:i+3] for i in range(0, len(nt_sequence), 3)]

	return used_codons

def map_ft_to_ntseq(id_aa_ft_dict, id_nt_aa_dict, not_equals, codon_stats):
	ft_to_nt_dict = dict({})
	codon_type = None
	for i in id_nt_aa_dict:
		if i not in not_equals and i in id_aa_ft_dict.keys():
			codons = split_in_codons(id_nt_aa_dict[i][0])
			if len(codons[-1]) != 3:
				continue
			ft_to_nt_dict[i] = []
			for idx, codon in enumerate(codons, start = 1):
				codon_features = []
				for feature in id_aa_ft_dict[i][1]:
					ft, position = feature.split()
					if ft == 'DISULFID' or ft == 'CROSSLNK':
						if len(position.split('..')) == 1:
							if idx == int(position):
								codon_features.append(ft)
						elif len(position.split('..')) == 2:
							start, end = position.split('..')
							if idx == int(start) or idx == int(end): 
								codon_features.append(ft)
					else:						
						if len(position.split('..')) == 1:
							if idx == int(position):
								codon_features.append(ft)
						elif len(position.split('..')) == 2:
							start, end = position.split('..')
							if idx in range(int(start), int(end)+1):
								codon_features.append(ft)

				ft_to_nt_dict[i].append([codon, str(idx), codon_stats[codon], codon_features])

	return ft_to_nt_dict

######################################################################################
## ARGPARSE
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("kegg_file", type=str, help="kegg_db file to be processed")
parser.add_argument("uniprot_file", type=str, help="uniprot_db file to be processed")
parser.add_argument("codon_usage_file", type=str, help="codon_usage_db file to be processed")
parser.add_argument('output_path', nargs='?', type=str, help="desired path for writting output file", default=os.getcwd())
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-o', '--organism', type=str, help='name of the organism to be analyzed')
args = parser.parse_args()

######################################################################################
## MAIN
######################################################################################

kegg_file = load_gz_file(args.kegg_file)
uniprot_file = load_file(args.uniprot_file)
codon_stats_file = load_tabular_file(args.codon_usage_file)
kegg_dict = parse_kegg_file(kegg_file)
uniprot_dict = parse_uniprot_file(uniprot_file)
codon_stats_dict = parse_codon_usage_file(codon_stats_file)
not_equal_aaseq = check_kegg2uniprot_aasq(kegg_dict, uniprot_dict)
codon2fts = map_ft_to_ntseq(uniprot_dict, kegg_dict, not_equal_aaseq, codon_stats_dict)
organism_name = args.organism
results_path = os.path.join(args.output_path, f'results/{organism_name}')
isExist_results = os.path.exists(results_path)

if not isExist_results:
	os.mkdir(results_path)

with open(f'{results_path}/annot_prot_by_codons.txt', 'w') as f:	
	for prot in codon2fts:
		for codon in codon2fts[prot][:-1]:
			f.write(prot + "\t" + codon[0] + "\t" + uniprot_dict[prot][0][(int(codon[1])-1)] + "\t" + str(codon[1]) + "\t" + codon[2][0] + "\t" + codon[2][1] + "\t" + ','.join(codon[3]) + "\n")