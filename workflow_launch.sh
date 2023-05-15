! /bin/env bash

ls data > organism_list
mkdir -p results
while read organism; do
	./seq2features.sh -c data/$organism/T* data/$organism/uni* data/$organism/nuclear_codon_statistics.tsv $organism
	python3 analyze_features.py results/$organism/annot_prot_by_codons.txt -f feature_list -d RSCU -w 2 -o $organism
	python3 analyze_features.py results/$organism/annot_prot_by_codons.txt -f feature_list -d codon_type -w 2 -o $organism
done <organism_list

mkdir -p results/RSCU_all_organism_profiles
mkdir -p results/codon_type_all_organism_profiles
mkdir -p results/RSCU_all_organism_profiles/graphs
mkdir -p results/codon_type_all_organism_profiles/graphs

while read feature; do
	cat results/*/RSCU/$feature* > results/RSCU_all_organism_profiles/"$feature"_all_organism_profiles
	cat results/*/codon_type/$feature* > results/codon_type_all_organism_profiles/"$feature"_all_organism_profiles
	cat results/*/RSCU/random_profiles.txt > results/RSCU_all_organism_profiles/all_random_profiles	
	cat results/*/codon_type/random_profiles.txt > results/codon_type_all_organism_profiles/all_random_profiles
	python3 perform_all_organism_ttest.py results/RSCU_all_organism_profiles/"$feature"_all_organism_profiles results/RSCU_all_organism_profiles/all_random_profiles -d RSCU -f $feature -w 2	
	python3 perform_all_organism_ttest.py results/codon_type_all_organism_profiles/"$feature"_all_organism_profiles results/codon_type_all_organism_profiles/all_random_profiles -d codon_type -f $feature -w 2
	python3 all_organism_graphs.py results/RSCU_all_organism_profiles/"$feature"_all_organism_profiles -d RSCU -f "$feature" -w 2	
	python3 all_organism_graphs.py results/codon_type_all_organism_profiles/"$feature"_all_organism_profiles -d codon_type -f "$feature" -w 2	
done <feature_list

python3 compare_data_types.py results/coli/annot_prot_by_codons.txt  results/coli/RSCU/HELIX_profiles.txt results/coli/codon_type/HELIX_profiles.txt data/coli/\[tRNA\]_data.txt -f 'HELIX' -w 2

python3 compare_data_types.py results/coli/annot_prot_by_codons.txt  results/coli/RSCU/HELIX_profiles.txt results/coli/codon_type/HELIX_profiles.txt data/coli/\[tRNA\]_data.txt -f 'TURN' -w 2
