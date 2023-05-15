#! /bin/bash

while getopts 'cv' OPTION; do
	case "$OPTION" in
		c)
			python3 get_annot_sequences.py $2 $3 $4 -v -o $5
			;;
		v)
			proteins=$2
			for i in ${proteins//,/ }
			do
        			sed -n '/AC   '"$i"'/{x;p;d;}; x' $4 > proteins_records.txt
        			sed -n '/AC   '"$i"'/,${p;/^\/\//q}' $4 >> proteins_records.txt
        			python3 visualize_proteins.py $3 proteins_records.txt -p $i -w 8
				rm proteins_records.txt
			done
			;;
	esac
done	
