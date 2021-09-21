#!/bin/bash

#this script was created to automate the process of optimization of the parameters m, M and n of Stacks programme
#the conda environment that hosts Stacks needs to be activated

#argument 1 is the path where the demultiplexed files are stored (output of process_radtags)
$1
#argument 2 is a text file of the population map. It should contain two columns (individual_id and population)
$2

#parameter values that are intended to be tested
#each parameter value is tested separately, as the other parameter stay on default (M=2 m=3 n=1)

m_values=( 2 3 4 5 6 )
M_values=( 1 2 3 4 5 6 7 8 )
n_values=( 2 3 4 5 6 7 )

#variables for the parameter values needed to create the final file
VARmv=""
VARMv=""
VARnv=""

#files needed to extract the metrics from the Stacks's output files
echo "var_ass=" > var_ass.tsv
echo "var_poly=" > var_poly.tsv
echo "var_snps" > var_snps.tsv

#the comments inside the following loop are equivalent for the next two loops, but referent to 

for m in "${m_values[@]}"; do
	
	#create folder for each parameter value
	mkdir m$m
	VARmv="m$m"
	#run denovo assembly of Stacks for each parameter value  
	denovo_map.pl -T 8 -M 2 -m "$m" -n 1 -o "./m$m" --samples $1 --popmap $2 --paired -X "populations: -r 0.80"
	
	#extract the number of assembled loci resulted from the assembly and create a file with that information
	cat ./m${m}/populations.haplotypes.tsv | grep -v ^"#" | wc -l > ./m${m}/m${m}_nr_assembled_loci.tsv
	paste -d'\0' ./var_ass.tsv ./m${m}/m${m}_nr_assembled_loci.tsv > ./m${m}/m${m}_assembled_loci.tsv
	rm ./m${m}/m${m}_nr_assembled_loci.tsv
	. ./m${m}/m${m}_assembled_loci.tsv
	echo -e  "m\t"$VARmv"\tassembled_loci\t"$var_ass >> m_assloci_final.tsv

	#extract the number of polymorphic loci resulted from the assembly and create a file with that information
	cat ./m${m}/populations.hapstats.tsv | grep -v ^"#" | wc -l > ./m${m}/m${m}_nr_polymorphic_loci.tsv
	paste -d'\0' ./var_poly.tsv ./m${m}/m${m}_nr_polymorphic_loci.tsv > ./m${m}/m${m}_polymorphic_loci.tsv
	rm ./m${m}/m${m}_nr_polymorphic_loci.tsv
	. ./m${m}/m${m}_polymorphic_loci.tsv
	echo -e  "m\t"$VARmv"\tpoly_loci\t"$var_poly >> m_polyloci_final.tsv
	
	#extract the number of SNPs resulted from the assembly and create a file with that information
	cat ./m${m}/populations.sumstats.tsv | grep -v ^"#" | wc -l > ./m${m}/m${m}_nr_total_snps.tsv
	paste -d'\0' ./var_snps.tsv ./m${m}/m${m}_nr_total_snps.tsv > ./m${m}/m${m}_total_snps.tsv
	rm ./m${m}/m${m}_nr_total_snps.tsv
	. ./m${m}/m${m}_total_snps.tsv
	echo -e  "m\t"$VARmv"\ttotal_snps\t"$var_snps >> m_snps_final.tsv
done

for M in "${M_values[@]}"; do
	mkdir M$M
	VARmv="M$M"
	denovo_map.pl -T 8 -M "$M" -m 3 -n 1 -o "./M$M" --samples $1 --popmap $2 --paired -X "populations: -r 0.80"
	
	cat ./M${M}/populations.haplotypes.tsv | grep -v ^"#" | wc -l > ./M${M}/M${M}_nr_assembled_loci.tsv
	paste -d'\0' ./var_ass.tsv ./M${M}/M${M}_nr_assembled_loci.tsv > ./M${M}/M${M}_assembled_loci.tsv
	rm ./M${M}/M${M}_nr_assembled_loci.tsv
	. ./M${M}/M${M}_assembled_loci.tsv
	echo -e  "M\t"$VARMv"\tassembled_loci\t"$var_ass >> M_assloci_final.tsv

	cat ./M${M}/populations.hapstats.tsv | grep -v ^"#" | wc -l > ./M${M}/M${M}_nr_polymorphic_loci.tsv
	paste -d'\0' ./var_poly.tsv ./M${M}/M${M}_nr_polymorphic_loci.tsv > ./M${M}/M${M}_polymorphic_loci.tsv
	rm ./M${M}/M${M}_nr_polymorphic_loci.tsv
	. ./M${M}/M${M}_polymorphic_loci.tsv
	echo -e  "M\t"$VARMv"\tpoly_loci\t"$var_poly >> M_polyloci_final.tsv
	
	cat ./M${M}/populations.sumstats.tsv | grep -v ^"#" | wc -l > ./M${M}/M${M}_nr_total_snps.tsv
	paste -d'\0' ./var_snps.tsv ./M${M}/M${M}_nr_total_snps.tsv > ./M${M}/M${M}_total_snps.tsv
	rm ./M${M}/M${M}_nr_total_snps.tsv
	. ./M${M}/M${M}_total_snps.tsv
	echo -e  "M\t"$VARMv"\ttotal_snps\t"$var_snps >> M_snps_final.tsv 
	
done 

for n in "${n_values[@]}"; do
	mkdir n$n
	VARnv="n$n"
	denovo_map.pl -T 8 -M 2 -m 3 -n "$n" -o "./n$n" --samples $1 --popmap $2 --paired -X "populations: -r 0.80" 

	cat ./n${n}/populations.haplotypes.tsv | grep -v ^"#" | wc -l > ./n${n}/n${n}_nr_assembled_loci.tsv
	paste -d'\0' ./var_ass.tsv ./M${M}/M${M}_nr_assembled_loci.tsv > ./M${M}/M${M}_assembled_loci.tsv 
	rm ./M${M}/M${M}_nr_assembled_loci.tsv
	. ./M${M}/M${M}_assembled_loci.tsv
	echo -e  "M\t"$VARMv"\tassembled_loci\t"$var_ass >> M_assloci_final.tsv

	cat ./n${n}/populations.hapstats.tsv | grep -v ^"#" | wc -l > ./n${n}/n${n}_nr_polymorphic_loci.tsv
	paste -d'\0' ./var_poly.tsv ./M${M}/M${M}_nr_polymorphic_loci.tsv > ./M${M}/M${M}_polymorphic_loci.tsv
	rm ./M${M}/M${M}_nr_polymorphic_loci.tsv
	. ./M${M}/M${M}_polymorphic_loci.tsv 
	echo -e  "M\t"$VARMv"\tpoly_loci\t"$var_poly >> M_polyloci_final.tsv

	cat ./n${n}/populations.sumstats.tsv | grep -v ^"#" | wc -l > ./n${n}/n${n}_nr_total_snps.tsv
	paste -d'\0' ./var_snps.tsv ./n${n}/n${n}_nr_total_snps.tsv > ./n${n}/n${n_total_snps.tsv
	rm ./n${n}/n${n}_nr_total_snps.tsv
	. ./n${n}/n${n}_total_snps.tsv
	echo -e  "n\t"$VARnv"\ttotal_snps\t"$var_snps >> n_snps_final.tsv
done

#create final files for input for R
#each file contains the results from only one metric and for all parameters

cat m_assloci_final.tsv M_assloci_final.tsv n_assloci_final.tsv > ass_loci_final.tsv
cat m_poly_final.tsv M_poly_final.tsv n_poly_loci.tsv > poly_loci_final.tsv
cat m_snps_final.tsv M_snps_final.tsv n_snps_final.tsv > snps_final.tsv


#run R script that builds bar plots of the results of parameter optimization

for i in ass_loci poly_loci snps; do
	./param_optimization_plots.R ${i}_final.tsv optimization_${i}_barplot.png
done
