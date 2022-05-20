#!/bin/bash

#this script was created to automate the process of optimization of the parameters m, M and n of Stacks programme

m_values=(2 3 4 5 6)
M_values=(1 2 3 4 5 6 7)
n_values=(1 2 3 4 5 6 7)


#argument 1 is the path where the demultiplexed files are stored (output of process_radtags)
$1
#argument 2 is a text file of the population map. It should contain two columns (individual_id and population)
$2


for m in ${m_values[@]}; do
	#create a folder for each parameter
	mkdir m$m
	#stacks de novo assembly for a certain set of parameters
	denovo_map.pl -T 8 -M 2 -m "$m" -n 1 -o "./m$m" --samples $1 --popmap $2 --paired -X "populations: -r 0.80" -X "populations: -p 4"
	#extract a table from stacks output that described the number of loci per number of snps
	sed -n '/BEGIN snps_per_loc_postfilters/,/END snps_per_loc_postfilters/ p' ./m$m/populations.log.distribs | grep -E '^[0-9]' > ./m$m/snps_per_loc.log
	#use vcftools to calculate the depth of coverage for every snp
	vcftools --vcf ./m${m}/populations.snps.vcf --out ./m$m/m${m} --depth
	#use an R script to calculate the mean depth of coverage
	./coverageplot.R ./m${m}/m${m}.idepth ./m$m/dep.csv
	#create the file for each parameter, for coverage
	echo -e  "m\t"m$m"\tcoverage" > ./m$m/temp_cov_m$m.tsv
	paste ./m$m/temp_cov_m$m.tsv ./m$m/dep.csv > ./m$m/m${m}_cov.tsv
	cat ./m$m/m${m}_cov.tsv >> m_coverage_final.tsv
done

#this loop does the same as the previous one but varying M parameter
for M in ${M_values[@]}; do
	denovo_map.pl -T 8 -M "$M" -m 3 -n 1 -o "./M$M" --samples $1 --popmap $2 --paired -X "populations: -r 0.80" -X "populations: -p 4"
	sed -n '/BEGIN snps_per_loc_postfilters/,/END snps_per_loc_postfilters/ p' ./M$M/populations.log.distribs | grep -E '^[0-9]' > ./M$M/snps_per_loc.log
	vcftools --vcf ./M${M}/populations.snps.vcf --out ./M$M/M${M} --depth
	./coverageplot.R ./M${M}/M${M}.idepth ./M$M/dep.csv
	 echo -e  "M\t"M$M"\tcoverage" > ./M$M/temp_cov_M$M.tsv
	paste ./M$M/temp_cov_M$M.tsv ./M$M/dep.csv > ./M$M/M${M}_cov.tsv
	cat ./M$M/M${M}_cov.tsv >> M_coverage_final.tsv
done

#this loop does the same as the previous one but varying the n parameter
for n in ${n_values[@]}; do
	denovo_map.pl -T 8 -M 2 -m 3 -n "$n" -o "./n$n" --samples $1 --popmap $2 --paired -X "populations: -r 0.80" -X "populations: -p 4"
	sed -n '/BEGIN snps_per_loc_postfilters/,/END snps_per_loc_postfilters/ p' ./n$n/populations.log.distribs | grep -E '^[0-9]' > ./n$n/snps_per_loc.log
	vcftools --vcf ./n${n}/populations.snps.vcf --out ./n$n/n${n} --depth
	./coverageplot.R ./n${n}/n${n}.idepth ./n$n/dep.csv
	 echo -e  "n\t"n$n"\tcoverage" > ./n$n/temp_cov_n$n.tsv
	paste ./n$n/temp_cov_n$n.tsv ./n$n/dep.csv > ./n$n/n${n}_cov.tsv
	cat ./n$n/n${n}_cov.tsv >> n_coverage_final.tsv
done


for m in ${m_values[@]}; do
	#use an R script to extract the metrics values for each parameter based on the output file of stacks
	./extract_metrics_files.R ./m$m/snps_per_loc.log ./m$m/nr_ass_loci.txt ./m$m/nr_poly_loci.txt ./m$m/nr_total_snps.txt
	#create the file for each parameter, for the assembled loci
	echo -e  "m\t"m$m"\tassembled_loci" > ./m$m/temp_ass_m$m.tsv
	paste ./m$m/temp_ass_m$m.tsv ./m$m/nr_ass_loci.txt >> m_assloci_final.tsv
	#create the file for each parameter, for the polymorphic loci
	echo -e  "m\t"m$m"\tpoly_loci" > ./m$m/temp_poly_m$m.tsv
	paste ./m$m/temp_poly_m$m.tsv ./m$m/nr_poly_loci.txt >> m_polyloci_final.tsv
	#create the file for each parameter, for the total number of snps
	echo -e  "m\t"m$m"\ttotal_snps" > ./m$m/temp_snps_m$m.tsv
	paste ./m$m/temp_snps_m$m.tsv ./m$m/nr_total_snps.txt >> m_totalsnps_final.tsv
	#remove the intermidiate, temporary files
	rm ./m$m/*temp*
done

#this loop does the same as the previous one but varying the M parameter
for M in ${M_values[@]}; do
	./extract_metrics_files.R ./M$M/snps_per_loc.log ./M$M/nr_ass_loci.txt ./M$M/nr_poly_loci.txt ./M$M/nr_total_snps.txt
	echo -e  "M\t"M$M"\tassembled_loci" > ./M$M/temp_ass_M$M.tsv
	paste ./M$M/temp_ass_M$M.tsv ./M$M/nr_ass_loci.txt >> M_assloci_final.tsv
	echo -e  "M\t"M$M"\tpoly_loci" > ./M$M/temp_poly_M$M.tsv
	paste ./M$M/temp_poly_M$M.tsv ./M$M/nr_poly_loci.txt >> M_polyloci_final.tsv
	echo -e  "M\t"M$M"\ttotal_snps" > ./M$M/temp_snps_M$M.tsv
	paste ./M$M/temp_snps_M$M.tsv ./M$M/nr_total_snps.txt >> M_totalsnps_final.tsv
	rm ./M$M/*temp*
done

#this loop does the same as the previous one but varying the n parameter
for n in ${n_values[@]}; do
	./extract_metrics_files.R ./n$n/snps_per_loc.log ./n$n/nr_ass_loci.txt ./n$n/nr_poly_loci.txt ./n$n/nr_total_snps.txt
	echo -e  "n\t"n$n"\tassembled_loci" > ./n$n/temp_ass_n$n.tsv
	paste ./n$n/temp_ass_n$n.tsv ./n$n/nr_ass_loci.txt >> n_assloci_final.tsv
	echo -e  "n\t"n$n"\tpoly_loci" > ./n$n/temp_poly_n$n.tsv
	paste ./n$n/temp_poly_n$n.tsv ./n$n/nr_poly_loci.txt >> n_polyloci_final.tsv
	echo -e  "n\t"n$n"\ttotal_snps" > ./n$n/temp_snps_n$n.tsv
	paste ./n$n/temp_snps_n$n.tsv ./n$n/nr_total_snps.txt >> n_totalsnps_final.tsv
	rm ./n$n/*temp*
done

#create the final files to plot the results. the is one final file for each metric. the final files have four columns. one with the parameter (m,M or n), one with the metric value (e.g. m3,M2,n4),
#one with the metric (e.g. coverage) and one with the metric_value(e.g. 27.12)
cat m_coverage_final.tsv M_coverage_final.tsv n_coverage_final.tsv > coverage_opt.tsv #coverage
cat m_assloci_final.tsv M_assloci_final.tsv n_assloci_final.tsv > assloci_opt.tsv #assembled loci
cat m_polyloci_final.tsv M_polyloci_final.tsv n_polyloci_final.tsv > polyloci_opt.tsv #polymorphic loci
cat m_totalsnps_final.tsv M_totalsnps_final.tsv n_totalsnps_final.tsv > totalsnps_opt.tsv #total number of snps
#remove the temporary files
rm *final*
