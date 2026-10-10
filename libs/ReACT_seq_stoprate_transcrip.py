#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
# ReACT_seq_stoprate_transcrip.py
# ---------------------------------------------------------------------------
# Purpose: compute per-site RT-stop rates from samtools mpileup output
#          (plus strand: 5'-start events '^'; minus strand: 3'-end events
#          '$'), restricted to gene/transcript intervals; feeds the
#          downstream high-stop-site calling step.
#
# Background:
#   Core quantification script of the RT-stop (ReACT-seq) analysis.
#
# Changes vs. the historical version kept in the monorepo:
#   1. gene-interval index switched from dict (keyed by gene name) to a
#      coordinate-sorted list with a sliding window, avoiding key collisions
#      when intervals overlap / multiple genes share coordinates;
#   2. a three-point sliding window (fir/sed/thr) confirms that a stop site
#      falls inside a gene interval (front/mid/rear check instead of a
#      single-point test), reducing boundary misassignment;
#   3. the old coverage < 10 filter was removed; the output layout
#      (7 columns) is unchanged, and stop_pos/stop_base are taken from the
#      preceding position of the three-point window.
#
# Usage:
#   python3 ReACT_seq_stoprate_transcrip.py <pileup> <plus.bed> <minus.bed> \
#       <plus_out.txt> <minus_out.txt>
#
# Inputs:
#   pileup    : samtools mpileup output (chr, pos, ref, cov, base...)
#   plus.bed  : plus-strand gene/transcript intervals (chr, start, end, gene)
#   minus.bed : minus-strand gene/transcript intervals (chr, start, end, gene)
# Output (tab separated):
#   chr_name stop_pos stop_base stop_cov total_cov stop_rate gene
# ---------------------------------------------------------------------------
import re
import sys

pileup_file = sys.argv[1]
plus_file = sys.argv[2]
minus_file =  sys.argv[3]
plus_out = sys.argv[4]
minus_out =  sys.argv[5]


plus_text = open(plus_file,'rt')
plus_site_dict = {}
for line in plus_text:
	linelist = re.split("\t",line.rstrip())[0:4]
	if not linelist[0] in plus_site_dict:
		plus_site_dict[linelist[0]] = []
	plus_site_dict[linelist[0]].append([int(linelist[1]),int(linelist[2]),linelist[3]])
plus_text.close()

minus_text = open(minus_file,'rt')
minus_site_dict = {}
for line in minus_text:
	linelist = re.split("\t",line.rstrip())[0:4]
	if not linelist[0] in minus_site_dict:
		minus_site_dict[linelist[0]] = []
	minus_site_dict[linelist[0]].append([int(linelist[1]),int(linelist[2]),linelist[3]])
minus_text.close()

pileup_text = open(pileup_file,'rt')
plus_out_text = open(plus_out,'wt')
minus_out_text = open(minus_out,'wt')
plus_out_text.write("chr_name\tstop_pos\tstop_base\tstop_cov\ttotal_cov\tstop_rate\tgene\n")
minus_out_text.write("chr_name\tstop_pos\tstop_base\tstop_cov\ttotal_cov\tstop_rate\tgene\n")
fir_linelist = []
sed_linelist = []
thr_linelist = []
for line in pileup_text:  
	now_linelist =  re.split("\t",line)[0:5]
	if not (now_linelist[0] in plus_site_dict or now_linelist[0] in minus_site_dict):
		continue
	#get three sites info
	if not (fir_linelist):
		fir_linelist = now_linelist
	elif not (sed_linelist):
		sed_linelist = now_linelist
	elif not (thr_linelist):
		thr_linelist = now_linelist
	else:
		fir_linelist = sed_linelist
		sed_linelist = thr_linelist
		thr_linelist = now_linelist
	if not (thr_linelist):
		continue
	#confirm stop site in gene
	gene_list_plus = []
	gene_list_minus = []
	for gene in plus_site_dict[sed_linelist[0]]:
		if int(fir_linelist[1]) <= gene[0]:
			break
		elif int(fir_linelist[1]) > gene[0] and int(thr_linelist[1]) < gene[1]:
			gene_list_plus.append(gene[2])
		elif int(thr_linelist[1]) >= gene[1]:
			plus_site_dict[sed_linelist[0]].remove(gene)
		else:
			print(fir_linelist[0],fir_linelist[1],"error",sep = "\t")
	if gene_list_plus:
		start_count = len(re.findall(r'\^',sed_linelist[4]))
		start_rate = round(start_count/int(sed_linelist[3])*100,2)
		start_trans = ','.join(gene_list_plus)
		print(fir_linelist[0],fir_linelist[1],fir_linelist[2],start_count,sed_linelist[3],start_rate,start_trans,end = "\n",sep = "\t",file = plus_out_text)
	for gene in minus_site_dict[sed_linelist[0]]:
		if int(fir_linelist[1]) < gene[0]:
			break
		elif int(fir_linelist[1]) >= gene[0] and int(thr_linelist[1]) < gene[1]-1:
			gene_list_minus.append(gene[2])
		elif int(thr_linelist[1]) >= gene[1]-1:
			minus_site_dict[sed_linelist[0]].remove(gene)
		else:
			print(thr_linelist[0],thr_linelist[1],"error",sep = "\t")
	if gene_list_minus:
		stop_count = len(re.findall(r'\$',sed_linelist[4]))
		stop_rate = round(stop_count/int(sed_linelist[3])*100,2)
		stop_trans = ','.join(gene_list_minus)
		print(thr_linelist[0],thr_linelist[1],thr_linelist[2],stop_count,sed_linelist[3],stop_rate,stop_trans,end = "\n",sep = "\t",file = minus_out_text)
pileup_text.close()
plus_out_text.close()
minus_out_text.close()
