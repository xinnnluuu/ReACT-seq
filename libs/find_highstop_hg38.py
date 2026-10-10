#!/usr/bin/env python3
# find_highstop_hg38.py — call high-stop sites for the human_trans preset:
# compares per-site RT-stop rates of treated vs control samples (plus/minus
# stoprate tables from ReACT_seq_stoprate_transcrip.py) and flags sites
# enriched in stops within candidate peak regions (hg38 track).
# Usage: python3 find_highstop_hg38.py <treat_plus.txt> <treat_minus.txt> \
#           <output.highstop.txt> <peaks.bed>
# Output columns: chr_name site_pos site_ref strand_type t_stopcount
#                 t_stoprate stoprate_dif_minus1 stoprate_dif_plus1 gene_id
import sys
import re

treat_plus_file = sys.argv[1]
treat_minus_file = sys.argv[2]
output_file = sys.argv[3]
peak_file = sys.argv[4]


peak_text = open(peak_file,"rt")
peak_dict = {}
for p_line in peak_text:
        linelist = re.split("\t",p_line.rstrip())
        if linelist[0] not in peak_dict:
                peak_dict[linelist[0]] = []
        peak_dict[linelist[0]].append([int(linelist[1]),int(linelist[2])])
peak_text.close()

output_text = open(output_file,"wt")
output_text.write("chr_name\tsite_pos\tsite_ref\tstrand_type\tt_stopcount\tt_stoprate\tstoprate_dif_minus1\tstoprate_dif_plus1\tgene_id\n")
#get plus strand site
t_plus_text = open(treat_plus_file,"rt")
fir_t_linelist = []
sed_t_linelist = []
thr_t_linelist = []
for t_line in t_plus_text:
        #get three site info
        if not (fir_t_linelist):
                fir_t_linelist = re.split("\t",t_line)
                continue
        elif not (sed_t_linelist):
                sed_t_linelist = re.split("\t",t_line)
                continue
        elif not (thr_t_linelist):
                thr_t_linelist = re.split("\t",t_line)
        else:
                fir_t_linelist = sed_t_linelist
                sed_t_linelist = thr_t_linelist
                thr_t_linelist = re.split("\t",t_line)
        #confirm base type
        if not (sed_t_linelist[2] == "A" or sed_t_linelist[2] == "a"):
                continue
        #confirm peak region
        for peak in peak_dict[sed_t_linelist[0]]:
                if (int(sed_t_linelist[1]) > peak[0] and int(sed_t_linelist[1]) < peak[1]):
                        break
        else:
                continue
        #confirm three site in one gene
        if not (fir_t_linelist[6] == sed_t_linelist[6] and sed_t_linelist[6] == thr_t_linelist[6]):     #skip normal site
                if not (int(fir_t_linelist[1]) + int(thr_t_linelist[1]) == 2 * int(sed_t_linelist[1])): #skip consecutive site
                        continue
                elif (len(re.split(",",fir_t_linelist[6])) == 1 and len(re.split(",",thr_t_linelist[6])) == 1): #skip special case
                        continue
        #confirm stop info
        sed_t_stopcount = int(sed_t_linelist[3])
        fir_t_stoprate = float(fir_t_linelist[5])
        sed_t_stoprate = float(sed_t_linelist[5])
        thr_t_stoprate = float(thr_t_linelist[5])
        if not (sed_t_stopcount >= 20  and sed_t_stoprate >= 40):
                continue
        stoprate_diff_minus1 = (sed_t_stoprate-fir_t_stoprate)/sed_t_stoprate
        stoprate_diff_plus1 = (sed_t_stoprate-thr_t_stoprate)/sed_t_stoprate
        if not (stoprate_diff_minus1 >= 0.4 and stoprate_diff_plus1 >= 0.4):
                continue
        print(sed_t_linelist[0],sed_t_linelist[1],sed_t_linelist[2],"+",sed_t_stopcount,sed_t_stoprate,stoprate_diff_minus1,stoprate_diff_plus1,sed_t_linelist[6],sep = "\t",file = output_text)
t_plus_text.close()


#get minus strand site
t_minus_text = open(treat_minus_file,"rt")
fir_t_linelist = []
sed_t_linelist = []
thr_t_linelist = []
for t_line in t_minus_text:
        #get three site info
        if not (fir_t_linelist):
                fir_t_linelist = re.split("\t",t_line)
                continue
        elif not (sed_t_linelist):
                sed_t_linelist = re.split("\t",t_line)
                continue
        elif not (thr_t_linelist):
                thr_t_linelist = re.split("\t",t_line)
        else:
                fir_t_linelist = sed_t_linelist
                sed_t_linelist = thr_t_linelist
                thr_t_linelist = re.split("\t",t_line)
        #confirm base type
        if not (sed_t_linelist[2] == "T" or sed_t_linelist[2] == "t"):
                continue
        #confirm peak region
        for peak in peak_dict[sed_t_linelist[0]]:
                if (int(sed_t_linelist[1]) > peak[0] and int(sed_t_linelist[1]) < peak[1]):
                        break
        else:
                continue
        #confirm three site in one gene
        if not (fir_t_linelist[6] == sed_t_linelist[6] and sed_t_linelist[6] == thr_t_linelist[6]):     #skip normal site
                if not (int(fir_t_linelist[1]) + int(thr_t_linelist[1]) == 2 * int(sed_t_linelist[1])): #skip consecutive site
                        continue
                elif (len(re.split(",",fir_t_linelist[6])) == 1 and len(re.split(",",thr_t_linelist[6])) == 1): #skip special case
                        continue
        #confirm stop info
        sed_t_stopcount = int(sed_t_linelist[3])
        fir_t_stoprate = float(fir_t_linelist[5])
        sed_t_stoprate = float(sed_t_linelist[5])
        thr_t_stoprate = float(thr_t_linelist[5])
        if not (sed_t_stopcount >= 20  and sed_t_stoprate >= 40):
                continue
        stoprate_diff_minus1 = (sed_t_stoprate-fir_t_stoprate)/sed_t_stoprate
        stoprate_diff_plus1 = (sed_t_stoprate-thr_t_stoprate)/sed_t_stoprate
        if not (stoprate_diff_minus1 >= 0.4 and stoprate_diff_plus1 >= 0.4):
                continue
        print(sed_t_linelist[0],sed_t_linelist[1],sed_t_linelist[2],"-",sed_t_stopcount,sed_t_stoprate,stoprate_diff_minus1,stoprate_diff_plus1,sed_t_linelist[6],sep = "\t",file = output_text)
t_minus_text.close()

output_text.close()

