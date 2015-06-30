#!/usr/bin/python
import sys
import re
import tabix

#set up to grab the CADD file
url = "/proj/famgen/resources/CADD/whole_genome_SNVs.tsv.gz"
tb = tabix.open(url)


file1 = sys.argv[1]
f  = open(file1, 'r')

def get_CADD_score(end, chrom, ref, alt):
    end = int(float(end))
    start = end - 1
    records = tb.query(chrom, start, end)
    cadd_raw = "NotPresent"
    cadd_scale = "NotPresent"
    for record in records:
        if record[2]==ref and record[3]==alt:
            cadd_raw = record[4]
            cadd_scale=record[5]
    return(cadd_raw, cadd_scale)


for line in f:
    line = line.rsplit('\t')
    if 'com' in file1:
        chr_pos = 11 
        ref_pos = 4
        var_pos = ref_pos + 1
    else:
        chr_pos = 7
        ref_pos = 2
        var_pos = ref_pos + 1

    if len(line[ref_pos])==1 and len(line[var_pos])==1:
        raw, scale = get_CADD_score(line[0], str(line[chr_pos]), line[ref_pos], line[var_pos])
        result = [str(line[chr_pos]), str(line[0]), line[ref_pos], line[var_pos], raw, scale]
    else:
        result = [str(line[chr_pos]), str(line[0]), line[ref_pos], line[var_pos], "NaN", "NaN"]
    
    print('\t'.join(map(str,result)))
        

#    pos = str(line[0])+ "_" + str(line[chr_pos]) + "_" +  line[ref_pos] + "_" + line[var_pos]
#    geneCord[pos] = 'random'
#    print pos

#for CADD in sys.stdin:
#    line = CADD.split('\t')
#    if 'Institute' in CADD or 'RawScore' in CADD:
#        continue
#    pos_key = line[1] + "_" + line[0] + "_" + line[2] + "_" + line[3]
#    if pos_key in geneCord.keys():
#        print line
#        del geneCord[pos_key]

#for k,v in geneCord.items():
#    print k.replace("_", "\t")
