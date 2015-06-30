def overlap_bedfile(fh1):
    import os
    import pandas as pd
    bedfiles = ['UCSC_3_UTR_genes.bed.gz', 'UCSC_5_UTR.bed', 'UCSC_Exons.tsv.gz', 'UCSC_Introns.tsv.gz', 'TFBS_Conserved.gz'] 
    cmd = 'awk \'{OFS="\t"; print $1,$2,$2}\' ' + fh1 + ' |grep -v pos > temp1.txt'
    overlap = [] # create empty values for which to place overlap results
    columns_test = ['UTR3', 'UTR5','EXON', 'INTRON', 'TFBPS','INTRAGENIC']
    os.system(cmd)

    for newbed in bedfiles:
        dir_bed = '/users/chumphri/Reference/'
        realBed = dir_bed + newbed
        command = 'intersectBed -a temp1.txt -b ' + realBed + '  -c | awk \'{ if ($(NF)>0) print}\' | wc -l'
        process = os.popen(command).readlines()
        overlap.append(process[0].rstrip('\n'))

    ##to obtain intergenic counts, I need get the counts that equal 0 or do not overlap the wholeGenes bedfile
    intergenic = '/users/chumphri/Reference/UCSC_WholeGenes.tsv.gz'
    command = 'intersectBed -a temp1.txt -b ' + intergenic + '  -c | awk \'{ if ($(NF)==0) print}\' | wc -l'
    process = os.popen(command).readlines()
    overlap.append(process[0].rstrip('\n'))

    values = pd.DataFrame( [overlap], columns=columns_test)
    last_command = 'rm -f temp1.txt'
    os.system(last_command)
    return(values)


def grab_dir_name(directory):
    filelist=[]
    for root, dirs, files in os.walk(directory, topdown=True):
        #print dirs
        for file in files:
            if file == 'SigRes_qcb_cm.txt':
                arguement1 = os.path.join(root,file)
                filelist.append(arguement1)
    return(filelist)


def rownames(files):
    rowids = []
    for fh in files:
        r =fh.split('/')
        rowids.append(r[10])
    return(rowids)


def nucleotide_to_num(line):
    #changes multisample VCF to allele counts
    import re
    mapping = [ ('0/0', '0'), ('0/1', '1'), ('1/0', '1'), ('1/1', '2')]
    for k,v in mapping:
        line = line.rstrip('\n')
        line = line.replace(k,v)
    
    #removes hemizygous and no calls from line as well
    line = re.sub("\./\.|\./1|\./0|1/\.|0/\.", "NA", line)
    return line

def header_M_to_FAM(first):
    import re

    first = first.rstrip('\n')
    first = re.sub("M-","", first)
    furs  = first.split('\t')

    for i,k in enumerate(furs):
        furs[i] = k+"-FAM"

    new_line = "\t".join(furs)
    return new_line
