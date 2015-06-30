__author__ = 'chumphri'

VCF_files = {
    'F'  :  '/isb/chumphri/PTBVarPheno/Overlap_Prune/Fathers5076_variants',
    'FAM' : '/isb/chumphri/PTBVarPheno/Overlap_Prune/Family5076_variants',
    'M'   : '/isb/chumphri/PTBVarPheno/Overlap_Prune/Mothers5076_variants',
    'NB'  : '/isb/chumphri/PTBVarPheno/Overlap_Prune/Newborns5076_variants'
}

def create_new(file, vcf):
    location = []
    f = open(file, 'r')
    for line in f:
        line = line.rstrip('\n')
        chr,pos = line.split('\t')[0:2]
        loc = chr + "_" + pos
        location.append(loc)
    no_snps = len(location)
    print "Opening the temp vcf file..."
    target = open("temp_vcf", 'w')

    with open(vcf, 'r') as fp:
        header  = next(fp)
        target.write(header)

        for line in fp:
            line = line.rstrip('\n')
            loc  = line.split('\t')
            if (loc[2] in location):
                target.write(line)
                target.write('\n')
    target.close()
    return('temp_vcf', no_snps)

def sample_file(f):
    import re
    r = f.split('/')[-1].split('_')[2]
    matching = [s for s in sample_list if r in s]
    file_id = matching[0]
    print "Opening the temp samples file..."
    target = open("temp_samples.txt", 'w')
    with open(file_id, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            new = re.sub('-[A-z]{1,3}', '', line)
            target.write(new)
            target.write('\n')
    target.close()
    return("temp_samples.txt", r)

def create_bed(bed, temp_vcf):
    print "Opening the temp samples file..."
    target = open("temp.bed", 'w')
    import os.path
    if os.path.isfile(bed):
        bed_real = bed
    else:
        command = 'cut -f3 ' + temp_vcf + ' > fake_bed'
        os.system(command)
        bed_real = 'fake_bed'
    count = 0
    with open(bed_real, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            new  = line.split('_')
            new_line = new[0]+'\t' + new[1] + '\t' + new[1]
            target.write(new_line)
            target.write('\n')
            count+=1
    target.close()
    return("temp.bed", count)

def overlap_bedfile(fh1):
    import os
    import pandas as pd
    intersectBed = '/isb/chumphri/BEDtools/bin/intersectBed'
    bedfiles = ['UCSC_3_UTR_genes.bed.gz', 'UCSC_5_UTR.bed', 'UCSC_Exons.tsv.gz', 'UCSC_Introns.tsv.gz', 'TFBS_Conserved.gz']
    overlap = [] # create empty values for which to place overlap results
    columns_test = ['UTR3', 'UTR5','EXON', 'INTRON', 'TFBPS','INTRAGENIC']

    for newbed in bedfiles:
        dir_bed = '/isb/chumphri/Reference/'
        realBed = dir_bed + newbed
        command = intersectBed + ' -a ' + fh1 + ' -b ' + realBed + '  -c | awk \'{ if ($(NF)>0) print}\' | wc -l'
        process = os.popen(command).readlines()
        overlap.append(process[0].rstrip('\n'))

    ##to obtain intergenic counts, I need get the counts that equal 0 or do not overlap the wholeGenes bedfile
    intergenic = '/isb/chumphri/Reference/UCSC_WholeGenes.tsv.gz'
    inter_comm = intersectBed + ' -a ' + fh1 + ' -b ' + intergenic + '  -c | awk \'{ if ($(NF)==0) print}\' | wc -l'
    process = os.popen(inter_comm).readlines()
    overlap.append(process[0].rstrip('\n'))

    values = pd.DataFrame([overlap], columns=columns_test)
    last_command = 'rm -f ' + fh1
    os.system(last_command)
    return(values)

import glob
files =  glob.glob("../all_sig_var/*.txt")
sample_list = glob.glob("/isb/chumphri/PTBVarPheno/SamplePhenotypes/*txt")

import os
import pandas as pd

#create empty data frame
columns = ['UTR3', 'UTR5','EXON', 'INTRON', 'TFBPS','INTRAGENIC', 'PHENO', 'MEMBER', 'OrigSNP', 'Prunned']
main_df = pd.DataFrame(columns=columns)

for f in files:
    type = f.split('/')[-1].split('_')[1]
    temp = VCF_files[type]
    temp_vcf, orig = create_new(f, temp)
    samples, pheno  = sample_file(f)
    random_output = 'random'
    output_temp = '/isb/chumphri/' + random_output
    command = 'plink --vcf {0} --vcf-half-call m --indep-pairwise 50 5 0.8 --keep-fam {1} --out {2}'.format(temp_vcf,
                                                                                                            samples,
                                                                                                            output_temp)
    os.system(command)
    new_file = '/isb/chumphri/random.prune.in'
    bedfile, prun = create_bed(new_file, temp_vcf)
    temp_df = overlap_bedfile(bedfile)
    remove_files='rm -f ' + new_file
    temp_df['PHENO']=pheno
    temp_df['MEMBER']=type
    temp_df['OrigSNP']=orig
    temp_df['Prunned']=prun
    print temp_df
    results = [ main_df, temp_df]
    main_df = pd.concat(results)
    print "Completed file: " + str(files.index(f)) + " of " + str(len(files))

print main_df
main_df.to_csv("final_stats_dataframe", sep="\t")
