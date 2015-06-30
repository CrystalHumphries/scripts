#########################
## connect to database ##
#########################
from impala.dbapi import connect
conn = connect(host=#####, port=####)
cur = conn.cursor()

########################################
## read in files and convert to lists ##
########################################
genes = open('list_of_genes', 'r').readlines()
gene_list = map(str.strip, genes)
ftb_mothers = open('FTB1v4_mothers.txt', 'r').readlines()
ftb_list = map(str.strip, ftb_mothers)
ptb_mothers = open('PTB1v4_mothers.txt', 'r').readlines()
ptb_list = map(str.strip, ptb_mothers)

##################
## create queries ##
##################
from impala.util import as_pandas
ptb_query = ("SELECT * "
       "FROM p7_ptb.illumina_variant as ill, public_hg19.ensembl_genes as ens "
       "WHERE ens.gene_name IN ('" + "','".join(gene_list) + "') " +
       "AND ill.sample_id IN ('" + "','".join(ptb_list) + "') " +
       "AND ens.chromosome = ill.chr AND (ill.pos >= ens.start AND ill.pos <= ens.stop)"
       )

cur.execute(ptb_query)
ptb_results = as_pandas(cur)
ptb_results.to_csv("./ptb_variants_Ill.tsv", sep="\t")

ptb_query = ("SELECT * "
                    "FROM p7_ptb.comgen_variant as com, public_hg19.ensembl_genes as ens "
                    "WHERE ens.gene_name IN ('" + "','".join(gene_list) + "') " +
                    "AND com.sample_id IN ('" + "','".join(ptb_list) + "') " +
                    "AND ens.chromosome = com.start AND (com.start >= ens.start AND com.start<= ens.stop)"
                    )

cur.execute(ptb_query)
ptb_results = as_pandas(cur)
ptb_results.to_csv("./ptb_variants_com.tsv", sep="\t")


ftb_query = ("SELECT * "
       "FROM p7_ptb.illumina_variant as ill, public_hg19.ensembl_genes as ens "
       "WHERE ens.gene_name IN ('" + "','".join(gene_list) + "') " +
       "AND ill.sample_id IN ('" + "','".join(ftb_list) + "') " +
       "AND ens.chromosome = ill.chr AND (ill.pos >= ens.start AND ill.pos <= ens.stop)"
       )

cur.execute(ftb_query)
ftb_results = as_pandas(cur)
ftb_results.to_csv("./ftb_variants_ill.tsv", sep="\t")
