#!/isb/chumphri/anaconda/bin/python
import sys
import os.path
import numpy
import scipy
import time
from scipy import stats

def stouffers_pval(pvals):
    if (len(pvals) > 0):
        z,s_pval = scipy.stats.combine_pvalues(pvals, method="stouffer")
    else:
        s_pval = "NA"
    return (s_pval)

def append_pval(effect, pval):
    global GenePvalPos
    global GenePvalNeg

    if ( effect >= 1 ):   # effect is positive                                                                                                                                                                   
            GenePvalPos.append(pval)
    elif( effect < 1 ):  #effect is negative                                                                                                                                                                  
            GenePvalNeg.append(pval)

def clean_variable(effect, pval):
    effect = float(effect)
    pval   = float(pval)
    return (effect, pval)

def print_results():
    nullvarP_nullGeneP_p = 0
    nullvarP_sigGeneP_p  = 0
    sigvarP_nullGeneP_p  = 0
    sigvarP_sigGeneP_p   = 0
    nullvarP_nullGeneP_n = 0
    nullvarP_sigGeneP_n  = 0
    sigvarP_nullGeneP_n  = 0
    sigvarP_sigGeneP_n   = 0

    for k,p in GenePosEffect.items():
        n = GeneNegEffect[k]
        E = GeneUni[k]
        a = GeneALLEffect[k]
        result_string = [ k, p, GeneNegEffect[k], GeneUni[k], a]
        new.write("\t".join("%s" % e1 for e1 in  result_string))
        new.write("\n")
        
        #calculate sig for positive genes
        if  ((p < 0.05) and  (E == "DE")):
            sigvarP_sigGeneP_p  += 1
        elif ((p < 0.05) and  (E != "DE")):
            sigvarP_nullGeneP_p += 1
        elif ((p >= 0.05) and (E =="DE")):
            nullvarP_sigGeneP_p += 1
        else:
            nullvarP_nullGeneP_p += 1

        #calculate sig for negative genes
        if ((n < 0.05) and (E == "DE")):
            sigvarP_sigGeneP_n += 1
        elif ((n < 0.05) and (E != "DE")):
            sigvarP_nullGeneP_n += 1
        elif ((n > 0.05) and (E == "DE")):
            nullvarP_sigGeneP_n += 1
        else:
            nullvarP_nullGeneP_n += 1

    print "Direction\tnullvarP_nullGeneP\tnullvarP_sigGeneP\tsigvarP_nullGeneP\tnullvarP_nullGeneP"
    print "Pos\t" + str(nullvarP_nullGeneP_p) + "\t" + str(nullvarP_sigGeneP_p) + "\t" +  str(sigvarP_nullGeneP_p) + "\t"+ str(nullvarP_nullGeneP_p) 
    print "Neg\t" + str(nullvarP_nullGeneP_n) + "\t" + str(nullvarP_sigGeneP_n) + "\t" +  str(sigvarP_nullGeneP_n) + "\t"+ str(nullvarP_nullGeneP_n)

def print_header():
    header_string = ["Chr","PosPval", "NegPval", "ExpDE"]
    Namefile1=sys.argv[1]
    Namefile2=sys.argv[2]
    line1= "##"+ time.strftime("%d/%m/%Y")
    line2= "## Exp Diff File:"+str(Namefile1)
    line3= "## Variant File:"+str(Namefile2)
    line4=("\t".join(header_string))
    new.write(line1)
    new.write("\n")
    new.write(line2)
    new.write("\n")
    new.write(line3)
    new.write("\n")
    new.write(line4)
    new.write("\n")

f = open(sys.argv[1],"r")
next(f) # skip header
GeneUni={}
for line in f:
    info = line.split()
    pval = float(info[9])
    gene = info[10]
    gene = gene.rstrip()
    gene = gene.upper()

    if (pval < 0.1):
        val = "DE"
    else:
        val = "Same"

    GeneUni.update({gene:val})
    

var = open(sys.argv[2], "r")
next (var)

oldgene    = "default"
GenePosEffect  = {}
GeneNegEffect  = {}
GeneALLEffect  = {}
GenePvalPos    = []
GenePvalNeg    = [] 

 
for line in var:
    info = line.split()
    vargene  = info[2]
    vargene  = vargene.upper()
    vargene  = vargene.rstrip()
    if vargene in GeneUni:
        if info[3]=="NA":
            var.next()
            continue
        effect, pval  = clean_variable(info[3], info[4])
        if oldgene == vargene:
            append_pval(effect, pval)
        elif oldgene == "default":
            append_pval(effect, pval)
        else: 
            GenePosEffect.update({oldgene:stouffers_pval(GenePvalPos)})
            GeneNegEffect.update({oldgene:stouffers_pval(GenePvalNeg)})
            combinedPval = GenePvalPos + GenePvalNeg
            GeneALLEffect.update({oldgene:stouffers_pval(combinedPval)})
            append_pval(effect, pval)
            del GenePvalPos[:]
            del GenePvalNeg[:]
        
        oldgene = vargene

GenePosEffect.update({oldgene:stouffers_pval(GenePvalPos)})
GeneNegEffect.update({oldgene:stouffers_pval(GenePvalNeg)})

combinedPval = GenePvalPos + GenePvalNeg
GeneALLEffect.update({oldgene:stouffers_pval(combinedPval)})

new_fileName = os.path.basename(sys.argv[1]) + "_" + os.path.basename(sys.argv[2]) + "_" + time.strftime("%d_%m_%Y") + ".txt"

new = open(new_fileName, 'w')
print_header()
print_results()
