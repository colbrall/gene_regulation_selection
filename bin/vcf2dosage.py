# vcf2dosage.py
# Laura Colbran 5/8/17
# convert altai neanderthal vcf genome into dosage files for PrediXcan
# not generalizable as of yet
# currently only filtering out places where genotype = ./.
#
# USAGE: python vcf2dosage.py TYPE PATH/TO/VCF
# TYPE can by 1kG or update
# python 2

# dosage output:
# chr | snp_id | pos | a1 | a2 | MAF | altai

import sys
import string
import gzip

ALLELES = ["A","C","G","T"]

def thousGens(path):
    comment_char = "##"
    delim = "\t"
    with gzip.open(sys.argv[-1], 'r') as f:
        for line in f:
            if line.startswith(comment_char): continue
            l = str.split(line.strip(),delim)
            pop = l[9:]
            if line.startswith("#CHROM"):
                print "#chr\tsnp_id\tpos\ta1\ta2\tMAF\t%s" % (delim.join(pop))
                continue
            chrm = l[0]
            pos = l[1]
            id = l[2]
            a1 = l[3]
            a2 = l[4]
            maf = "."
            ac = []
            for person in pop: #iterate through genotypes
                if person == "0|0":
                    ac.append("0")
                elif person == "0|1" or person == "1|0":
                    ac.append("1")
                else:
                    ac.append("2")
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrm,id,pos,a1,a2,maf,delim.join(ac))

def update_vcf(path): #prufer 2017 VCFs
    comment_char = "##"
    delim = "\t"
    with gzip.open(sys.argv[-1], 'r') as f:
        for line in f:
            if line.startswith(comment_char): continue
            l = str.split(line.strip(),delim)
            if line.startswith("#CHROM"):
                print "#chr\tsnp_id\tpos\ta1\ta2\tMAF\t%s" % ('\t'.join(l[9:]))
                continue
            chrm = l[0]
            pos = l[1]
            id = l[2]
            a1 = l[3]
            a2 = l[4]
            maf = "."
            ac = []
            try:
                a1 = ALLELES[int(a1)-1]
                a2 = ALLELES[int(a2)-1]
            except:
                continue
            types = [str.split(x,":")[0] for x in l[9:]]
            for gt in types:
                if gt == "0/0":
                    ac.append("0")
                elif gt == "0/1" or gt == "1/0":
                    ac.append("1")
                elif gt == "1/1":
                    ac.append("2")
                else:
                    ac.append("NA")
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrm,id,pos,a1,a2,maf,'\t'.join(ac))

def old_vcf(line): #original altai, denisovan
    chrm = line[0]
    pos = line[1]
    id = line[2]
    a1 = line[3]
    a2 = line[4]
    maf = "."
    ac = "."
    for entry in str.split(line[7],";"):
        if entry.startswith("AC="):
            ac = str.split(entry,"=")[1]
    if ac != ".": #ac=. when locus didn't pass some sort of filtering
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrm,id,pos,a1,a2,maf,ac)

def main():
    if sys.argv[-2] == "1kG":
        thousGens(sys.argv[-1])
    elif sys.argv[-2] == "update":
        update_vcf(sys.argv[-1])
    else:
        comment_char = "#"
        print "#chr\tsnp_id\tpos\ta1\ta2\tMAF\taltai"
        with gzip.open(sys.argv[-1], 'r') as f:
            for line in f:
                if line.startswith(comment_char): continue
                old_vcf(str.split(line.strip(),"\t"))

if __name__ == "__main__":
    main()
