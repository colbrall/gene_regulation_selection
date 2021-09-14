# vcf2dosage.py
# Laura Colbran 5/8/17
# updated 09/14/21 to be python 3, and do mean imputing
#
# USAGE: python vcf2dosage.py PATH/TO/VCF
# confirmed functional on NYGC VCFs. double-check formatting before running others
# python 3

# dosage output:
# chr | snp_id | pos | a1 | a2 | MAF | [genotypes]

import sys
import string
import gzip

ALLELES = ["A","C","G","T"]
COMMENT_CHAR = "##"
DELIM = "\t"

def NYFC_1kG(path):
    with gzip.open(sys.argv[-1], 'r') as f:
        for line in f:
            if line.decode("utf-8").startswith(COMMENT_CHAR): continue
            l = str.split(line.decode("utf-8").strip(),DELIM)
            pop = l[9:]
            if l[0] == "#CHROM":
                print( "#chr\tsnp_id\tpos\ta1\ta2\tMAF\t%s" % (DELIM.join(pop)))
                continue
            # skip multiallelic sites
            if "," in l[4]: continue
            chrm = l[0]
            pos = l[1]
            id = l[2]
            a1 = l[3]
            a2 = l[4]
            maf = str.split(str.split(l[7],";")[1],"=")[1]
            ac = []
            for person in pop: #iterate through genotypes
                gt = str.split(person,":")[0]
                if gt == "0/0":
                    ac.append("0")
                elif gt == "0/1" or gt == "1/0":
                    ac.append("1")
                elif gt == "1/1":
                    ac.append("2")
                else: #should only catch missing GTs
                    ac.append("NA")
            # mean impute missing dosages
            mean_dosage = sum([float(e) for e in ac if e != "NA" ])/len([e for e in ac if e != "NA" ])
            ac = [str(mean_dosage) if e == "NA" else e for e in ac]
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrm,id,pos,a1,a2,maf,DELIM.join(ac)))

def main():
    NYFC_1kG(sys.argv[-1])

if __name__ == "__main__":
    main()
