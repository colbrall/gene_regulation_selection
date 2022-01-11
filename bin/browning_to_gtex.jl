# browning_to_gtex.jl
#
# converts browning introgressed SNPs to dbSNP151 and hg38 coordinates/IDs
# works by chromosome. outputs 1 file per chromosome
#
# requires installed liftOver
# @author Laura Colbran

using ArgParse
using GZip

COMPS = Dict{String,String}("A"=>"T","C"=> "G","G"=>"C","T"=>"A")

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--indir","-i"
            help = "directory with all the input files you want to use. will search for all files with chr22 (eg)"
            arg_type = String

        "--chromosomes","-c"
            nargs='*'
            help = "chromosomes to run for. default all autosomes"
            arg_type = Int64
            default = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

        "--mapchain","-m"
            help = "map.chain file for liftOver. needs hg19 as the target, hg38 as query."
            arg_type = String

        "--snp_list","-s"
            help = "file with position and rsIDs to use. assumes GTEx freeze file."
            arg_type = String

        "--column","-l"
            help = "column from SNP list with rsID (default 7). assumes chr/pos is column 2/3."
            default=7
            arg_type = Int64

        "--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String
            default = "./"
    end
    return parse_args(s)
end

# reads Browning supplemental Files, combines into one BED file per chromosome
function combineFiles(dir_path::String,chr::Int64)
    snps = Dict{String,String}()
    for item in readdir(dir_path)
        if contains(item,"chr$(chr).")
            pop = split(item,".")[1]
            open("$(dir_path)$(item)") do inf
                for line in eachline(inf)
                    l = split(chomp(line),"\t")
                    if l[2] == "POS" continue end
                    pos = parse(Int64,l[2])
                    if "chr$(chr)\t$(pos-1)\t$(pos)" in keys(snps)
                        snps["chr$(chr)\t$(pos-1)\t$(pos)"] = "$(snps["chr$(chr)\t$(pos-1)\t$(pos)"]),$(pop)"
                    else
                        snps["chr$(chr)\t$(pos-1)\t$(pos)"] = "$(l[4])\t$(l[5])\t$(l[6])\t$(l[8])\t$(l[7])\t$(l[9])\t$(l[10])\t$(pop)"
                    end
                end
            end
        end
    end
    open("tmp_chr$(chr)_hg19.bed","w") do outf
        open("tmp_chr$(chr).anno","w") do anno
            write(outf,"#chr\tpos-1\tpos\n")
            # write(anno,"#chr\tpos-1\tpos\tref\talt\tsegment\tsprimeScore\tallele\tNMatch\tDMatch\tpopulations\n")
            for snp in keys(snps)
                write(outf,"$(snp)\n")
                write(anno,"$(snp)\t$(snps[snp])\n")
            end
        end
    end
end

# wrapper for liftOver to convery hg19 version to hg38
function liftOver(hg19_path::String,chain_path::String,chr::Int64)
    hg38_path = "tmp_chr$(chr)_hg38.bed"
    cmd = `liftOver $(hg19_path) $(chain_path) $(hg38_path) chr$(chr).unmapped`
    run(cmd)

    cmd = `wc -l chr$(chr).unmapped`
    open(cmd) do inf
        line = readlines(inf)[1]
        n = parse(Int64,split(chomp(line)," ")[1])
        if n > 0
            cmda = `cut -f 3 chr$(chr).unmapped`
            cmdb = `grep -v "Deleted"`
            run(pipeline(pipeline(cmda,cmdb),stdout="rm.txt"))


            cmda = `grep -v -f rm.txt tmp_chr$(chr).anno`
            cmdb = `cut -f 4-11`
            run(pipeline(pipeline(cmda,cmdb),stdout="anno.txt"))
        else
        cmdb = `cut -f 4-11 tmp_chr$(chr).anno`
        run(pipeline(cmdb,stdout="anno.txt"))
        end
    end

    cmda = `paste -d"\t" $(hg38_path) anno.txt`
    cmdb = `sort -n -k 3`
    run(pipeline(pipeline(cmda,cmdb),stdout="tmp_chr$(chr)_hg38_anno.bed"))
    # run(pipeline(cmda,stdout="tmp_chr$(chr)_hg38_anno.bed"))

    cmd = `sed -i '1s/^/#chr\tpos-1\tpos\tref\talt\tsegment\tsprimeScore\tallele\tNMatch\tDMatch\tpopulations\n/' tmp_chr$(chr)_hg38_anno.bed`
    run(cmd)
end

# reads in GTEx snps
function readSNPs(path::String,col::Int64,chr::Int64)
    snps = Dict{String,Array{SubString,1}}()
    GZip.open(path) do inf
        for line in eachline(inf)
            l = split(chomp(line),'\t')
            if l[2] != "chr$(chr)" continue end
            snps["$(l[2])\t$(l[3])"] = [l[4],l[5],l[col]]
        end
    end
    return snps
end

# uses GTEx freeze file to pull correct rsID, based on position, ref/alt
function convertRsIDs(hg38_path::String,snp_path::String,col::Int64,out_path::String,chr::Int64)
    gtex_snps = readSNPs(snp_path,col,chr)
    open(hg38_path) do inf
        open("$(out_path)chr$(chr)_updated.txt", "w") do outf
            write(outf,"#chr\thg38_pos\tdbSNP151_rsID\tref\talt\tsegment\tsprimeScore\tallele\tNMatch\tDMatch\tpopulations\n")
            for line in eachline(inf)
                l = split(chomp(line),"\t")
                var = "$(l[1])\t$(l[3])"
                if !in(var,keys(gtex_snps)) continue end
                #skips if either the gtex ref or alt is not present
                gtex_ref = gtex_snps[var][1]
                gtex_alt = gtex_snps[var][2]
                if (gtex_ref != l[4]) & (gtex_ref != l[5])
                    if (gtex_ref != COMPS[l[4]]) & (gtex_ref != COMPS[l[5]])
                        continue
                    end
                end
                if (gtex_alt != l[4]) & (gtex_alt != l[5]) & (gtex_alt != COMPS[l[4]]) & (gtex_alt != COMPS[l[5]])
                    continue
                end
                write(outf,"$(var)\t$(gtex_snps[var][3])\t$(join(l[4:end],"\t"))\n")
            end
        end
    end
end

function main()
    parsed_args = parseCommandLine()
    for chr in parsed_args["chromosomes"]
        combineFiles(parsed_args["indir"],chr)
        liftOver("tmp_chr$(chr)_hg19.bed",parsed_args["mapchain"],chr)
        convertRsIDs("tmp_chr$(chr)_hg38_anno.bed",parsed_args["snp_list"],parsed_args["column"],parsed_args["out_dir"],chr)
    end
end

main()
