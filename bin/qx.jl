# qx.jl
#
# @author Laura Colbran 2021-01-14
# implementation of Qx test for polygenic selection
# as described in Berg & Coop 2014
# requires sqlite and tabix
# calls qx_per_gene.jl in a bsub job for N genes
# julia 1.5


using ArgParse
using SQLite
using GZip

SNP_COL = :rsid #saving variable model DB uses
GENE_COL = :gene #saving variable model DB uses
WEIGHT_COL = :weight #saving variable model DB uses
REF_COL = :ref_allele
ALT_COL = :eff_allele
DBSNP_FILE = "/project/mathilab/colbranl/gene_regulation_selection/data/jti_snp_coordinates.txt.gz"

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--db","-d"
            help = "sqlite database with genes and snp effect sizes"
            arg_type = String
            required=true
        "--pop_freqs","-f"
            help = "path to files with population AFs for calculating Qx. Assumes that they are split by chr, and that you've put a * in for the chr number."
            arg_type = String
            required = true
        "--matched_snps","-m"
            help = "file with list of IDs of matched SNPs for run. if left empty, script will match SNPs itself using --snps option"
            arg_type =  String
            default = "matched_snps.txt"
        "--pop_tags","-t"
            help = "list of population INFO tags from VCF to pull AFs for. default matches 1kG continental pops"
            arg_type = String
            nargs='*'
            default = ["AFR", "AMR","EUR","EAS","SAS"]
        "--genes_per_job","-j"
            help = "number of genes to run per bsub job"
            arg_type = Int64
            default = 100
    end
    return parse_args(s)
end

function mapIDs(snp_path)
    ids = Dict{String,String}()
    GZip.open(snp_path) do f
        for line in eachline(f)
            if startswith(line,"#") continue end
            l = split(chomp(line),"\t")
            if occursin("_",l[1]) continue end #skip non-standard chromosomes
            ids[l[4]] = join([join([l[1],l[3],l[5],x,"b38"],"_") for x in split(l[7],",")[1:end-1]],",")
        end
    end
    return ids
end

function coordinateID(ids::Dict{String,String},rsids::Array{String,1},ref_arr::Array{String,1},alt_arr::Array{String,1})
    coord_ids = String[]
    for index in 1:length(rsids)
        new_id = ""
        if rsids[index] in keys(ids)
            id = ids[rsids[index]]
            if occursin(",",id)
                new_id = join([split(id,"_")[1],split(id,"_")[2],ref_arr[index],alt_arr[index],"b38"],"_")
            else
                new_id = id
            end
        end
        coord_ids = vcat(coord_ids,[new_id])
    end
    return coord_ids
end

# read in database, then for each gene pull snps and effect sizes
function parseDB(path::String)
    weights = Dict{String,Array{Tuple{String,Float64,String,String},1}}() # Dict{gene => [(id,weight)]}
    # open SQLite connection, read weights table
    q = "SELECT * FROM 'weights'"
    for snp in DBInterface.execute(SQLite.DB(path),q)
        if !(snp[GENE_COL] in keys(weights))
            weights[snp[GENE_COL]] = [(snp[SNP_COL],snp[WEIGHT_COL],snp[REF_COL],snp[ALT_COL])]
        else
            weights[snp[GENE_COL]] = vcat(weights[snp[GENE_COL]],[(snp[SNP_COL],snp[WEIGHT_COL],snp[REF_COL],snp[ALT_COL])])
        end
    end
    return weights
end

function runBSUB(outdir::String,commands::Array{String,1},n::Int64,genes_per_job::Int64)
    bsub_path = "$(outdir)genes$(n-(genes_per_job)+1)-$(n).bsub"
    open(bsub_path,"w") do outf
        write(outf,"#!/bin/bash\n")
        write(outf,"#BSUB -J qx$(n-genes_per_job+1)-$(n)\n")
        write(outf,"#BSUB -o $(outdir)genes$(n+1-genes_per_job)-$(n).out\n")
        write(outf,"#BSUB -e $(outdir)genes$(n+1-genes_per_job)-$(n).err\n\n")
        for comm in commands
            write(outf,comm)
        end
    end
    run(pipeline(bsub_path,`bsub`))
end

#organizes input parsing for Qx calculation
function QxByGene(db_path::String,match_path::String,freq_path::String,pop_tags::String,genes_per_job::Int64)
    outdir = "$(splitext(db_path)[1])/"
    if !isdir(outdir) mkdir(outdir) end
    genes = parseDB(db_path) # Dict{gene => [(id,weight)]}
    n = 1
    commands = String[]
    coord_ids = mapIDs(DBSNP_FILE)
    for gene in keys(genes)
        if n%genes_per_job == 0
            runBSUB(outdir,commands,n,genes_per_job)
            commands = String[]
        end
        snp_arr = [snp[1] for snp in genes[gene]]
        if !occursin("_",snp_arr[1]) #if IDs are not in coordinate ID form
            snp_arr = coordinateID(coord_ids,snp_arr,[snp[3] for snp in genes[gene]],[snp[4] for snp in genes[gene]])
        end
        unmatched = findall(x -> x=="",snp_arr)
        deleteat!(snp_arr,unmatched)
        snps = join([snp for snp in snp_arr]," ")
        betas = [snp[2] for snp in genes[gene]]
        deleteat!(betas,unmatched)
        command = "julia ./bin/qx_per_gene.jl -g $(gene) -l $(snps) -b $(join(betas,' ')) -o $(outdir)$(gene)_qx.txt -f '$(freq_path)' -m $(match_path) -t $(pop_tags)\n"
        # println(command)
        # exit()
        commands = vcat(commands,[command])
        n+=1
    end
    runBSUB(outdir,commands,n,length(commands)) #to catch the last few genes
end

function main()
    parsed_args = parseCommandLine()
    QxByGene(parsed_args["db"],parsed_args["matched_snps"],parsed_args["pop_freqs"],join(parsed_args["pop_tags"]," "),parsed_args["genes_per_job"])
end

main()
