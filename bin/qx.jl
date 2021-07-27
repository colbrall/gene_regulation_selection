# qx.jl
#
# @author Laura Colbran 2021-01-14
# implementation of Qx test for polygenic selection
# as described in Berg & Coop 2014
# requires sqlite
# calls qx_per_gene.jl in a bsub job for each gene
# julia 1.5

# julia bin/qx.jl -d ~/1kGvar_models/dbs/Whole_Blood_full_alpha0.5_window1e6_filtered.db -v "../../data/1000g/phase3_vcfs/ALL.chr*.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" -p ../skin_pigmentation_regulation/data/1kG_1240k_dosage/1kG_phase3_superpops.txt -s ../data/gtex_1kGvariants_mafbins.txt

using ArgParse
using SQLite

SNP_COL = :rsid #saving variable model DB uses
GENE_COL = :gene #saving variable model DB uses
WEIGHT_COL = :weight #saving variable model DB uses
DBSNP_FILE = "/project/mathilab/colbranl/gene_regulation_selection/data/snp150.txt.gz"

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--snps","-s"
            help = "path to file with ascertainment SNP AFs and any other bin info. assumes snp ids are in form chr_pos_ref_alt_b38"
            arg_type = String
            default = ""
        "--db","-d"
            help = "sqlite database with genes and snp effect sizes"
            arg_type = String
            required=true
        "--pop_vcfs","-v"
            help = "path to vcfs with population AFs for calculating Qx. Assumes that they are split by chr, and that you've put a * in for the chr number, and that they've been indexed by tabix"
            arg_type = String
            required = true
        "--num_match","-n"
            help = "number of neutral SNPs to pull to match each test SNP"
            arg_type = Int64
            default = 10
        "--matched_snps","-m"
            help = "file with list of IDs of matched SNPs for run. if left empty, script will match SNPs itself using --snps option"
            arg_type =  String
            default = "_matched_snps.txt"
        "--populations","-p"
            help = "tab-delim file with population assignments for all individuals in vcf files. assumes col1 is id, col2 is pop"
            arg_type = String
            required = true
        "--pop_tags","-t"
            help = "list of population INFO tags from VCF to pull AFs for. default matches 1kG continental pops"
            arg_type = String
            default = "AFR AMR EUR EAS SAS"
        "--genes_per_job","-j"
            help = "number of genes to run per bsub job"
            arg_type = Int64
            default = 100
    end
    return parse_args(s)
end

function coordinateID(rsids::Array{String,1})
    # read in mapping file and pull correct IDs
    return join(coord_ids," ")
end

# read in database, then for each gene pull snps and effect sizes
function parseDB(path::String)
    weights = Dict{String,Array{Tuple{String,Float64},1}}() # Dict{gene => [(id,weight)]}
    # open SQLite connection, read weights table
    q = "SELECT * FROM 'weights'"
    for snp in DBInterface.execute(SQLite.DB(path),q)
        if !(snp[GENE_COL] in keys(weights))
            weights[snp[GENE_COL]] = [(snp[SNP_COL],snp[WEIGHT_COL])]
        else
            weights[snp[GENE_COL]] = vcat(weights[snp[GENE_COL]],[(snp[SNP_COL],snp[WEIGHT_COL])])
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
function QxByGene(db_path::String,match_path::String,bin_path::String,vcf_path::String,pop_path::String,pop_tags::String,num_to_match::Int64,genes_per_job::Int64)
    outdir = "$(splitext(db_path)[1])/"
    if !isdir(outdir) mkdir(outdir) end
    genes = parseDB(db_path) # Dict{gene => [(id,weight)]}
    n = 1
    commands = String[]
    for gene in keys(genes)
        if n > 10 break end
        if n%genes_per_job == 0
            # runBSUB(outdir,commands,n,genes_per_job)
            commands = String[]
        end
        snps = join([snp[1] for snp in genes[gene]]," ")
        println(snps)
        if !occursin("_",snps) #if IDs are not in coordinate ID form
            snps = coordinateID([snp[1] for snp in genes[gene]])
        end
        println(snps)
        betas = join([snp[2] for snp in genes[gene]]," ")
        command = "julia ./bin/qx_per_gene.jl -g $(gene) -l $(snps) -b $(betas) -o $(outdir)$(gene)_qx.txt -s $(bin_path) -v '$(vcf_path)' -t $(pop_tags) -p $(pop_path) -n $(num_to_match)\n"
        commands = vcat(commands,[command])
        n+=1
    end
    # runBSUB(outdir,commands,n,length(commands)) #to catch the last few genes
end

function main()
    parsed_args = parseCommandLine()
    QxByGene(parsed_args["db"],parsed_args["matched_snps"],parsed_args["snps"],parsed_args["pop_vcfs"],parsed_args["populations"],parsed_args["pop_tags"],parsed_args["num_match"],parsed_args["genes_per_job"])
end

main()
