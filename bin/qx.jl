# qx.jl
#
# @author Laura Colbran 2021-01-14
# implementation of Qx test for polygenic selection
# as described in Berg & Coop 2014
# requires sqlite and tabix
# calls qx_per_gene.jl in a bsub job for N genes
# julia 1.5


using ArgParse
using Dates,StatsBase
using SQLite,GZip,DelimitedFiles

SNP_COL = :rsid #saving variable model DB uses
GENE_COL = :gene #saving variable model DB uses
WEIGHT_COL = :weight #saving variable model DB uses
REF_COL = :ref_allele
ALT_COL = :eff_allele
DBSNP_FILE = "/project/mathilab/colbranl/data/jti_snp_coordinates.txt.gz"

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
        "--pvalue"
            help = "if you want to calculate p-value by resampling effect sizes"
            action=:store_true
        "--resample_n","-r"
            help = "number of times to resample for p-value"
            default = 100
            arg_type = Int64
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
    return weights,[snp[2] for gene in keys(weights) for snp in weights[gene]]
end

#reads file with matched SNPs for each gene (output of matchSNPs())
function readMatch(path::String,gene::String)
    snps = String[]
    found_gene = false
    open(path) do inf
        for line in eachline(inf)
            if startswith(line,gene)
                found_gene = true
                snps = vcat(snps,"$(split(line,"\t")[2])")
                continue
            end
            if found_gene break end #keeps from reading whole file if it doesn't have to
        end
    end
    return snps
end

#pulls population freqencies for SNPs
function pullFreqs(snps::Array{String,1},frq_path::String,pop_tags::Array{String,1})
    # initiate fqcy matrix
    freqs = zeros(Float64,length(pop_tags),length(snps))
    chrs = [split(split(snp,"_")[1],"r")[2] for snp in snps]
    positions = [join(split(snp,"_")[1:2],"_") for snp in snps]
    alts = [split(snp,"_")[4] for snp in snps]
    #build tabix regions file
    tmp_file = "tmp.regions.$(Dates.now())" #naming so that parallel jobs don't collide
    open(tmp_file,"w") do outf
        for ind in 1:length(snps)
            write(outf,"$(replace(positions[ind],'_'=>'\t'))\n")
        end
    end
    for chr in unique(chrs)
        chr_path = replace(frq_path,"*" => chr)
        # id indices for pops
        pop_inds = Int64[]
        GZip.open(chr_path) do chrf
            for line in eachline(chrf)
                l = split(chomp(line),"\t")
                if startswith(l[1],"CHROM")
                    for pop in pop_tags
                        pop_inds = vcat(pop_inds,findfirst([isequal(pop,x) for x in l]))
                    end
                    break
                end
            end
        end
        #query SNPs for tabix
        command = `tabix -R $(tmp_file) $(chr_path)`
        open(command) do inf
            for line in eachline(inf)
                l = split(chomp(line),"\t")
                snp_ind = findfirst(isequal("$(l[1])_$(l[2])"),positions)
                # println(l)
                # println(snps[snp_ind])
                allele_ind = findfirst(isequal(alts[snp_ind]),[split(x,":")[1] for x in split(l[pop_inds[1]],",")])
                # println(allele_ind)
                if typeof(allele_ind) == Nothing continue end #skip if alt allele isn't present
                snp_freqs = [split(split(pop,",")[allele_ind],":")[2] for pop in l[pop_inds]]
                # modify freqs so correct column has snp_freq values
                freqs[:,snp_ind] = parse.(Float64,snp_freqs)
            end
        end
    end
    file_rm = `rm $(tmp_file)`
    run(file_rm)
    indices = Int64[]
    for col in 1:size(freqs)[2] #record indices of any SNPs we couldn't find frequencies for, so we can remove them
        if freqs[:,col] == repeat([0.0],size(freqs)[1])
            indices = vcat(indices,col)
        end
    end
    return freqs,indices
end

function runBSUB(outdir::String,commands::Array{String,1},n::Int64,genes_per_job::Int64,name::String)
    bsub_path = "$(outdir)$(name)$(n-(genes_per_job)+1)-$(n).bsub"
    open(bsub_path,"w") do outf
        write(outf,"#!/bin/bash\n")
        write(outf,"#BSUB -M 1000")
        write(outf,"#BSUB -J qx$(n-genes_per_job+1)-$(n)\n")
        write(outf,"#BSUB -o $(outdir)$(name)$(n+1-genes_per_job)-$(n).out\n")
        write(outf,"#BSUB -e $(outdir)$(name)$(n+1-genes_per_job)-$(n).err\n\n")
        for comm in commands
            write(outf,comm)
        end
    end
    run(pipeline(bsub_path,`bsub`))
end

#organizes input parsing for Qx calculation
function QxByGene(db_path::String,match_path::String,freq_path::String,pop_tags::Array{String,1},genes_per_job::Int64,pvalue::Bool,res_n::Int64)
    outdir = "$(splitext(db_path)[1])/"
    if !isdir(outdir) mkdir(outdir) end
    genes,eff_sizes = parseDB(db_path) # Dict{gene => [(id,weight)]}
    n = 1
    commands = String[]
    coord_ids = mapIDs(DBSNP_FILE)

    for gene in keys(genes)
        println(gene)
        if n%genes_per_job == 0
            runBSUB(outdir,commands,n-1,genes_per_job,"genes")
            commands = String[]
        end
        snp_arr = [snp[1] for snp in genes[gene]]
        if !occursin("_",snp_arr[1]) #if IDs are not in coordinate ID form
            snp_arr = coordinateID(coord_ids,snp_arr,[snp[3] for snp in genes[gene]],[snp[4] for snp in genes[gene]])
        end
        no_id = findall(x -> x=="",snp_arr)
        deleteat!(snp_arr,no_id)
        snps = join([snp for snp in snp_arr]," ")
        betas = ["$(snp[2])" for snp in genes[gene]]
        deleteat!(betas,no_id)

        matched_snps = readMatch(match_path,gene)
        all_freqs,zero_inds = pullFreqs(vcat(snp_arr,matched_snps),freq_path,pop_tags)
        all_freqs = all_freqs[:,setdiff(1:end, zero_inds)] #remove all-zero-freq SNPs

        betas = betas[setdiff(1:end,zero_inds[findall(<(length(betas)+1),zero_inds)])]

        #write the freq matrix so that it can be passed to qx_per_gene
        freq_file = "$(outdir)frqs/tmp.freqs.$(gene)"
        writedlm(freq_file,all_freqs)

        command = "julia ./bin/qx_per_gene.jl -g $(gene) -o $(outdir)$(gene)_model_qx.txt -b $(join(betas,' ')) -f $(freq_file)\n"
        commands = vcat(commands,[command])
        n+=1
        if pvalue
            #resample from eff_sizes and add to command
            m=1
            rn=1
            p_coms = String[]
            for samp in 1:res_n
                if n%genes_per_job == 0
                    m+=1
                    runBSUB(outdir,p_coms,rn-1,genes_per_job,"$(gene)_rand")
                    p_coms = String[]
                end
                rand_betas = sample(eff_sizes,length(betas))
                println(rand_betas)
                command = "julia ./bin/qx_per_gene.jl -g $(gene) -o $(outdir)$(gene)_rand$(m).txt -b $(join(rand_betas,' ')) -f $(freq_file)\n"
                p_coms = vcat(p_coms,[command])
                rn += 1
            end
            runBSUB(outdir,p_coms,rn-1,length(p_coms),"$(gene)_rand")
        end
    end
    runBSUB(outdir,commands,n-1,length(commands),"genes") #to catch the last few genes
end

function main()
    parsed_args = parseCommandLine()
    QxByGene(parsed_args["db"],parsed_args["matched_snps"],parsed_args["pop_freqs"],parsed_args["pop_tags"],parsed_args["genes_per_job"],parsed_args["pvalue"],parsed_args["resample_n"])
end

main()
