# qx_per_gene.jl
#
# @author Laura Colbran 2021-07-16
# utility script called by qx.jl.
# julia 1.5

using ArgParse,SQLite,GZip
using Dates
using StatsBase,LinearAlgebra,Random

RNG = MersenneTwister(1234)
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
        "--pvalue"
            help = "if you want to calculate p-value by resampling effect sizes"
            action=:store_true
        "--resample_n","-r"
            help = "number of times to resample for p-value"
            default = 100
            arg_type = Int64
        "--gene","-g"
            help = "gene id for this run"
            arg_type = String
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

# makes sure IDs are in coordinateID for the model SNPS. Also arranges it such that the ref/alt match what the model thinks they should be.
# this ensures the effect size and the frequency are referring to the same allele.
function coordinateID(ids::Dict{String,String},rsids::Array{String,1},ref_arr::Array{String,1},alt_arr::Array{String,1})
    coord_ids = String[]
    for index in 1:length(rsids)
        new_id = ""
        if rsids[index] in keys(ids) # if IDs are starting as rsIDs
            id = ids[rsids[index]]
            new_id = join([split(id,"_")[1],split(id,"_")[2],ref_arr[index],alt_arr[index],"b38"],"_") #always use ref/alt from model
        elseif occursin("_",rsids[1]) # if we just need to confirm the alleles are the correct way around
            new_id = join([split(rsids[index],"_")[1],split(rsids[index],"_")[2],ref_arr[index],alt_arr[index],"b38"],"_")
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
                allele_ind = findfirst(isequal(alts[snp_ind]),[split(x,":")[1] for x in split(l[pop_inds[1]],",")])
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

# returns matrix to center and scale Z matrix of AFs
function centerScaleAFMatrix(M::Int64)
    T = fill(-1/M,(M-1,M))# M-1 x M matrix. (M-1)/M on diagonal,-1/M everywhere else
    T[diagind(T)] .= (M-1)/M
    ### Iain's code does this, then he thinks has an option to center on subset of pops instead of all
    return T
end

# returns the scaling factor for the genetic values
function scalingFactor(betas::Transpose{Float64,Array{Float64,1}},freqs::Transpose{Float64,Array{Float64,2}})
    epsilons = [mean(freqs[i,:]) for i in 1:size(freqs)[1]]
    VA = 2*sum([betas[i]^2*epsilons[i]*(1-epsilons[i]) for i in 1:size(freqs)[1]]) ## iain has times 4 in his function for this bc he sets his va variable to equal 2*VA
    return VA
end

# returns M by M positive definite matrix describing the correlation structure of allele frequencies across populations relative to the mean/ancestral frequency.
function neutralCovarianceMatrix(freqs::Array{Float64,2})
    T = centerScaleAFMatrix(size(freqs)[1])
    S = fill(0.0, (size(freqs)[2],size(freqs)[2]))
    for i in 1:size(freqs)[2]
        epsilon = mean(freqs[:,i])
        S[i,i] = 1/(epsilon*(1-epsilon))
    end
    F=(1/(size(freqs)[2]-1)) *T*freqs*S*transpose(freqs)*transpose(T)
    return F
end

# calculates Qx statistic
function Qx(snp_betas::Transpose{Float64,Array{Float64,1}},snp_freqs::Transpose{Float64,Array{Float64,2}},c::LowerTriangular{Float64,Array{Float64,2}})
    z = [2*sum(snp_betas * snp_freqs[:,i]) for i in 1:size(snp_freqs)[2]] # calculating mean genetic values for each population
    z_prime = centerScaleAFMatrix(length(z)) * z # estimated genetic values for the first M-1 populations, centered at the mean
    va = scalingFactor(snp_betas,snp_freqs)

    x = 1/(sqrt(2*va))*inv(c)*z_prime
    qx=transpose(x) * x
    return qx
end

function readInput(gene::String,db_path::String,match_path::String,freq_path::String,pop_tags::Array{String,1},pvalue::Bool,res_n::Int64)
    outdir = "$(splitext(db_path)[1])/"
    if !isdir(outdir) mkdir(outdir) end
    genes,all_eff_sizes = parseDB(db_path) # Dict{gene => [(id,weight)]}
    n = 1
    commands = String[]
    coord_ids = mapIDs(DBSNP_FILE)

    # println(gene)
    snp_arr = [snp[1] for snp in genes[gene]]
    snp_arr = coordinateID(coord_ids,snp_arr,[snp[3] for snp in genes[gene]],[snp[4] for snp in genes[gene]])
    no_id = findall(x -> x=="",snp_arr)
    deleteat!(snp_arr,no_id)
    snps = join([snp for snp in snp_arr]," ")
    betas = parse.(Float64,["$(snp[2])" for snp in genes[gene]])
    deleteat!(betas,no_id)

    matched_snps = readMatch(match_path,gene)
    all_freqs,zero_inds = pullFreqs(vcat(snp_arr,matched_snps),freq_path,pop_tags)
    all_freqs = all_freqs[:,setdiff(1:end, zero_inds)] #remove all-zero-freq SNPs

    snp_betas = transpose(betas[setdiff(1:end,zero_inds[findall(<(length(betas)+1),zero_inds)])])
    num_snps = length(snp_betas)

    snp_freqs = transpose(all_freqs[1:end,1:num_snps])
    matched_freqs = all_freqs[1:end,(num_snps+1):end]

    f = neutralCovarianceMatrix(matched_freqs)
    c = cholesky(Hermitian(f)).L #was crashing initially due to rounding error in f making it look non-Hermitian
    qx_path = "$(outdir)$(gene)_qx.txt"
    open(qx_path,"w") do outf
        qx = Qx(snp_betas,snp_freqs,c)
        outline = "$(gene)\t$(num_snps)\t$qx"

        if pvalue
            #resample from all_eff_sizes and add to command
            all_eff_sizes = abs.(all_eff_sizes)
            signs = sign.(snp_betas)
            rand_qx = zeros(Float64,res_n)
            for samp in 1:res_n
                rand_betas = sample(RNG,all_eff_sizes,num_snps)
                rand_betas = rand_betas .* transpose(signs)
                rand_qx[samp] = Qx(transpose(rand_betas[:]),snp_freqs,c)
            end
            pval = length(filter(x-> x >= qx,rand_qx))/res_n
            if pval < 1/(res_n/10) #if it's a very small p-value, give it an order of magnitude more precision
                prec_n = 10*res_n
                rand_qx = vcat(rand_qx,zeros(Float64,prec_n-res_n))
                for samp in (res_n+1):prec_n
                    rand_betas = zeros(Float64,num_snps)
                    rand_betas[pos_inds] = sample(RNG,pos_effs,length(pos_inds))
                    rand_betas[neg_inds] = sample(RNG,neg_effs,length(neg_inds))
                    rand_qx[samp] = Qx(transpose(rand_betas[:]),snp_freqs,c)
                end
                pval = length(filter(x-> x >= qx,rand_qx))/prec_n
            end
            outline = "$(outline)\t$(pval)"
        end
        # println(outline )
        write(outf,"$(outline)\n")
    end
end

function main()
    parsed_args = parseCommandLine()
    readInput(parsed_args["gene"],parsed_args["db"],parsed_args["matched_snps"],parsed_args["pop_freqs"],parsed_args["pop_tags"],parsed_args["pvalue"],parsed_args["resample_n"])
end

main()