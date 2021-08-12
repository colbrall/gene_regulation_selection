# qx_per_gene.jl
#
# @author Laura Colbran 2021-07-16
# utility script called by qx.jl.
# julia 1.5

using ArgParse
using Dates
using GZip
using StatsBase
using LinearAlgebra

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--gene","-g"
            help = "gene ID"
            arg_type = String
        "--pop_frqs","-f"
            help = "path to vcfs with population AFs for calculating Qx. Assumes that they are split by chr, and that you've put a * in for the chr number, and that they've been indexed by tabix"
            arg_type = String
            required = true
        "--matched_snps","-m"
            help = "file with list of IDs of matched SNPs for run. if left empty, script will match SNPs itself using --snps option"
            arg_type =  String
            default = ""
        "--pop_tags","-t"
            nargs='*'
            help = "list of population tags from frq file to pull AFs for. default matches 1kG continental pops"
            arg_type = String
            default = ["AFR","AMR","EUR","EAS","SAS"]
        "--model_snps","-l"
            nargs='*'
            help = "list of SNP IDs in model"
            arg_type = String
        "--betas","-b"
            nargs='*'
            help = "list of betas for SNPs in model"
            arg_type = Float64
        "--out_path","-o"
            help = "temp file to save gene to. idea is to concatenate all temp files once run is finished"
            arg_type = String
    end
    return parse_args(s)
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
function Qx(snp_betas::Transpose{Float64,Array{Float64,1}},snp_freqs::Transpose{Float64,Array{Float64,2}},matched_freqs::Array{Float64,2})
    z = [2*sum(snp_betas * snp_freqs[:,i]) for i in 1:size(snp_freqs)[2]] # calculating mean genetic values for each population
    # println("Z: $z")
    z_prime = centerScaleAFMatrix(length(z)) * z # estimated genetic values for the first M-1 populations, centered at the mean
    # println("Z': $z_prime")
    va = scalingFactor(snp_betas,snp_freqs)
    # println("Va: $va")
    f = neutralCovarianceMatrix(matched_freqs)
    # display(f)
    # println()
    c = cholesky(Hermitian(f)).L #was crashing initially due to rounding error in f making it look non-Hermitian
    # display(c)
    # println()
    x = 1/(sqrt(2*va))*inv(c)*z_prime
    # println("x: $x")
    qx=transpose(x) * x
    # println("Qx: $qx")
    return qx
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

# reads population and ids into a dictionary
function popDict(path::String)
    dict = Dict{String,Array{String,1}}() #pop => id list
    GZip.open(path) do inf
        for line in eachline(inf)
            e = split(chomp(line),"\t")
            if !in("$(e[2])",keys(dict))
                dict["$(e[2])"] = ["$(e[1])"]
            else
                dict["$(e[2])"] = vcat(dict["$(e[2])"],["$(e[1])"])
            end
        end
    end
    return dict
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
        println(chr_path)
        GZip.open(chr_path) do chrf
            for line in eachline(chrf)
                l = split(chomp(line),"\t")
                if startswith(l[1],"CHROM")
                    for pop in pop_tags
                        pop_inds = vcat(pop_inds,findfirst([isequal(pop,x) for x in l]))
                    end
                    # println(pop_inds)
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
        println("\tchr$(chr) freqs: $(Dates.now())")
    end
    file_rm = `rm $(tmp_file)`
    run(file_rm)
    indices = Int64[]
    for col in 1:size(freqs)[2] #record indices of any SNPs we couldn't find frequencies for, so we can remove them
        if freqs[:,col] == repeat([0.0],size(freqs)[1])
            indices = vcat(indices,col)
        end
    end
    println("\tupdate freqs matrix: $(Dates.now())")
    return freqs,indices
end

function main()
    parsed_args = parseCommandLine()
    println(parsed_args["gene"])
    println("Start: $(Dates.now())")
    qx_path = parsed_args["out_path"]
    gene = parsed_args["gene"]
    num_snps = length(parsed_args["model_snps"])

    matched_snps = readMatch(parsed_args["matched_snps"],parsed_args["gene"])
    println("Read Matches: $(Dates.now())")
    # println(size(matched_snps))
    all_freqs,zero_inds = pullFreqs(vcat(parsed_args["model_snps"],matched_snps),parsed_args["pop_frqs"],parsed_args["pop_tags"])
    # println(size(all_freqs))
    # println(all_freqs[1:end,1:10])

    println("pull all Freqs: $(Dates.now())")
    snp_freqs = all_freqs[1:end,1:num_snps]
    # println(size(snp_freqs))
    snp_freqs = snp_freqs[:,setdiff(1:end, zero_inds[findall(<(num_snps+1),zero_inds)])] #remove all-zero-freq SNPs
    # println(size(snp_freqs))
    # println(size(parsed_args["betas"]))
    snp_betas = transpose(parsed_args["betas"][setdiff(1:end,zero_inds[findall(<(num_snps+1),zero_inds)])])
    snp_freqs = transpose(snp_freqs)

    matched_freqs = all_freqs[1:end,(num_snps+1):end]
    # println(size(matched_freqs))
    matched_freqs = matched_freqs[:,setdiff(1:end, zero_inds[findall(>(num_snps),zero_inds)] .- num_snps)]  #remove all-zero-freq SNPs
    # println(size(matched_freqs))

    println("Split Freqs: $(Dates.now())")
    qx = Qx(snp_betas,snp_freqs,matched_freqs)
    open(qx_path,"w") do outf
        write(outf,"$(gene)\t$(num_snps)\t$qx\n")
    end
    println("Finish: $(Dates.now())")
end

main()
