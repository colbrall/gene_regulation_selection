# qx_per_gene.jl
#
# @author Laura Colbran 2021-07-16
# utility script called by qx.jl. requires bcftools.
# julia 1.5

using ArgParse
using Dates
using GZip
using StatsBase
using LinearAlgebra

BIN_COL = 3 #column in bin file that contains bin id

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--gene","-g"
            help = "gene ID"
            arg_type = String
        "--snps","-s"
            help = "path to file with ascertainment SNP AFs and any other bin info. assumes snp ids are in form chr_pos_ref_alt_b38"
            arg_type = String
            default = ""
        "--pop_frqs","-f"
            help = "path to vcfs with population AFs for calculating Qx. Assumes that they are split by chr, and that you've put a * in for the chr number, and that they've been indexed by tabix"
            arg_type = String
            required = true
        "--num_match","-n"
            help = "number of neutral SNPs to pull to match each test SNP. Default 100"
            arg_type = Int64
            default = 100
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
    println("Z: $z")
    z_prime = centerScaleAFMatrix(length(z)) * z # estimated genetic values for the first M-1 populations, centered at the mean
    println("Z': $z_prime")
    va = scalingFactor(snp_betas,snp_freqs)
    println("Va: $va")
    f = neutralCovarianceMatrix(matched_freqs)
    display(f)
    println()
    c = cholesky(Hermitian(f)).L #was crashing initially due to rounding error in f making it look non-Hermitian
    display(c)
    println()
    x = 1/(sqrt(2*va))*inv(c)*z_prime
    println("x: $x")
    qx=transpose(x) * x
    println("Qx: $qx")
    return qx
end

#reads file with matched SNPs for each gene (output of matchSNPs())
function readMatch(path::String,gene::String)
    snps = String[]
    open(path) do inf
        for line in eachline(inf)
            if startswith(line,gene)
                snps = vcat(snps,"$(split(line,"\t")[2])")
            end
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

# identifies matched sets of SNPs based on bins
function matchSNPs(gene::String,snps::Array{String,1},bin_path::String,num_to_match::Int64,qx_path::String)
    target_bins = String[]
    bin_snps = Dict{String,Array{String,1}}()
    matched_snps = String[]
    open(bin_path) do inf
        for snp in snps
            s = split(snp,",")
            for line in eachline(inf)
                id = split(chomp(line),"\t")[1]
                bin = split(chomp(line),"\t")[BIN_COL]
                if id in s
                    target_bins = vcat(target_bins,bin)
                    break
                end
            end
        end
        seek(inf,0)
        # search for and save all SNPs for each bin we need. is faster to do this for small models rather than always saving everything
        for line in eachline(inf)
            id = split(chomp(line),"\t")[1]
            bin = split(chomp(line),"\t")[BIN_COL]
            if in(bin,target_bins)
                if in(bin,keys(bin_snps))
                    bin_snps[bin] = vcat(bin_snps[bin],"$id")
                else
                    bin_snps[bin] = ["$id"]
                end
            end
        end
    end
    for bin in keys(bin_snps)
        b_snps = bin_snps[bin][findall(x->!in(x,snps),bin_snps[bin])] #ensure we don't include the target snp
        matched_snps = vcat(matched_snps,sample(b_snps,num_to_match*count(isequal(bin),target_bins),replace=false))
    end
    out_path= "$(splitdir(qx_path)[1])/$(gene)_matched_snps.txt"
    open(out_path,"a") do outf
        for snp in matched_snps
            write(outf,"$gene\t$snp\n")
        end
    end
    return matched_snps
end

#pulls population freqencies for SNPs
function pullFreqs(snps::Array{String,1},frq_path::String,pop_tags::Array{String,1})
    # initiate fqcy matrix
    freqs = zeros(Float64,length(pop_tags),length(snps))
    chrs = [split(split(snp,"_")[1],"r")[2] for snp in snps]
    regions = join([join([split(split(snp,"_")[1],"r")[2],split(snp,"_")[2]],":") for snp in snps],",")
    alts = [split(snp,"_")[4] for snp in snps]
    for chr in unique(chrs)
        chr_path = replace(frq_path,"*" => chr)
        pop_inds = Int64[]
        open(chr_path) do inf
            for line in eachline(inf)
                l = split(chomp(line),"\t")
                if startswith(l[1],"CHROM")
                    for pop in pop_tag
                        pop_inds = vcat(pop_inds,findfirst(isequal(pop,l)))
                    end
                    continue
                end
                snp_ind = findfirst(isequal("$(l[1])_$(l[2])"),[join(split(snp,"_")[1:2],"_") for snp in snps])
                if typeof(snp_ind) == Nothing continue end #skip SNPs not of interest
                allele_ind = findfirst(isequal(alts[snp_ind]),[split(x,":")[1] for x in split(l[pop_inds[1]],",")])
                snp_freqs = [split(split(pop,",")[allele_ind],":")[2] for pop in l[pop_inds]]
                # modify freqs so correct column has snp_freq values
                freqs[:,snp_ind] = parse.(Float64,snp_freq)
            end
        end
    end
    indices = Int64[]
    for col in 1:size(freqs)[2]
        if freqs[:,col] == repeat([0.0],size(freqs)[1])
            indices = vcat(indices,col)
        end
    end
    return freqs,indices
end

function main()
    parsed_args = parseCommandLine()
    println(parsed_args["gene"])
    println(Dates.now())
    qx_path = parsed_args["out_path"]
    match_suffix = parsed_args["matched_snps"]
    gene = parsed_args["gene"]
    num_snps = length(parsed_args["model_snps"])
    if isfile("$(splitdir(qx_path)[1])/$(gene)$(match_suffix)") > 0
        matched_snps = readMatch("$(splitdir(qx_path)[1])/$(gene)$(match_suffix)",parsed_args["gene"])
    else
        matched_snps = matchSNPs(parsed_args["gene"],parsed_args["model_snps"],parsed_args["snps"],parsed_args["num_match"],qx_path)
    end

    all_freqs,zero_inds = pullFreqs(vcat(parsed_args["model_snps"],matched_snps),parsed_args["pop_frqs"],parsed_args["pop_tags"])

    snp_freqs = all_freqs[1:end,1:num_snps]
    snp_freqs = snp_freqs[:,setdiff(1:end, zero_inds[findall(<(num_snps+1),zero_inds)])] #remove all-zero-freq SNPs

    snp_betas = transpose(parsed_args["betas"][setdiff(1:end,zero_inds[findall(<(num_snps+1),zero_inds)])])
    snp_freqs = transpose(snp_freqs)

    matched_freqs = all_freqs[1:end,(num_snps+1):end]
    matched_freqs = matched_freqs[:,setdiff(1:end, zero_inds[findall(>(num_snps),zero_inds)] .- num_snps)]  #remove all-zero-freq SNPs

    qx = Qx(snp_betas,snp_freqs,matched_freqs)
    open(qx_path,"w") do outf
        write(outf,"$(gene)\t$(num_snps)\t$qx\n")
    end
    println(Dates.now())
end

main()
