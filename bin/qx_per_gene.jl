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
GT_DELIM = "|" #splitter for genotypes in VCF

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
        "--pop_vcfs","-v"
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
        "--populations","-p"
            help = "tab-delim file with population assignments for all individuals in vcf files. assumes col1 is id, col2 is pop"
            arg_type = String
            required = true
        "--pop_tags","-t"
            nargs='*'
            help = "list of population INFO tags from VCF to pull AFs for. default matches 1kG continental pops"
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
            for line in eachline(inf)
                id = split(chomp(line),"\t")[1]
                bin = split(chomp(line),"\t")[BIN_COL]
                if snp == id
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

#calculates population-level frequencies
function calcFreq(ids::Array{String,1},gts::Array{String,1},pop_path::String,pop_tags::Array{String,1},flipped::Bool,alt_gt::Int64)
    freq = zeros(Float64,length(pop_tags))
    pops = popDict(pop_path)
    i = 1
    for pop in pop_tags
        pop_gts = [split(gt,":")[1] for gt in gts[findall(ind -> in(ind,pops[pop]),ids)]] #pull gts for inds in this pop
        pop_gts = pop_gts[findall(!isequal(".$(GT_DELIM)."),pop_gts)] # remove missing entries
        alleles = Int64[]
        for gt in pop_gts
            alleles = vcat(alleles,parse.(Int64,split(gt,GT_DELIM)))
        end
        # handle multiallelic sites
        alleles[findall(!isequal(alt_gt),alleles)] .= 0
        f = sum(alleles)/(2*length(pop_gts))
        if flipped
            f = 1.0-f
        end
        freq[i] = f
        i += 1
    end
    return freq
end

#pulls population freqencies for SNPs
function pullFreqs(snps::Array{String,1},vcf_path::String,pop_path::String,pop_tags::Array{String,1})
    # initiate fqcy matrix
    freqs = zeros(Float64,length(pop_tags),length(snps))
    chrs = [split(split(snp,"_")[1],"r")[2] for snp in snps]
    # regions = join([join(split(snp,"_")[1:2],":") for snp in snps],",")
    regions = join([join([split(split(snp,"_")[1],"r")[2],split(snp,"_")[2]],":") for snp in snps],",")
    refs = [split(snp,"_")[3] for snp in snps]
    alts = [split(snp,"_")[4] for snp in snps]
    for chr in unique(chrs)
        chr_path = replace(vcf_path,"*" => chr)
        ids = String[]
        command = `bcftools view $chr_path -r $regions`
        open(command) do inf
            for line in eachline(inf)
                flipped = false
                if startswith(line,"##") continue end
                l = ["$x" for x in split(chomp(line),"\t")]
                if startswith(l[1],"#")
                    ids = l[9:end]
                    continue
                end
                index = findfirst(isequal("chr$(l[1])_$(l[2])"),[join(split(snp,"_")[1:2],"_") for snp in snps])
                if typeof(index) == Nothing continue end #skip extra indels that get picked up if there's an indel that overlaps the target
                alt = alts[index]
                if in(refs[index],split(l[5],",")) #if ref is backwards
                    flipped = true
                    alt = refs[index]
                end
                if !in(alt,split(l[5],",")) #skip if different alt
                    continue
                end
                alt_gt = findfirst(isequal(alt),split(l[5],","))
                snp_freq = calcFreq(ids,l[9:end],pop_path,pop_tags,flipped,alt_gt)
                if snp_freq == zeros(Float64,length(pop_tags)) println("chr$(l[1])_$(l[2])") end
                # modify freqs so correct column has snp_freq values
                freqs[:,index] = snp_freq
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
    if isfile("$(splitdir(qx_path)[1])/$(gene)$(match_suffix)") > 0
        matched_snps = readMatch("$(splitdir(qx_path)[1])/$(gene)$(match_suffix)",parsed_args["gene"])
    else
        matched_snps = matchSNPs(parsed_args["gene"],parsed_args["model_snps"],parsed_args["snps"],parsed_args["num_match"],qx_path)
    end

    snp_freqs,zero_inds = pullFreqs(parsed_args["model_snps"],parsed_args["pop_vcfs"],parsed_args["populations"],parsed_args["pop_tags"])
    snp_freqs = snp_freqs[:,setdiff(1:end, zero_inds)] #remove all-zero-freq SNPs

    snp_betas = transpose(parsed_args["betas"][setdiff(1:end,zero_inds)])
    snp_freqs = transpose(snp_freqs)

    matched_freqs,zero_inds = pullFreqs(matched_snps,parsed_args["pop_vcfs"],parsed_args["populations"],parsed_args["pop_tags"])
    matched_freqs = matched_freqs[:,setdiff(1:end, zero_inds)]  #remove all-zero-freq SNPs

    qx = Qx(snp_betas,snp_freqs,matched_freqs)
    open(qx_path,"w") do outf
        num_snps = length(parsed_args["model_snps"])
        gene = parsed_args["gene"]
        write(outf,"$(gene)\t$(num_snps)\t$qx\n")
    end
    println(Dates.now())
end

main()
