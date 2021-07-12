# qx.jl
#
# @author Laura Colbran 2021-01-14
# implementation of Qx test for polygenic selection
# as described in Berg & Coop 2014
# requires bcftools

using ArgParse
using GZip
using StatsBase
using LinearAlgebra
using SQLite

SNP_COL = :rsid #saving variable model DB uses
GENE_COL = :gene #saving variable model DB uses
WEIGHT_COL = :weight #saving variable model DB uses
BIN_COL = 3 #column in bin file that contains bin id

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
function matchSNPs(gene::String,snps::Array{String,1},bin_path::String,num_to_match::Int64)
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
    open("matched_snps.txt","w") do outf
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
        # pop_gts = pop_gts[findall(!isequal("./."),pop_gts)] # remove missing entries
        pop_gts = pop_gts[findall(!isequal(".|."),pop_gts)] # remove missing entries
        alleles = Int64[]
        for gt in pop_gts
            alleles = vcat(alleles,parse.(Int64,split(gt,"|")))
            # alleles = vcat(alleles,parse.(Int64,split(gt,"/")))
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
        # if in(chr,["1","16"]) continue end
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
                # println(index)
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
                # println(l[1:7])
                # println(alt)
                # println(alt_gt)
                snp_freq = calcFreq(ids,l[9:end],pop_path,pop_tags,flipped,alt_gt)
                if snp_freq == zeros(Float64,length(pop_tags)) println("chr$(l[1])_$(l[2])") end
                println("$(l[1])_$(l[2]): $snp_freq")
                # modify freqs so correct column has snp_freq values
                freqs[:,index] = snp_freq
                # exit()
            end
        end
    end
    return freqs
end

#organizes input parsing for Qx calculation
function QxByGene(db_path::String,match_path::String,bin_path::String,vcf_path::String,pop_path::String,pop_tags::Array{String,1},num_to_match::Int64)
    outp = "$(splitext(db_path)[1])_qx.txt"
    # snp_betas = transpose([-0.0510446249,-0.0586809545,0.3037922667])
    # snp_freqs = transpose([0.02 0.3 0.4;
    #                         0.02 0.2 0.4;
    #                         0.42 0.42 0.3;
    #                         0.04 0.14 0.3;
    #                         0.02 0.14 0.5])
    # matched_freqs = [0.02 0.3 0.4 0.03 0.05 0.0; #this is gets us a real value for Qx
    #                 0.02 0.3 0.4 0.03 0.05 0.09;
    #                 0.42 0.32 0.3 0.04 0.05 0.1;
    #                 0.04 0.04 0.3 0.03 0.03 0.06;
    #                 0.02 0.04 0.5 0.2 0.02 0.06]
    # matched_freqs = [0.02 0.3 0.4 0.0 0.05 0.03; #this is results in NaNs for F matrix onwards, and NaN Qx
    #                 0.02 0.3 0.4 0.0 0.05 0.09;
    #                 0.42 0.32 0.3 0.0 0.05 0.1;
    #                 0.04 0.04 0.3 0.0 0.03 0.06;
    #                 0.02 0.04 0.5 0.0 0.02 0.06]
    # matched_freqs = [0.02 0.3 0.4 0.03 0.05 0.01; #this breaks the cholesky factorization bc F is not positive definite
    #             0.02 0.3 0.4 0.03 0.05 0.01;
    #             0.42 0.32 0.3 0.04 0.05 0.0;
    #             0.04 0.04 0.3 0.03 0.03 0.01;
    #             0.02 0.04 0.5 0.2 0.02 0.02]
    # println(size(snp_freqs))
    # println(size(matched_freqs))
    # qx = Qx(snp_betas,snp_freqs,matched_freqs)
    # exit()
    genes = parseDB(db_path) # Dict{gene => [(id,weight)]}
    n = 1
    open(outp,"w") do outf
        write(outf,"gene\tnum_snps\tqx\n")
        for gene in keys(genes)
            println(gene)
            snp_betas = transpose([snp[2] for snp in genes[gene]])
            display([snp for snp in genes[gene]])
            println()
            println(snp_betas)
            # then pull AFs, and matched AFs
            if length(match_path) > 0
                matched_snps = readMatch(match_path,gene)
            else
                matched_snps = matchSNPs(gene,[snp[1] for snp in genes[gene]],bin_path,num_to_match)
            end
            snp_freqs = transpose(pullFreqs([snp[1] for snp in genes[gene]],vcf_path,pop_path,pop_tags))
            display(snp_freqs)
            println()
            exit()
            matched_freqs = pullFreqs(matched_snps,vcf_path,pop_path,pop_tags)
            display(matched_freqs[1:end,1:15])
            println()
            qx = Qx(snp_betas,snp_freqs,matched_freqs)
            write(outf,"$gene\t$(length(genes[gene]))\t$qx\n")
            exit()
        end
    end
end

function main()
    parsed_args = parseCommandLine()
    QxByGene(parsed_args["db"],parsed_args["matched_snps"],parsed_args["snps"],parsed_args["pop_vcfs"],parsed_args["populations"],parsed_args["pop_tags"],parsed_args["num_match"])
end

main()
