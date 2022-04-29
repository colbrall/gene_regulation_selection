# qx_per_gene.jl
#
# @author Laura Colbran 2021-07-16
# utility script called by qx.jl.
# julia 1.5

using ArgParse
using StatsBase,LinearAlgebra

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--gene","-g"
            help = "gene ID"
            arg_type = String
        "--frequencies","-f"
            help = "file with list of IDs of matched SNPs for run. if left empty, script will match SNPs itself using --snps option"
            arg_type =  String
            default = ""
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

function readFloatDelim(f::String,delim::String)
   arr=Float64[]
   open(f) do inf
       for line in eachline(inf)
           if length(arr) == 0
               arr = parse.(Float64,split(chomp(line),delim))
           else
               arr = hcat(arr,parse.(Float64,split(chomp(line),delim)))
           end
       end
   end
   return transpose(arr)
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

function main()
    parsed_args = parseCommandLine()
    qx_path = parsed_args["out_path"]
    gene = parsed_args["gene"]
    num_snps = length(parsed_args["betas"])


    snp_betas = transpose(parsed_args["betas"])

    #read in freqs, and transpose the ones for the model snps
    all_freqs = readFloatDelim(parsed_args["frequencies"],"\t")

    snp_freqs = transpose(all_freqs[1:end,1:num_snps])
    # println(size(snp_freqs))
    matched_freqs = all_freqs[1:end,(num_snps+1):end]
    # println(size(matched_freqs))

    qx = Qx(snp_betas,snp_freqs,matched_freqs)
    open(qx_path,"a") do outf
        write(outf,"$(gene)\t$(num_snps)\t$qx\n")
    end
end

main()
