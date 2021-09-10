# qx_stats.jl
# @author Laura colbran
#
# julia1.5.3. reqires R with MASS installed

using ArgParse
using CSV,DataFrames
using StatsBase
using RCall
using NMF
using Plots,Seaborn

COLOURS = [:firebrick,:tomato,:sienna,:tan,:darkblue,:dodgerblue,:darkgreen,:mediumseagreen,:black,:silver]
Seaborn.set(style="white", palette="muted")
HEADER = ["gene","numsnps","qx"]
ENV["GKSwstype"] = "100" #fixes a segfault bug in GR backend for Plots.jl: https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--qx","-q"
            help = "path to file with qx stats in it. assumes is concatenated output of qx_per_gene.jl"
            arg_type = String
		"--deg_f","-d"
			help = "number of degrees of freedom"
			arg_type = Int64
			default = 4 #matches for the default pops in qx.jl
        "--qx_matrix","-m"
            help = "path to gene x tissue matrix"
            arg_type = String
        "--empirical_p"
            help = "to calculate empirical p values for set of Qx stats in -q flag"
            action = :store_true
        "--nmf"
            help = "to use nmf to dimension reduce a gene x tissue matrix of qx stats"
            action = :store_true
        "--num_groups","-k"
            nargs='*'
            arg_type = Int64
            help = "number of dimensions for NMF. can be multiple."
            default= [1,3,5,10,15,20,25,30,35,40]
    end
    return parse_args(s)
end

function calcP(chi_stat::Array{Float64,1},deg_f::Int64)
    @rput chi_stat
    @rput deg_f
    R"""
        library(fitdistrplus)
        raw_p <- pchisq(chi_stat, df=deg_f, lower.tail=F) #raw p
        lambda <- median(chi_stat,na.rm=TRUE)/qchisq(0.5, df=deg_f) #inflation factor
        gc_p <- pchisq(chi_stat/lambda, df=deg_f, lower.tail=F)  #gc-correction
        # print(length(chi_stat))
        # print(length(chi_stat[!is.na(chi_stat)]))
        params <- fitdist(chi_stat[!is.na(chi_stat)],distr= "gamma",lower=c(0,0),start=list(shape=1,rate=1)) #fit gamma dist to qx
        gamma_p <- pgamma(chi_stat, shape=params$estimate[1], rate=params$estimate[2], lower.tail=FALSE) # gamma p
    """
    @rget raw_p
    @rget gc_p
    @rget gamma_p
    return raw_p,gc_p,gamma_p
end

# calls shell wc -l to count rows without reading whole file
function numRows(path::String)
    n = 0
    cmd = `wc -l $(path)`
    open(cmd) do inf
        for line in eachline(inf)
            n = parse(Int64,split(chomp(line)," ")[1])
        end
    end
    return n
end

# plot distribution and calculate empirical p-value, given list of qx values
function summaryStat(df::DataFrames.DataFrame,qx_path::String,deg_f::Int64)
    println("Qx $(summarystats(filter(!isnan,df[!,:qx])))")
    println("missing: $(length(filter(isnan,df[!,:qx])))")
    df[:,:raw_pval],df[:,:gc_pval],df[:,:gamma_p] = calcP(df[!,:qx],deg_f)
    CSV.write("$(splitext(qx_path)[1])_pvals.txt",df;delim="\t")

# plot qqplot
    df = df[findall(x -> !isnan(x),df[!,:qx]),:]
    sort!(df,:gamma_p)
    df[:,:exp] .= 0.0
    for i in 1:nrow(df)
        df[i,:exp] = -log10(i/nrow(df))
    end
    x = [0,maximum(df[!,:exp])]

    qq = Plots.plot(x,x,color = :grey,xlabel="-log10(Expected P)",ylabel = "-log10(observed P)",margin=7Plots.mm,grid=false) #,lims=(0,maximum(-log10.(results[!,:corr_pval])))
	qq = Plots.scatter!(df[!,:exp],-log10.(df[!,:gamma_p]), legend = false,color = :black, alpha = 0.5,markersize=3)
    Plots.savefig("$(splitext(qx_path)[1])_qqplot.pdf")
end

# use NMF to reduce dims of gene x tissue matrix
# calls summaryStat for first k
function dimReduceNMF(mat_path::String,k::Array{Int64},deg_f::Int64)
    # read in qx_matrix
    header = String[]
    genes = String[]
    qx_mat = Float64[]
    ngenes = numRows(mat_path) - 1
    open(mat_path) do inf
        for (ind,line) in enumerate(eachline(inf))
            if startswith(line,"gene")
                header = split(chomp(line),'\t')
                qx_mat = zeros(Float64,ngenes,length(header)-1)
                genes = repeat(["0"],ngenes)
                continue
            end
            l = split(chomp(line),"\t")
            genes[ind-1] = l[1]
            #replace missing data with zeroes. this is dominated by expression patterns, so while will affect the grouping, not in a bad way
            qx_mat[ind-1,1:end] = parse.(Float64,replace(replace(l[2:end],"NaN"=>"0"),"NA"=>"0"))
        end
    end
    # for a range of k, run NMF and calculate MSE
    results = NMF.Result{Float64}[]
    mses = Float64[]
    for rank in k
        println("k: $rank")
        r = nnmf(qx_mat,rank)
        mses = vcat(mses,msd(qx_mat,r.W*r.H))
        results = vcat(results,r)
    end
    # plot MSE
    # println(mses)
    mse = Plots.scatter(k,mses, xlabel="k",ylabel = "MSE",legend = false,color = :black, alpha = 0.5,markersize=3,margin=7Plots.mm,grid=false)
    Plots.savefig("$(splitext(mat_path)[1])_nmf_MSEs.pdf")
    # plot H as heatmap, since is effectively weights for tissues
    for i in 1:length(k)
        h = Seaborn.heatmap(transpose(results[i].H),xticklabels=collect(1:k[i]),yticklabels=header[2:end],square=true,cmap="OrRd")
        ylabel("Tissue")
        xlabel("Group")
        Seaborn.savefig("$(splitext(mat_path)[1])_k$(k[i])_H.pdf")
        clf()
    end
    W = results[1].W
    # run summaryStat for each column (ie tissue group) in W
    for i in 1:size(W)[2]
        distplot = Seaborn.kdeplot(W[:,i],color=COLOURS[i])
    end
    xlabel("W Matrix value")
    Seaborn.savefig("$(splitext(mat_path)[1])_k$(k[1])_dist.pdf")
    clf()
    for i in 1:size(W)[2]
        println("Group $(i):")
        df = DataFrames.DataFrame(gene = genes,qx=W[:,i])
        summaryStat(df[df[!,:qx].!=0,:],"$(splitext(mat_path)[1])_k$(k[1])_group$(i)",deg_f)
    end
end

function main()
    parsed_args = parseCommandLine()
    if parsed_args["empirical_p"]
        df = CSV.read(parsed_args["qx"],DataFrame;header=HEADER)

        # plot Qx distribution
        distplot = Seaborn.kdeplot(filter(!isnan,df[!,:qx]))
        xlabel("Qx")
        Seaborn.savefig("$(parsed_args["qx"])_dist.pdf")
        clf()
        summaryStat(df,parsed_args["qx"],parsed_args["deg_f"])
    end
    if parsed_args["nmf"]
        dimReduceNMF(parsed_args["qx_matrix"],parsed_args["num_groups"],parsed_args["deg_f"])
    end
end

main()
