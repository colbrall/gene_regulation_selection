# qx_stats.jl
# @author Laura colbran
#
# julia1.5.3. reqires R with MASS installed

using ArgParse
using CSV,DataFrames
using StatsBase
using RCall
using LinearAlgebra,NMF
using Plots,Seaborn

COLOURS = ["#b93e0f","#a1c6d8","#5f3912","#feae02","#6b878f","#c6d737","#171b17","#fe8515","#b7c1c2","#9a5611"]
Seaborn.set(style="white", palette="muted")
HEADER = ["gene","tissue","qx","genename"]
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
            help = "number of dimensions for NMF. can be multiple. First is the k used for summary stats."
            default= [1,3,5,10,15,20,25,30,35,40]
        "--impute"
            action=:store_true
            help = "for NMF, if you want to impute missing values rather than filling in 1"
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

# returns bonferroni multiple testing correction significance threshold
function bonferroni(m::Int64)
  return 0.05/m
end

# returns Benjamini Hochberg FDR multiple testing correction significance threshold
function FDR(a::Array{Any,1})
  fdr = 0.0
  f = DataFrame(values = sort(a), rank = collect(1:length(a)))
  for i in 1:nrow(f)
    if (0.05*f[i,2]/nrow(f)) < f[i,1]
      fdr = f[i,1]
      break
    end
  end
  return fdr
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

function fillColMeans(mat::Array{Float64,2})
    for i in 1:size(mat)[2]
        col_mean = mean(mat[findall(x-> x!=-100,mat[:,i]),i])
        mat[findall(x-> x==-100,mat[:,i]),i] .= col_mean
    end
    return mat
end

# code adapted from iain's R version
function qqCI(N::Int64)
    @rput N
    R"""
    lo <- -log10(sapply(1:N,function(x) qbeta(.025,x,N-x+1)))
    up <- -log10(sapply(1:N,function(x) qbeta(.975,x,N-x+1)))
    """
    @rget lo
    @rget up
    x_vals = -log10.((collect(1:N).-0.5)./N)
    x_coords = vcat(x_vals,reverse(x_vals))
    y_coords = vcat(lo,reverse(up))
    return x_coords,y_coords
end

# uses NMF to impute missing data, at whatever k was first (so same k as that plotted)
function imputeNMF(mat::Array{Float64,2},k::Int64)
    # record indices of missing values (-100)
    missing_inds = findall(x -> x==-100,mat)
    # set them to column means for first iteration
    mat = fillColMeans(mat)
    converged = false
    num_its = 0
    e = -100
    while !converged
        # do NMF with k
        r = nnmf(mat,k)
        # multiply W and H
        impute_mat = r.W*r.H
        # replace old values at missing indices  with multiplication values
        new_mat = copy(mat)
        new_mat[missing_inds] = impute_mat[missing_inds]
        # check convergence using frobenius norm of difference
        new_e = norm(new_mat - mat)
        if abs(e-new_e) <= 1
            converged = true
        end
        num_its += 1
        mat = copy(new_mat)
        e = copy(new_e)
    end
    println("Imputation:\nNum. Iterations: $(num_its); final diff. norm: $(e)")
    return mat
end

# plot distribution and calculate empirical p-value, given list of qx values
function summaryStat(df::DataFrames.DataFrame,qx_path::String,deg_f::Int64)
    println("Qx $(summarystats(filter(!isnan,df[!,:qx])))")
    println("missing: $(length(filter(isnan,df[!,:qx])))")
    df[:,:raw_pval],df[:,:gc_pval],df[:,:gamma_p] = calcP(df[!,:qx],deg_f)
    CSV.write("$(splitext(qx_path)[1])_pvals.txt",df;delim="\t")

# calculate stats for qq plot
    df = df[findall(x -> !isnan(x),df[!,:qx]),:]
    sort!(df,:gamma_p)
    df[:,:exp] .= 0.0
    for i in 1:nrow(df)
        df[i,:exp] = -log10(i/nrow(df))
    end
    return df[:,[:gamma_p,:exp]]
end

# use NMF to reduce dims of gene x tissue matrix
# calls summaryStat for first k
function dimReduceNMF(mat_path::String,k::Array{Int64},deg_f::Int64,to_impute::Bool)
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
            qx_mat[ind-1,1:end] = parse.(Float64,replace(replace(l[2:end],"NaN"=>"-100"),"NA"=>"-100"))
        end
    end
    if to_impute
        qx_mat = imputeNMF(qx_mat,k[1])
    else
        qx_mat = fillColMeans(qx_mat)
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
        h.tick_params(axis="x",labelsize="4")
        h.tick_params(axis="y",labelsize="4")
        Seaborn.savefig("$(splitext(mat_path)[1])_k$(k[i])_H.pdf")
        clf()
    end
    W = results[1].W
    # run summaryStat for each column (ie tissue group) in W
    for i in 1:size(W)[2]
        distplot = Seaborn.kdeplot(W[:,i],color=COLOURS[i],legend = true,)
    end
    xlabel("W Matrix value")
    Seaborn.savefig("$(splitext(mat_path)[1])_k$(k[1])_dist.pdf")
    clf()
    all_p = DataFrames.DataFrame(gamma_p = Float64,exp=Float64[],group = Int64[])
    max_N = 0
    for i in 1:size(W)[2]
        println("Group $(i):")
        df = DataFrames.DataFrame(gene = genes,qx=W[:,i])
        df = summaryStat(df[df[!,:qx].!=0,:],"$(splitext(mat_path)[1])_k$(k[1])_group$(i)",deg_f)
        # df = summaryStat(df,"$(splitext(mat_path)[1])_k$(k[1])_group$(i)",deg_f)
        df[!,:group] .= i
        if nrow(df) > max_N
            max_N = nrow(df)
        end
        all_p = vcat(all_p,df)
    end

    #combined qqplot
    ci_poly_x,ci_poly_y = qqCI(max_N)
    x = [0,maximum(ci_poly_x)]
    bonf_y = repeat([-log10(bonferroni(nrow(all_p)))],2)
    fdr_y = repeat([-log10(FDR(all_p[!,:gamma_p]))],2)
    println("Bonferroni: $(bonf_y[1])")
    println("FDR: $(fdr_y[1])")
    qq = Plots.plot(Shape(ci_poly_x,ci_poly_y),color=:grey,alpha=0.5,xlabel="-log10(Expected P)",ylabel = "-log10(observed P)",margin=7Plots.mm,grid=false)
    qq = Plots.plot!(x,x,color = :grey) #,lims=(0,maximum(-log10.(results[!,:corr_pval])))
    qq = Plots.plot!(x,bonf_y,color = :red)
    qq = Plots.plot!(x,fdr_y,color = :orange)
    for i in 1:size(W)[2]
        qq = Plots.scatter!(all_p[all_p[!,:group] .== i,:exp],-log10.(all_p[all_p[!,:group] .== i,:gamma_p]), legend = :topleft,color = COLOURS[i], alpha = 0.5,markersize=3)
    end
    Plots.savefig("$(splitext(mat_path)[1])_k$(k[1])_qqplot.pdf")
end

function main()
    parsed_args = parseCommandLine()
    if parsed_args["empirical_p"]
        df = CSV.read(parsed_args["qx"],DataFrame;header=HEADER,comment="#")
        println(first(df,6))
        # plot Qx distribution
        distplot = Seaborn.kdeplot(filter(!isnan,df[!,:qx]))
        xlabel("Qx")
        Seaborn.savefig("$(parsed_args["qx"])_dist.pdf")
        clf()
        df = summaryStat(df,parsed_args["qx"],parsed_args["deg_f"])

        x = [0,maximum(df[!,:exp])]
        qq = Plots.plot(x,x,color = :grey,xlabel="-log10(Expected P)",ylabel = "-log10(observed P)",margin=7Plots.mm,grid=false) #,lims=(0,maximum(-log10.(results[!,:corr_pval])))
    	qq = Plots.scatter!(df[!,:exp],-log10.(df[!,:gamma_p]), legend = false,color = :black, alpha = 0.5,markersize=3)
        Plots.savefig("$(splitext(parsed_args["qx"])[1])_qqplot.pdf")
    end
    if parsed_args["nmf"]
        dimReduceNMF(parsed_args["qx_matrix"],parsed_args["num_groups"],parsed_args["deg_f"],parsed_args["impute"])
    end
end

main()
