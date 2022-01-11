# qx_plots.jl
# @author Laura colbran
#
# just the plotting part of qx_stats.jl
# ie takes the written output of that and regenerates the plots (rather than needing to rerun the NMF, for example)
# julia1.5.3. requires R installation

using ArgParse
using DataFrames
using RCall
using Plots,Seaborn

COLOURS = ["#b93e0f","#a1c6d8","#5f3912","#feae02","#6b878f","#c6d737","#171b17","#fe8515","#b7c1c2","#9a5611"]
Seaborn.set(style="white", palette="muted")
HEADER = ["gene","numsnps","qx"]
ENV["GKSwstype"] = "100" #fixes a segfault bug in GR backend for Plots.jl: https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--h_matrix","-m"
            help = "path to file with H matrix info to put in heatmap"
            default = "."
            arg_type = String
        "--p_values","-p"
            help = "path(s) to p-values for each group to put in qq plot"
            arg_type = String
            nargs='*'
            default = String[]
        "--column","-c"
            help = "column for the p-value you want to plot"
            arg_type = Int64
    end
    return parse_args(s)
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

# returns bonferroni multiple testing correction significance threshold
function bonferroni(m::Int64)
  return 0.05/m
end

# returns Benjamini Hochberg FDR multiple testing correction significance threshold
function FDR(a::Array{Float64,1})
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

# plot H as heatmap, since is effectively weights for tissues
function hHeatmap(h_path::String)
    h = Seaborn.heatmap(h_mat,xticklabels=collect(1:size(h_mat)[1]),yticklabels=header[2:end],square=true,cmap="OrRd")
    ylabel("Tissue")
    xlabel("Group")
    h.tick_params(axis="x",labelsize="4")
    h.tick_params(axis="y",labelsize="4")
    Seaborn.savefig("H.pdf")
    clf()
end

function qqPlot(p_paths::Array{String,1},col::Int64)
    all_p = DataFrames.DataFrame(gamma_p = Float64[],group = Int64[],exp = Float64[])

    max_N = 0
    for i in 1:length(p_paths)
        p_values = Float64[]
        open(p_paths[i]) do inf
            for line in eachline(inf)
                if startswith(line,"gene") continue end
                p = split(chomp(line),"\t")[col]
                append!(p_values,parse(Float64,p))
            end
        end
        df = DataFrames.DataFrame(gamma_p = p_values)
        sort!(df,:gamma_p)
        df[:,:group] .= i
        df[:,:exp] .= 0.0
        for j in 1:nrow(df)
            df[j,:exp] = -log10(j/nrow(df))
        end
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

    for i in 1:length(unique(all_p[!,:group]))
        qq = Plots.scatter!(all_p[all_p[!,:group] .== i,:exp],-log10.(all_p[all_p[!,:group] .== i,:gamma_p]), legend = :topleft,color = COLOURS[i], alpha = 0.5,markersize=3)
    end
    Plots.savefig("qqplot.pdf")
end

function main()
    parsed_args = parseCommandLine()
    if parsed_args["h_matrix"] != "."
        hHeatmap(parsed_args["h_matrix"])
    end
    if length(parsed_args["p_values"]) > 0
        qqPlot(parsed_args["p_values"],parsed_args["column"])
    end
end

main()
