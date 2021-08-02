# qx_stats.jl
# @author Laura colbran
#
# julia1.5.3

using ArgParse
using CSV,DataFrames
using StatsBase
using RCall
using Seaborn,Plots

Seaborn.set(style="white", palette="muted")
HEADER = ["gene","numsnps","qx"]

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
    end
    return parse_args(s)
end

function genomicControl(chi_stat::Array{Float64,1},deg_f::Int64)
    @rput chi_stat
    @rput deg_f
    R"""
        corrected_p <- pchisq(chi_stat, df=deg_f, lower.tail=F)
        # lambda <- median(chi_stat,na.rm=TRUE)/qchisq(0.5, df=deg_f) #inflation factor
        # corrected_p <- pchisq(chi_stat/lambda, df=deg_f, lower.tail=F)
    """
    @rget corrected_p
    return corrected_p
end

# plot distribution and calculate empirical p-value, given list of qx values
function summaryStat(qx_path::String,deg_f::Int64)
    df = CSV.read(qx_path,DataFrame;header=HEADER)
    println("Qx $(summarystats(filter(!isnan,df[!,:qx])))")
    println("missing: $(length(filter(isnan,df[!,:qx])))")
    df[:,:pval] = genomicControl(df[!,:qx],deg_f)
    CSV.write("$(splitext(qx_path)[1])_rawpval.txt",df;delim="\t")

# plot Qx distribution
    dis_plot = hist(filter(!isnan,df[!,:qx]);bins=50)
    # dis_plot.set_xlabel("Qx")
    Seaborn.savefig("qx_dist.pdf")
    clf()

# plot qqplot
    df = df[findall(x -> !isnan(x),df[!,:qx]),:]
    sort!(df,:pval)
    df[:,:exp] .= 0.0
    for i in 1:nrow(df)
        df[i,:exp] = -log10(i/nrow(df))
    end
    x = [0,maximum(df[!,:exp])]

    qq = Plots.plot(x,x,color = :grey,xlabel="-log10(Expected P)",ylabel = "-log10(observed P)",margin=7Plots.mm,grid=false) #,lims=(0,maximum(-log10.(results[!,:corr_pval])))
	qq = Plots.scatter!(df[!,:exp],-log10.(df[!,:pval]), legend = false,color = :black, alpha = 0.5,markersize=3)
    Plots.savefig("qx_qqplot.pdf")
end

function main()
    parsed_args = parseCommandLine()
    summaryStat(parsed_args["qx"],parsed_args["deg_f"])
end

main()
