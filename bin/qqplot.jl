# qqplot.jl
# given a group of p-values, plots a qqplot, and runs a genomic control
# @author Laura Colbran
# julia 1.4

using ArgParse
using DataFrames, CSV
using RCall, MultipleTesting, HypothesisTests
using Plots

ENV["GKSwstype"] = "100"
GENES = []

# parses command-line arguments
function parseCommandLine()
	s = ArgParseSettings()
	# @add_arg_table! s begin
	@add_arg_table! s begin
		"--p_file", "-f"
			help="path to file with p-values"
			arg_type = String
			required = true
		"--out", "-o"
            help = "path to write output files"
            arg_type = String
            default = "./"
		"--column","-c"
			help = "1-indexed column containing p-values to plot"
			arg_type=Int64
			required = true
        "--gc"
            help = "if you want to do genomic control"
            action=:store_true
		"--multi"
			help = "if you want to  FDR correct p-values"
			action=:store_true
		"--deg_f","-d"
			help = "number of degrees of freedom"
			arg_type = Int64
			default = 1
		"--ymax","-y"
			help = "maximum value of y axis to plot."
			default = -1
			arg_type = Int64
		end
	return parse_args(s)
end

function genomicControl(p_values::Array{Float64,1},degFree::Int64)
    @rput p_values
    @rput degFree
    R"""
        chi_stat <- qchisq(p_values,df = degFree, lower.tail=F)
        lambda <- median(chi_stat)/qchisq(0.5, df=degFree) #inflation factor
        corrected_p <- pchisq(chi_stat/lambda, df=degFree, lower.tail=F)
    """
    @rget corrected_p
    return corrected_p
end

function qqPlot(path::String,col::Int64,out::String,gc::Bool,multi::Bool,degFree::Int64,y_max::Int64)
    results = CSV.read(path, DataFrame; delim='\t')
    rename!(results,col => :pval)
    results[isnan.(results.pval), :pval] .= 1 #only happens if literally all the predictions are the exact same value (4 genes in melanocytes)
    sort!(results,:pval)
    results[!,:exp] .= 1.0

    if gc
        results[!,:corr_pval] = genomicControl(results[!,:pval],degFree)
    end

	if multi
		if gc
			results[!,:bh_Q] = adjust(results[!,:corr_pval],BenjaminiHochberg())
		else
			results[!,:bh_Q] = adjust(results[!,:pval],BenjaminiHochberg())
		end
	end

    for i in 1:nrow(results)
        results[i,:exp] = -log10(i/nrow(results))
    end
    x = [0,maximum(results[!,:exp])]

#	println(ApproximateTwoSampleKSTest(results[[x in GENES for x in results[:,1]],:corr_pval],results[[!(x in GENES) for x in results[:,1]],:corr_pval]))
#	println(pvalue(ApproximateTwoSampleKSTest(results[[x in GENES for x in results[:,1]],:corr_pval],results[[!(x in GENES) for x in results[:,1]],:corr_pval])))
	if y_max == -1
	    qq = Plots.plot(x,x,color = :grey,xlabel="-log10(Expected P)",ylabel = "-log10(observed P)",margin=7Plots.mm,grid=false) #,lims=(0,maximum(-log10.(results[!,:corr_pval])))
		qq = Plots.scatter!(results[[!in(x,GENES) for x in results[:,1]],:exp],-log10.(results[[!in(x,GENES) for x in results[:,1]],:pval]),
	            legend = false,color = :black, alpha = 0.5,markersize=3)
	    qq = Plots.scatter!(results[[x in GENES for x in results[:,1]],:exp],-log10.(results[[x in GENES for x in results[:,1]],:pval]),
	            legend = false,color = :dodgerblue, alpha = 1,markershape=:vline,markersize=9)
	    Plots.savefig("$(out)qqplot.pdf")
	else
		qq = Plots.plot(x,x,color = :grey,xlabel="-log10(Expected P)",ylabel = "-log10(observed P)",margin=7Plots.mm,grid=false,ylims=(0,y_max))
		qq = Plots.scatter!(results[[!in(x,GENES) for x in results[:,1]],:exp],-log10.(results[[!in(x,GENES) for x in results[:,1]],:pval]),
				legend = false,color = :black, alpha = 0.5,markersize=3)
		qq = Plots.scatter!(results[[x in GENES for x in results[:,1]],:exp],-log10.(results[[x in GENES for x in results[:,1]],:pval]),
				legend = false,color = :dodgerblue, alpha = 1,markershape=:vline,markersize=9)
		Plots.savefig("$(out)qqplot.pdf")
	end
	CSV.write("$(out)corrected.txt",results;delim="\t")
end

function main()
    parsed_args = parseCommandLine()
    # println(typeof(parsed_args["p_file"]),typeof(parsed_args["col"]),typeof(parsed_args["out"]),typeof(parsed_args["gc"]))
    qqPlot(parsed_args["p_file"],parsed_args["column"],parsed_args["out"],parsed_args["gc"],parsed_args["multi"],parsed_args["deg_f"],parsed_args["ymax"])
end

main()
