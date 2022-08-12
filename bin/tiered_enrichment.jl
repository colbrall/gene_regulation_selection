# tiered_enrichment.jl
#
# @author Laur a Colbran
# calculates  Scaled ORA across the top X genes where X is allowed to vary

using ArgParse
using DataFrames,CSV
using HypothesisTests,Plots

ENV["GKSwstype"] = "100"

function parseCommandLine()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--gene_list","-l"
            help = "path to list of genes you want to test for enrichments in. assumes genes are col 1"
            arg_type = String
        "--column","-c"
            help = "column in gene list you want to sort by. if 0, will assume it's already sorted"
            arg_type = Int64
            default = 0
        "--reverse"
            help ="to sort largest to smallest instead of smallest to largest"
            action = :store_true
        "--gene_set","-s"
            help = "path to set of genes you're testing for enrichment. assumes genes are col 1"
            arg_type = String
        "--test_points","-p"
            arg_type = Int64
            nargs = '*'
            default = [10,20,30,50,75,100,150,200,300,500,1000,2000,3000,5000,10000]
        "--out_path","-o"
            help = "path prefix for output"
            arg_type = String
            default = "./"
    end
    return parse_args(s)
end

function readSet(path::String)
    gs = Set(String[])
    open(path) do inf
        for line in eachline(inf)
            if startswith(line,"#") continue end
            gene = split(chomp(line),"\t")[1]
            push!(gs,gene)
        end
    end
    return gs
end

function enrichment(list_path::String,col::Int64,reverse::Bool,set_path::String,points::Array{Int64,1},out_dir::String)
    gene_df = CSV.read(list_path,DataFrame)
    # println(first(gene_df,10))
    if col != 0
        if reverse
            sort!(gene_df,col,rev=true)
        else
            sort!(gene_df,col)
        end
    end
    # println(first(gene_df,10))
    println("Total genes: $(nrow(gene_df))")
    g_set = readSet(set_path)
    println("Genes in set: $(length(g_set))")
    ints = nrow(gene_df[findall(x->in(x,g_set),gene_df[!,1]),:])
    prop_overall = ints/nrow(gene_df)
    println("Total Intersection: $(ints) ($(100*prop_overall)%)")
    enr = DataFrames.DataFrame(cutoffs = vcat([nrow(gene_df)],points),prop=repeat([prop_overall],length(points)+1),enr = repeat([1.0],length(points)+1),pval=repeat([-1.0],length(points)+1))
    for cutoff in points
        tmp = gene_df[1:cutoff,:]
        n_succ = nrow(tmp[findall(x->in(x,g_set),tmp[!,1]),:])
        enr[findall(x->x==cutoff,enr[:,:cutoffs]),:prop] .= n_succ/nrow(tmp)
        enr[findall(x->x==cutoff,enr[:,:cutoffs]),:enr] .= (n_succ/nrow(tmp))/prop_overall
        enr[findall(x->x==cutoff,enr[:,:cutoffs]),:pval] .= pvalue(BinomialTest(n_succ,nrow(tmp),prop_overall))
    end
    sort!(enr,:cutoffs)
    println(enr)
    delete!(enr,nrow(enr))
    enr_plot = Plots.plot(enr[:,:cutoffs],repeat([1.0],nrow(enr)),color = :grey,xlabel="Cutoff",ylabel = "Enrichment",margin=7Plots.mm,grid=false)
    enr_plot = Plots.scatter!(enr[:,:cutoffs],enr[:,:enr],legend = false,color = :black, alpha = 0.5,markersize=5)
    Plots.savefig("$(out_dir)enrplot.pdf")
end

function main()
    parsed_args = parseCommandLine()
    enrichment(parsed_args["gene_list"],parsed_args["column"],parsed_args["reverse"],parsed_args["gene_set"],parsed_args["test_points"],parsed_args["out_path"])
end

main()
