#model_summary.jl

using SQLite,DBInterface,DataFrames
using ArgParse
using StatsBase
using Seaborn

Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

# parses command-line arguments
function parseCommandLine()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--databases", "-d"
            help="path to directory with databases in it"
            arg_type = String
            required = true
        "--models","-m"
            arg_type = String
            help = "Path with models to summarize (ensembl id and tissue/database name)"
    end
    return parse_args(s)
end

function readModels(mod_path::String)
    mods = Dict{String,Array{String,1}}()
    open(mod_path) do inf
        for line in eachline(inf)
            if startswith(line,"gene") continue end
            l=split(chomp(line),"\t")
            if !haskey(mods,l[2])
                mods[l[2]] = [l[1]]
            else
                mods[l[2]] = vcat(mods[l[2]],l[1])
            end
        end
    end
    return mods
end

function summModels(db_dir::String,t_dict::Dict{String,Array{String,1}})
    summs = DataFrames.DataFrame(gene=String[],nSNPs = Int64[],meanEff=Float64[],sdEff = Float64[])
    for tiss in keys(t_dict)
        db = "$(db_dir)$(tiss).db"
        for gene in t_dict[tiss]
            q = "SELECT * FROM 'weights' WHERE gene='$(gene)'"
            snps = DataFrame(DBInterface.execute(SQLite.DB(db),q))
            push!(summs,[gene,nrow(snps),mean(snps[!,:weight]),mean_and_std(snps[!,:weight])[2]])
        end
    end
    return summs
end

function summarize(db_dir::String,mod_path::String)
    tiss_dict = readModels(mod_path)
    summ_stats = summModels(db_dir,tiss_dict)

    s_plot = distplot(summ_stats[:,:nSNPs],color=:black,kde=true)
    Seaborn.savefig("nsnps.pdf")
    clf()
    println(describe(summ_stats[:,:nSNPs]))
    s_plot = distplot(summ_stats[:,:meanEff],color=:black,kde=true)
    Seaborn.savefig("mean_effectsize.pdf")
    clf()
    println(describe(summ_stats[:,:meanEff]))
    s_plot = distplot(summ_stats[:,:sdEff],color=:black,kde=true)
    Seaborn.savefig("std_eff_size.pdf")
    clf()
    println(describe(filter(:sdEff => x -> !isnan(x),summ_stats)[:,:sdEff]))
end

function main()
    parsed_args = parseCommandLine()
    summarize(parsed_args["databases"],parsed_args["models"])
end

main()
