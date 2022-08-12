# sel_summary.jl
#
# summarizes and plots selection statistics across sets genes
# julia 1.5.3

using ArgParse,GZip
using StatsBase
using Plots,DataFrames

COLOURS = ["#b93e0f","#a1c6d8","#5f3912","#feae02","#6b878f","#c6d737","#171b17","#fe8515","#b7c1c2","#9a5611"]
ENV["GKSwstype"] = "100" # prevents segfault in Plots
default(color=COLOURS[2],leg=false,grid=false,fontfamily="arial",alpha=0.5)

function parseCommandLine()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tables","-t"
            help = "path to results tables"
            arg_type = String
            nargs='*'
        "--genes","-g"
            arg_type = String
            nargs='*'
            help = "gene(s) to plot"
        "--column","-c"
            arg_type = Int64
            help = "column to look for genes in (so can be by name or ensembl id)"
        "--populations","-p"
            arg_type = String
            nargs='*'
            help = "Populations to calculate summary in. match column headers"
            default = ["YRI","ASW","MSL","ESN","LWK","ACB","GWD","CEU","FIN","TSI",
                        "IBS","GBR","CHB","JPT","CHS","CDX","KHV","GIH","PJL","BEB",
                        "ITU","STU","PUR","PEL","MXL","CLM"]
        "--output","-o"
            arg_type = String
            help = "prefix for output files"
            default = "./sel"
    end
    return parse_args(s)
end


function summarize(sel_files::Array{String,1},genes::Array{String,1},pops::Array{String,1},col::Int64,out_path::String)
    gene_dict = Dict{String,Dict{String,Array{Float64,1}}}()
    if isfile(genes[1])
        open(genes[1]) do inf
            for gene in eachline(inf)
                gene_dict[chomp(gene)] = Dict{String,Array{Float64,1}}()
                for pop in pops
                    gene_dict[chomp(gene)][pop] = Float64[]
                end
            end
            # println(collect(keys(gene_dict)))
        end
    else
        for gene in genes
            gene_dict[gene] = Dict{String,Array{Float64,1}}()
            for pop in pops
                gene_dict[gene][pop] = Float64[]
            end
        end
    end
    for file in sel_files
        GZip.open(file) do inf
            header = String[]
            for line in eachline(inf)
                if startswith(line,"#")
                    header = split(chomp(line),"\t")
                    continue
                end
                l = split(chomp(line),"\t")
                gene = l[col]
                if startswith(gene,"ENS")
                    gene = split(gene,".")[1]
                end
                if !in(gene,collect(keys(gene_dict))) continue end
                # println(gene)
                for pop in  pops
                    stat = l[findfirst(x->x==pop,header)]
                    # println(stat)
                    if in(stat,["NA","inf","-inf"]) continue end
                    gene_dict[gene][pop] = vcat(gene_dict[gene][pop],parse(Float64,stat))
                end
            end
        end
    end
    GZip.open("$(out_path)_gene_means.txt.gz","w") do outf
        write(outf,"#gene\t$(join(pops,"\t"))\n")
        for gene in keys(gene_dict)
            outl = "$gene"
            exists = false
            for pop in pops
                # println(gene_dict[gene][pop])
                if length(gene_dict[gene][pop]) == 0
                    outl = "$(outl)\tNA"
                else
                    exists = true
                    outl = "$(outl)\t$(round(mean(gene_dict[gene][pop]),digits=4))"
                end
            end
            if exists
                write(outf,"$outl\n")
            end
        end
    end
end

function main()
    parsed_args = parseCommandLine()
    summarize(parsed_args["tables"],parsed_args["genes"],parsed_args["populations"],parsed_args["column"],parsed_args["output"])
end

main()
