# gtex_expression.jl
# @author Laura Colbran
#
# compare predictions between groups of individuals
#
# plots GTEx expression of a gene by population
#
# julia1.1

using ArgParse
using GZip
using DataFrames
using CSV
using Seaborn

#set plotting defaults
# default(color=:black,leg=false,grid=false,fontfamily="arial")
Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

GROUPS = ["AFR","EAS","EUR"]


function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--expression","-e"
            help = "path to expression matrix. GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
            arg_type = String

        "--attributes","-a"
            help = "file with sample information. Use to identify which IDs are which tissue. eg gtex_sample_tissues.txt"
            arg_type = String

        "--tissue","-t"
            help = "tissue to build plot for. Must match label in attributes file."
            arg_type = String

        "--sample_list","-s"
            help = "file containing sample IDs and population annotations. required."
            arg_type = String

        "--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String
            default = "./"

        "--genes","-n"
            nargs='*'
            help = "gene id(s) to plot."
            arg_type = String
    end
    return parse_args(s)
end

# returns dict of group=>samp_ids
function sampDict(s_path::String,groups::Array{String,1},col::Int64)
    dict = Dict{SubString,Array{SubString,1}}()
    for g in groups
        dict[g] = Array{SubString,1}()
    end
    dict["others"] = Array{SubString,1}()
    open(s_path) do f
        for line in eachline(f)
            l = split(chomp(line),'\t')
            if startswith(line,"#")
                println("$(l[col]) groups to compare:")
                continue
            end
            try
                append!(dict[l[col]],[l[1]])
            catch
                append!(dict["others"],[l[1]])
            end
        end
    end
    if length(groups) != 1
        delete!(dict,"others")
    end
    return dict
end

function pullInds(d::Array{String,1},l::Array{SubString{String},1})
    indices = Int64[]
    tmp = [findfirst(x->x==i,l) for i in d]
    for e in tmp[tmp .!= nothing] #account for there being samples in group that aren't in predictions
        append!(indices,[e])
    end
    return indices
end

function pullInds(d::Array{SubString,1},l::Array{SubString{String},1})
    indices = Int64[]
    tmp = [findfirst(x->x==i,l) for i in d]
    for e in tmp[tmp .!= nothing] #account for there being samples in group that aren't in predictions
        append!(indices,[e])
    end
    return indices
end

# returns array with all samples for a particular tissue
function tissSamp(file::String,tiss::String)
    samples = String[]
    open(file) do inf
        for line in eachline(inf)
            l = split(chomp(line),"\t")
            if l[2] == tiss
                append!(samples,[l[1]])
            end
        end
    end
    return samples
end

function geneTargets(file)::Array{SubString,1}
    arr = SubString[]
    if typeof(file) == Nothing return arr end
    open(file) do f
        for line in eachline(f)
            if startswith(line,"#") continue end
            append!(arr,[split(chomp(line),"\t")[4]])
        end
    end
    return arr
end

function plotSwarm(expr_path::String,samp_path::String,attr_path::String,tiss::String,out_path::String,gene_ids::Array{String,1})
    groups = sampDict(samp_path,GROUPS,2)
    target_samples = tissSamp(attr_path,tiss)
    indices = Int64[]
    println("Number of samples in tissue: $(length(target_samples))")
    gene_res = DataFrames.DataFrame()
    GZip.open("$(expr_path)") do f
        for line in eachline(f)
            l = split(chomp(line),'\t')
            if startswith(line,"Name") #identify columns to pull-- sample ID -> tissues, individual
                indices = pullInds(target_samples,l) #pull all columns for that sample
                println("Number of samples found: $(length(indices))")
                gene_res[!,:samp_id] = l[indices]
                gene_res[!,:groups] = fill("Other",nrow(gene_res))
                continue
            end
            g = split(l[1],".")[1]
            if !in(g, gene_ids) continue end
            # println(first(gene_res,6))
            # gene_res[!,Symbol(l[1])] = [parse(Float64,i) for i in l[indices]]
            gene_res[!,Symbol(l[1])] = [log(parse(Float64,i)+1,10) for i in l[indices]]
            # println(first(gene_res,6))
        end
    end
    # convert sample ID to individual GTEx ID
    gene_res[:,:samp_id] = [join(split(i,"-")[1:2],"-") for i in gene_res[:,:samp_id]]
    # assign groups
    for g in keys(groups)
        inds = pullInds(groups[g],gene_res[:,:samp_id])
        gene_res[inds,:groups] .= g
    end
    for gene in names(gene_res)[3:end]
        s_plot = Seaborn.violinplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
                        y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],gene])
        # s_plot = Seaborn.boxplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
        #                 y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],gene],showfliers=false)
        # s_plot = swarmplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
        #                 y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],gene],color=:black,alpha=0.5,size=1)
        s_plot.set_title("$(gene)")
        s_plot.set_ylabel("log10(TPM+1)")
        Seaborn.savefig("$(gene)_$(join(GROUPS,'_')).pdf")
        clf()
    end
end

function main()
    parsed_args = parseCommandLine()
    plotSwarm(parsed_args["expression"],parsed_args["sample_list"],parsed_args["attributes"],parsed_args["tissue"],parsed_args["out_dir"],parsed_args["genes"])
end

main()
