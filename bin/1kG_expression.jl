# gtex_expression.jl
# @author Laura Colbran
#
# compare predictions between groups of individuals
#
# plots GTEx expression of a gene by population
#
# julia1.1

using ArgParse,GZip,SQLite
using DataFrames,CSV
using StatsBase,HypothesisTests
using Seaborn

#set plotting defaults
# default(color=:black,leg=false,grid=false,fontfamily="arial")
Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

GROUPS = ["CEU","TSI","FIN","GBR","YRI"]


function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--expression","-e"
            help = "path to expression matrix. genes.tpm.featurecounts.tsv"
            arg_type = String

        "--sample_list","-s"
            help = "1kG id to population/group"
            arg_type = String

        "--id_map","-i"
            help = "map of sample ids to 1kG ids"
            arg_type = String

        "--predictions","-p"
            help = "file of predictions to compare to"
            arg_type = String

        "--database","-d"
            help = "path to Database containing model stats"
            arg_type = String

        "--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String
            default = "./"

        "--plot"
            help = "to plot expression for each gene"
            action=:store_true

        "--genes","-n"
            nargs='*'
            help = "gene id(s) to plot."
            arg_type = String
    end
    return parse_args(s)
end

# returns dict of group=>samp_ids
function sampDict(s_path::String,groups::Array{String,1},col::Int64)
    dict = Dict{String,Array{String,1}}()
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

function pullInds(d::Array{String,1},l::Array{String,1})
    indices = Int64[]
    tmp = [findfirst(x->x==i,l) for i in d]
    for e in tmp[tmp .!= nothing] #account for there being samples in group that aren't in predictions
        append!(indices,[e])
    end
    return indices
end

function pullInds(d::Array{String,1},l::Array{SubString{String},1})
    indices = Int64[]
    tmp = [findfirst(x->x==i,l) for i in d]
    for e in tmp[tmp .!= nothing] #account for there being samples in group that aren't in predictions
        append!(indices,[e])
    end
    return indices
end

function pullInds(d::Array{SubString{String},1},l::Array{String,1})
    indices = Int64[]
    tmp = [findfirst(x->x==i,l) for i in d]
    for e in tmp[tmp .!= nothing] #account for there being samples in group that aren't in predictions
        append!(indices,[e])
    end
    return indices
end

function readGenes(file)::Array{String,1}
    genes = String[]
    open(file) do inf
        for line in eachline(inf)
            genes = vcat(genes,chomp(line))
        end
    end
    return genes
end

function idMap(file::String)
    dict = Dict{String,String}()
    open(file) do inf
        for line in eachline(inf)
            l = split(chomp(line),'\t')
            dict[l[1]] = l[2]
        end
    end
    return dict
end

function pullModels(db::String,genes::Array{String,1})
    db_name = join(split(basename(db),".")[1:(end-1)],".")
    q = "SELECT * FROM 'extra'" # WHERE gene='$(gene)'
    models = DataFrame(DBInterface.execute(SQLite.DB(db),q))
    if "n.snps.in.model" in names(models)
        models = models[[in(x,genes) for x in models[!,Symbol("gene")]],[Symbol("gene"),Symbol("n.snps.in.model"),Symbol("pred.perf.R2")]]
        rename!(models, "n.snps.in.model"=> "nsnps")
        rename!(models, "pred.perf.R2"=> "predR2")
    else
        models = models[[in(x,genes) for x in models[!,Symbol("gene")]],[Symbol("gene"),Symbol("snps"),Symbol("R2")]]
        rename!(models, "snps"=> "nsnps")
        rename!(models, "R2"=> "predR2")
    end
    return models
end

function pullBestModels(db_dir::String,genes::Array{String,1})
    mod_dict = Dict{String,Tuple{String,Int64,Float64}}()
    for item in readdir(db_dir)
        if !endswith(item,"db") continue end
        db_name = join(split(basename(item),".")[1:(end-1)],".")
        q = "SELECT * FROM 'extra'"
        tmp = DataFrame(DBInterface.execute(SQLite.DB("$(db_dir)$(item)"),q))
        tmp = tmp[[in(x,genes) for x in tmp[!,Symbol("gene")]],[Symbol("gene"),Symbol("n.snps.in.model"),Symbol("pred.perf.R2")]]
        for row in 1:nrow(tmp)
            if !haskey(mod_dict,tmp[row,:gene])
                mod_dict[tmp[row,:gene]] = (db_name,tmp[row,Symbol("n.snps.in.model")],tmp[row,Symbol("pred.perf.R2")])
            else
                if mod_dict[tmp[row,:gene]][3] < tmp[row,Symbol("pred.perf.R2")]
                    mod_dict[tmp[row,:gene]] = (db_name,tmp[row,Symbol("n.snps.in.model")],tmp[row,Symbol("pred.perf.R2")])
                end
            end
        end
    end
    models = DataFrames.DataFrame(gene=String[],nsnps=Int64[],predR2=Float64[])
    tiss = Dict{String,Array{String,1}}()
    #convert dict to df plus dict
    for gene in keys(mod_dict)
        push!(models, [gene,mod_dict[gene][2],mod_dict[gene][3]])
        if haskey(tiss,mod_dict[gene][1])
            tiss[mod_dict[gene][1]] = vcat(tiss[mod_dict[gene][1]],[gene])
        else
            tiss[mod_dict[gene][1]] = [gene]
        end
    end
    return models,tiss
end

function plotSwarm(expr_path::String,samp_path::String,id_map::String,pred_path::String,db_path::String,out_path::String,gene_ids::Array{String,1},plot::Bool)
    if isfile(gene_ids[1])
        gene_ids = readGenes(gene_ids[1])
    end
    groups = sampDict(samp_path,GROUPS,2)
    ids = idMap(id_map) # sample -> 1kG
    # println(ids)
    indices = Int64[]
    gene_res = DataFrames.DataFrame()
    open(expr_path) do f
        for line in eachline(f)
            l = split(chomp(line),'\t')
            if startswith(line,"Gene") #identify columns to pull-- those that have a sample -> 1kG mapping
                indices = pullInds(collect(keys(ids)),l)
                # println(indices)
                println("Number of samples found: $(length(indices))")
                gene_res[!,:samp_id] = [ids[id] for id in l[indices]] #convert sample ids to 1kG ids before saving
                # gene_res[!,:expr_id] = l[indices]
                continue
            end
            g = split(l[1],".")[1]
            if !in(g, gene_ids) continue end
            gene_res[!,Symbol(l[1])] = [parse(Float64,i) for i in l[indices]]
        end
    end
    # println(gene_res)

    #resolve duplicates
    gene_res = combine(groupby(gene_res,:samp_id),names(gene_res)[2:end] .=> mean,renamecols=false)
    # println(gene_res)
    println("Number of unique individuals: $(nrow(gene_res))")
    # assign groups
    if plot
        gene_res[!,:groups] = fill("Other",nrow(gene_res))
        for g in keys(groups)
            inds = pullInds(groups[g],gene_res[:,:samp_id])
            gene_res[inds,:groups] .= g
        end
    end
    # println(gene_res)

    model_stats = DataFrames.DataFrame()
    if isdir(pred_path)
        model_stats,tissues = pullBestModels(db_path,gene_ids)
        #read in for best models
        for tiss in keys(tissues)
            GZip.open("$(pred_path)$(tiss)_elasticNet0_0.5.full.gz") do f
                for line in eachline(f)
                    l = split(chomp(line),'\t')
                    if startswith(line,"gene") #identify columns to pull-- those that have a true expression value
                        indices = pullInds(gene_res[:,:samp_id],l)
                        # println("Number of samples found with predictions: $(length(indices))")
                        df_inds = pullInds(l,gene_res[:,:samp_id])
                        gene_res = gene_res[df_inds,:]
                        continue
                    end
                    g = split(l[1],".")[1]
                    if !in(g, tissues[tiss]) continue end
                    # println("$(tiss),$g")
                    gene_res[!,Symbol("$(l[1])_pred")] = [parse(Float64,i) for i in l[indices]]
                end
            end
        end
    else
        model_stats = pullModels(db_path,gene_ids)
        # println(model_stats)
        # calculate correlations
        GZip.open(pred_path) do f
            for line in eachline(f)
                l = split(chomp(line),'\t')
                if startswith(line,"gene") #identify columns to pull-- those that have a true expression value
                    df_inds = pullInds(l,gene_res[:,:samp_id])
                    # println(df_inds)
                    gene_res = gene_res[df_inds,:]
                    indices = pullInds(gene_res[:,:samp_id],l)
                    # println(indices)
                    println("Number of samples found with predictions: $(length(indices))")
                    continue
                end
                g = split(l[1],".")[1]
                if !in(g, gene_ids) continue end
                gene_res[!,Symbol("$(l[1])_pred")] = [parse(Float64,i) for i in l[indices]]
            end
        end
    end

    # println(gene_res)
    if plot
        for gene in gene_ids
            s_plot = Seaborn.violinplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
                            y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],gene])
            # s_plot = Seaborn.boxplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
            #                 y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],gene],
            #                 showfliers=false,width=0.5,color=:white)
            s_plot = swarmplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
                            y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],gene],color=:black,alpha=0.5,size=1)
            s_plot.set_title("$(gene)")
            s_plot.set_ylabel("TPM")
            Seaborn.savefig("$(out_path)$(gene)_$(join(GROUPS,'_'))_TPM.pdf")
            clf()
            for grp in GROUPS
                println("mean $(grp), $(gene): $(mean(gene_res[gene_res[:,:groups] .== grp,gene]))")
            end

            if in("$(gene)_pred",names(gene_res))
                s_plot = Seaborn.violinplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
                                y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],Symbol("$(gene)_pred")])
                s_plot = swarmplot(x=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],:groups],
                                y=gene_res[[in(i,GROUPS) for i in gene_res[!,:groups]],Symbol("$(gene)_pred")],color=:black,alpha=0.5,size=1)
                s_plot.set_title("$(gene)")
                s_plot.set_ylabel("Pred. Expr.")
                Seaborn.savefig("$(out_path)$(gene)_$(join(GROUPS,'_'))_predexp.pdf")
                clf()
                for grp in GROUPS
                    println("mean $(grp), $(gene): $(mean(gene_res[gene_res[:,:groups] .== grp,Symbol("$(gene)_pred")]))")
                end
            end
        end
    end

    rhos = DataFrames.DataFrame(gene=String[],rho=Float64[])
    open("$(out_path)corr.txt","w") do outf
        write(outf,"gene\trho\tpvalue\n")
        for gene in gene_ids
            try
                rho = corspearman(gene_res[!,Symbol("$(gene)_pred")],gene_res[!,gene])
                if !isnan(rho)
                    rhos = push!(rhos,[gene,rho])
                    p = pvalue(OneSampleZTest(atanh(rho), 1, nrow(gene_res)))
                    write(outf,"$(gene)\t$(rho)\t$(p)\n")
                else
                    write(outf,"$(gene)\tNA\tNA\n")
                end
            catch
                write(outf,"$(gene)\tNA\tNA\n")
            end
        end
    end
    # println(first(rhos,6))
    s_plot = distplot(rhos[:,:rho],color=:black,kde=true)
    Seaborn.savefig("$(out_path)corr_dist.pdf")
    clf()

    # println(first(model_stats,6))
    println("Num. Models: $(nrow(model_stats))")
    println("Num. Measured Genes: $(nrow(rhos))")
    println(describe(rhos[:,:rho]))
    rhos = innerjoin(rhos,model_stats,on=:gene)
    println("Num. with both: $(nrow(rhos))\n")
    rho = corspearman(rhos[!,:rho],rhos[!,:nsnps])
    p = pvalue(OneSampleZTest(atanh(rho), 1, nrow(rhos)))
    println("NumSNPs corr: rho=$(rho), p=$(p)")

    rho = corspearman(rhos[!,:rho],rhos[!,:predR2])
    p = pvalue(OneSampleZTest(atanh(rho), 1, nrow(rhos)))
    println("R2 corr: rho=$(rho), p=$(p)")
    s_plot = scatterplot(x=rhos[!,:predR2],y=rhos[!,:rho],color=:black,alpha=0.5,size=1)
    s_plot.set_ylabel("Expr. Rho")
    s_plot.set_ylabel("Training R2")
    Seaborn.savefig("$(out_path)expr_rho_vs_r2.pdf")
    clf()
end

function main()
    parsed_args = parseCommandLine()
    plotSwarm(parsed_args["expression"],parsed_args["sample_list"],parsed_args["id_map"],parsed_args["predictions"],
            parsed_args["database"],parsed_args["out_dir"],parsed_args["genes"],parsed_args["plot"])
end

main()
