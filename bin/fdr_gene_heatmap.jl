# gene_heatmap.jl
# @author Laura colbran
#
# julia1.5.3

using GZip
using Statistics,StatsBase,Random
using Seaborn


Seaborn.set(style="white", palette="muted")
COL_MAP = "coolwarm"
COLOURS = ["#b93e0f","#a1c6d8","#5f3912","#feae02","#6b878f","#c6d737","#171b17","#fe8515","#b7c1c2","#9a5611"]
NUM_PERMS = 10000

# ENV["GKSwstype"] = "100" #fixes a segfault bug in GR backend for Plots.jl: https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988

# returns dict of group=>samp_ids
function sampDict(s_path::String,groups::Array{String,1})
    dict = Dict{SubString,Array{SubString,1}}()
    for g in groups
        dict[g] = Array{SubString,1}()
    end
    dict["others"] = Array{SubString,1}()
    open(s_path) do f
        for line in eachline(f)
            l = split(chomp(line),'\t')
            if startswith(line,"#")
                println("$(l[2]) groups to compare:")
                continue
            end
            try
                append!(dict[l[2]],[l[1]])
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

function pullInds(d::Dict{SubString,Array{SubString,1}},l::Array{SubString{String},1})
    dict = Dict{SubString,Array{Int64,1}}()
    for g in keys(d)
        tmp = [findfirst(x->x==i,l) for i in d[g]]
        dict[g] = tmp[tmp .!= nothing] #account for there being samples in group that aren't in predictions
    end
    return dict
end

function idMap(file::String)
    dict = Dict{String,String}()
    open(file) do inf
        for line in eachline(inf)
            l = split(chomp(line),'\t')
            dict[l[2]] = l[1]
        end
    end
    return dict
end

function genoID(file::String)
    dict = Dict{String,String}()
    open(file) do inf
        for line in eachline(inf)
            l = split(chomp(line),'\t')
            dict[l[1]] = l[2]
        end
    end
    return dict
end

function pullInds(d::Dict{SubString,Array{SubString,1}},l::Array{SubString{String},1},id_path::String)
    dict = Dict{SubString,Array{Int64,1}}()
    ids = idMap(id_path) # 1kG -> sample
    for grp in keys(d)
        dict[grp] = [findfirst(x->x==ids[i],l) for i in filter(x->in(x,collect(keys(ids))),d[grp])]
    end
    return dict
end

function calcMatEmpP(expr::Array{Float64,2},pred::Array{Float64,2},num_perms::Int64)
    rhos = [corspearman(expr[i,:],pred[i,:]) for i in 1:size(expr,1)]
    true_rho = sum([rhos[i] for i in 1:size(expr,1) if length(filter(x->!isnan(x),expr[i,:])) != 0])
    num_greater = 0
    for i in 1:num_perms
        shuffled = [corspearman(shuffle(expr[i,:]),shuffle(pred[i,:])) for i in 1:size(expr,1)]
        perm_rho = sum([shuffled[i] for i in 1:size(expr,1) if length(filter(x->!isnan(x),expr[i,:])) != 0])
        if abs(perm_rho) >= abs(true_rho) num_greater += 1 end
    end
    return true_rho,num_greater/num_perms
end

function main()
    gene_order = ["ENSG00000115677","ENSG00000145979","ENSG00000130590","ENSG00000166263",
                "ENSG00000116903","ENSG00000142621","ENSG00000084072","ENSG00000048140"]
    gene_names = ["HDLBP","TBC1D7","SAMD10","STXBP4","EXOC8","FHAD1","PPIE","TSPAN17"] #effPerm FDR<0.05 genes
    # gene_order = ["ENSG00000103510","ENSG00000167394","ENSG00000169896","ENSG00000099365",
    #             "ENSG00000167395"]
    # gene_names = ["KAT8","ZNF668","ITGAM","STX1B", "ZNF646"] #KAT8 peak genes that replicate (FDR < 0.1) in 1kG and HGDP
    # gene_order = ["ENSG00000188322","ENSG00000197165"]
    # gene_names = ["SBK1","SULT1A2"] #sbk1 peak genes that replicate
    # gene_order = ["ENSG00000177414","ENSG00000117697","ENSG00000138030","ENSG00000152128","ENSG00000174123",
    #             "ENSG00000204967","ENSG00000179344","ENSG00000127399","ENSG00000253958","ENSG00000160886","ENSG00000151746",
    #             "ENSG00000149485","ENSG00000188243","ENSG00000140386","ENSG00000103510","ENSG00000140995",
    #             "ENSG00000141736","ENSG00000141452","ENSG00000104972","ENSG00000100412"]
    # gene_names = ["UBE2U","NSL1","KHK","TMEM163","TLR10","PCDHA4","HLA-DQB1","LRRC61","CLDN23","LY6K",
    #             "BICD1","FADS1","COMMD6","SCAPER","KAT8","DEF8","ERBB2","RMC1","LILRB1","ACO2"] #gamma top gene from each peak
    # acc_id_path = "../../data/1000g/geuvadis_expression/run_accession_to_1kG_id.txt"
    # preds_dir = "results/1kG_v8_JTI_predictions/"
    # preds_suff = "_elasticNet0_0.5.full.gz"
    # preds_order = repeat(["BestModels"],length(gene_names))
    # samp_path = "../data/1kG_pops/1kG_populations.txt"
    # expr_path = "../../data/1000g/geuvadis_expression/genes.tpm.featurecounts.tsv"
    # pop_order= ["YRI","ASW","MSL","ESN","LWK","ACB","GWD","CEU","FIN","TSI","IBS","GBR",
    #                 "CHB","JPT","CHS","CDX","KHV","GIH","PJL","BEB","ITU","STU",
    #                 "PUR","PEL","MXL","CLM"]
    # pop_label_colors= vcat(repeat(["#946bb3"],7),vcat(repeat(["#4f7ac9"],5),
    #     vcat(repeat(["#d48c63"],5),vcat(repeat(["#c86f6f"],5),repeat(["#74c46e"],4)))))

    preds_dir = "results/HGDP_v8_JTI_predictions/"
    preds_suff = "_elasticNet0_0.5.full.gz"
    preds_order = repeat(["BestModels"],length(gene_names))
    samp_path = "../data/HGDP_pops/hgdp_pops.txt"
    expr_path = "../../data/HGDP/merged_sample_gene_fpkm_table.txt"
    pop_order= ["AFRICA","EURO","M_EAST","E_ASIA","CS_ASIA","OCEAN","AMER"]
    pop_label_colors= ["#946bb3","#4f7ac9","#4fc9c5","#d48c63","#c86f6f","#c8c06f","#74c46e"]

    groups = sampDict(samp_path,pop_order) # group -> sample IDs
    samples = String[]
    genePredDict = Dict{String,Array{Float64,1}}() #to hold predictions
    for g in gene_order
        genePredDict[g] = [-100.0]
    end
    indices = Dict{String,Array{Int64,1}}() #group -> indices of its samples
    pred_samples = String[]
    for tiss in unique(preds_order)
        println(tiss)
        GZip.open("$(preds_dir)$(tiss)$(preds_suff)") do f
            for line in eachline(f)
                l = split(chomp(line),'\t')
                if startswith(line,"gene")
                    pred_samples = l[2:end]
                    indices = pullInds(groups,pred_samples)
                    continue
                end
                genes = gene_order[findall(x->x==tiss,preds_order)]
                if !in(l[1],genes) continue end
                # println(gene_names[findall(x->x==l[1],gene_order)])
                genePredDict[l[1]] = [parse(Float64,i) for i in l[2:end]]
            end
        end
    end
    # println(indices)

    # assemble matrix of tiss x pop in correct order
    mat = repeat([-100.0],length(gene_order),length(pop_order))
    for i in 1:length(gene_order)
        if genePredDict[gene_order[i]] == [-100.0]
            continue
        else
            for j in 1:length(pop_order)
                mat[i,j] = median(genePredDict[gene_order[i]][indices[pop_order[j]]])
            end
            mat[i,:] = (mat[i,:] .- mean(mat[i,:])) ./ std(mat[i,:]) #standardize plot
        end
    end
    # println(mat)

    # first, plot an empty heatmap with the right axis labels, cbar, etc.
    h = Seaborn.heatmap(mat,xticklabels=pop_order,yticklabels=gene_names,square=true,cmap=COL_MAP,center = 0)
    [t.set_color(i) for (i,t) in zip(pop_label_colors,h.xaxis.get_ticklabels())]
    h.tick_params(axis="x",labelsize="8")
    h.tick_params(axis="y",labelsize="8")
    ylabel("Gene")
    xlabel("Population")
    Seaborn.savefig("top_FDR_preds_heatmap.pdf")
    clf()

    expr_indices = Dict{String,Array{Int64,1}}() #group -> indices of its samples
    expr_samples = String[]
    geneExpDict = Dict{String,Array{Float64,1}}() #to hold predictions
    for g in gene_order
        geneExpDict[g] = [-100.0]
    end
    start_ind = 1
    open(expr_path) do f
        for line in eachline(f)
            l = split(chomp(line),'\t')
            if startswith(line,"Gene") #identify columns to pull-- those that have a sample -> 1kG mapping
                start_ind = 2
                sampMap = genoID(acc_id_path)
                expr_samples = [sampMap[sample] for sample in l[start_ind:end]]
                expr_indices = pullInds(groups,l[start_ind:end],acc_id_path)
                continue
            elseif startswith(line,"tracking_id") #identify for HGDP file
                start_ind = 5
                expr_samples = l[start_ind:end]
                expr_indices = pullInds(groups,expr_samples)
                continue
            end
            g = split(l[1],".")[1]
            if !in(g, gene_order) continue end
            geneExpDict[g] = [parse(Float64,i) for i in l[start_ind:end]]
        end
    end
    # println(expr_indices)
    # println(geneExpDict)

    # combine expr and preds
    tmp_pred = repeat([-100.0],size(mat)[1],size(mat)[2])
    tmp_expr = repeat([-100.0],size(mat)[1],size(mat)[2])
    plot_pops = String[] #to pull right pop_labels for combined
    for i in 1:length(gene_order)
        if geneExpDict[gene_order[i]] == [-100.0] #if the gene's expression is missing
            continue
        else
            for j in 1:length(pop_order)
                if length(expr_indices[pop_order[j]]) == 0 continue end #if there are no individuals in this pop with expression
                # filter to be same individuals
                plot_pops = vcat(plot_pops,pop_order[j])
                pred_inds = intersect(indices[pop_order[j]],findall(x->in(x,expr_samples),pred_samples)) #from indices, ones that are in expr sample
                expr_inds = intersect(expr_indices[pop_order[j]],findall(x->in(x,pred_samples),expr_samples)) #from expr_inds, ones that are in pred sample
                tmp_expr[i,j] = median(geneExpDict[gene_order[i]][expr_inds])
                tmp_pred[i,j] = median(genePredDict[gene_order[i]][pred_inds])
            end
            real_vals = filter(x->x!=-100,tmp_expr[i,:])
            tmp_expr[i,:] = [ifelse(exp != -100.0,(exp - mean(real_vals)) / std(real_vals),-100.0) for exp in tmp_expr[i,:]]
            real_vals = filter(x->x!=-100,tmp_pred[i,:])
            tmp_pred[i,:] = [ifelse(exp != -100.0,(exp - mean(real_vals)) / std(real_vals),-100.0) for exp in tmp_pred[i,:]]
        end
    end
    tmp_pred = hcat([tmp_pred[:,i] for i in 1:size(tmp_pred)[2] if tmp_pred[:,i] != repeat([-100.0],size(tmp_pred)[1])]...)
    tmp_expr = hcat([tmp_expr[:,i] for i in 1:size(tmp_expr)[2] if tmp_expr[:,i] != repeat([-100.0],size(tmp_expr)[1])]...)

    rho,p = calcMatEmpP(tmp_expr,tmp_pred,NUM_PERMS)
    println("Sum of Rhos: $(rho) (p = $(p))")

    comb_mat = hcat(hcat(tmp_pred,repeat([-100.0],size(tmp_pred)[1]),tmp_expr))
    # println(comb_mat)

    plot_pops = unique(plot_pops)
    comb_x = vcat(vcat(["$(x)_Pred" for x in plot_pops],[""]),["$(x)_Expr" for x in plot_pops])
    mask_inds = repeat([false],size(comb_mat)[1],size(comb_mat)[2]) # clear previous mask
    mask_inds[findall(x->x==-100.0,comb_mat)] .= true #mask missing tissues

    h = Seaborn.heatmap(comb_mat,xticklabels=comb_x,yticklabels=gene_names,square=true,cmap=COL_MAP,mask=mask_inds,center = 0)
    h.tick_params(axis="x",labelsize="8")
    h.tick_params(axis="y",labelsize="8")
    ylabel("Gene")
    xlabel("Population")
    Seaborn.savefig("top_FDR_combined_heatmap.pdf")
    clf()
    exit()


    samd10_expr = Float64[]
    samd10_pred = Float64[]
    TBC1D7_expr = Float64[]
    TBC1D7_pred = Float64[]
    for j in 1:length(pop_order)
        pred_inds = intersect(sort(indices[pop_order[j]]),findall(x->in(x,expr_samples),pred_samples)) #from indices, ones that are in expr sample
        expr_inds = intersect(sort(expr_indices[pop_order[j]]),findall(x->in(x,pred_samples),expr_samples)) #from expr_inds, ones that are in pred sample
        samd10_expr = vcat(samd10_expr,geneExpDict["ENSG00000130590"][expr_inds])
        samd10_pred = vcat(samd10_pred,genePredDict["ENSG00000130590"][pred_inds])
        TBC1D7_expr = vcat(TBC1D7_expr,geneExpDict["ENSG00000145979"][expr_inds])
        TBC1D7_pred = vcat(TBC1D7_pred,genePredDict["ENSG00000145979"][pred_inds])
    end

    real_vals = filter(x->x!=-100,samd10_expr)
    samd10_expr = [ifelse(exp != -100.0,(exp - mean(real_vals)) / std(real_vals),-100.0) for exp in samd10_expr]
    real_vals = filter(x->x!=-100,samd10_pred)
    samd10_pred = [ifelse(exp != -100.0,(exp - mean(real_vals)) / std(real_vals),-100.0) for exp in samd10_pred]

    real_vals = filter(x->x!=-100,TBC1D7_expr)
    TBC1D7_expr = [ifelse(exp != -100.0,(exp - mean(real_vals)) / std(real_vals),-100.0) for exp in TBC1D7_expr]
    real_vals = filter(x->x!=-100,TBC1D7_pred)
    TBC1D7_pred = [ifelse(exp != -100.0,(exp - mean(real_vals)) / std(real_vals),-100.0) for exp in TBC1D7_pred]

    p_plot = Seaborn.lineplot(x=[minimum(vcat(samd10_expr,samd10_pred)),maximum(vcat(samd10_expr,samd10_pred))],y=[minimum(vcat(samd10_expr,samd10_pred)),maximum(vcat(samd10_expr,samd10_pred))],color=:blue)
    p_plot = Seaborn.scatter(samd10_expr,samd10_pred,color=:black,alpha=0.4)
    Seaborn.xlabel("SAMD10 Observed")
    Seaborn.ylabel("SAMD10 Predicted")
    Seaborn.savefig("samd10_obspred_scatter.pdf")
    clf()

    p_plot = Seaborn.lineplot(x=[minimum(vcat(TBC1D7_expr,TBC1D7_pred)),maximum(vcat(TBC1D7_expr,TBC1D7_pred))],y=[minimum(vcat(TBC1D7_expr,TBC1D7_pred)),maximum(vcat(TBC1D7_expr,TBC1D7_pred))],color=:blue)
    p_plot = Seaborn.scatter(TBC1D7_expr,TBC1D7_pred,color=:black,alpha=0.4)

    Seaborn.xlabel("TBC1D7 Observed")
    Seaborn.ylabel("TBC1D7 Predicted")
    Seaborn.savefig("tbc1d7_obspred_scatter.pdf")
    clf()
end

main()
