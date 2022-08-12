using GZip,SQLite,DataFrames
using StatsBase,HypothesisTests
using Seaborn

Seaborn.set(style="white", palette="muted")
COL_MAP = "coolwarm"

NUM_OTHER_SCORES = 6
ORDER = ["Best_Qx_P","gene_length","R2","NSNPs","O/E_LoF","PhyloP_100way","iHS","nSL"]
SCORE_DIRS = ["200kb"] #,"20kb","2kb","2Mb"

function readQx()
    scores = Dict{String,Array{Float64,1}}() #gene -> [qx,ihs*SCORE_DIRS,nsl*SCORE_DIRS]. take the largest across tissues/pops
    open("results/best_model_qx/all_qx_matchPosNeg_FDR.txt") do qx
      for line in eachline(qx)
        if startswith(line, "#gene") continue end
        l = split(chomp(line),"\t")
        # qx = l[3]
        qx = l[4]
        # println("$(l[1]) $qx")
        if startswith(qx,"N")
            scores[l[1]] = repeat([-Inf],NUM_OTHER_SCORES+(2*length(SCORE_DIRS)))
        else
            scores[l[1]] = vcat([parse(Float64,qx)],repeat([-Inf],NUM_OTHER_SCORES -1 + 2*length(SCORE_DIRS)))
        end
      end
    end
    return scores
end

function readiHSnSL(scores::Dict{String,Array{Float64,1}})
    for dir in 1:length(SCORE_DIRS)
        for chr in 1:22
            GZip.open("results/selection/gene_summaries/$(SCORE_DIRS[dir])/ihs_chr$(chr)_all_gene_means.txt.gz") do ihs
                for line in eachline(ihs)
                    if startswith(line, "#") continue end
                    l = split(chomp(line),"\t")
                    ihs = abs.(parse.(Float64,[x for x in l[2:end] if x!="NA"])) #since sign doesn't matter for significance
                    if haskey(scores,l[1])
                        scores[l[1]][NUM_OTHER_SCORES+dir] = maximum(ihs)
                    end
                end
            end
        end

        for chr in 1:22
            GZip.open("results/selection/gene_summaries/$(SCORE_DIRS[dir])/nsl_chr$(chr)_all_gene_means.txt.gz") do nsl
                for line in eachline(nsl)
                    if startswith(line, "#") continue end
                    l = split(chomp(line),"\t")
                    nsl = abs.(parse.(Float64,[x for x in l[2:end] if x!="NA"])) #since sign doesn't matter for significance
                    if haskey(scores,l[1])
                        scores[l[1]][NUM_OTHER_SCORES+dir+length(SCORE_DIRS)] = maximum(nsl)
                    end
                end
            end
        end
    end
    return scores
end

function readLength(scores::Dict{String,Array{Float64,1}})
    open("../data/gencode.v26.GRCh38.all_genes.bed") do inf
        for line in eachline(inf)
            l = split(chomp(line),'\t')
            gene = split(l[5],".")[1]
            if haskey(scores,gene)
                scores[gene][2] = parse(Float64,l[3]) - parse(Float64,l[2])
            end
        end
    end
    return scores
end

function readModels(scores::Dict{String,Array{Float64,1}})
    db_path = "../data/zhou2020_JTI_models/Best_Models_ProteinCoding.db"
        q = "SELECT * FROM 'extra'" # WHERE gene='$(gene)'
        models = DataFrame(DBInterface.execute(SQLite.DB(db_path),q))
        for row in 1:nrow(models)
            if !haskey(scores,models[row,:gene]) continue end
            if models[row,Symbol("R2")] > scores[models[row,:gene]][4]
                scores[models[row,:gene]][3] = models[row,Symbol("R2")]
            end
            if models[row,Symbol("snps")] > scores[models[row,:gene]][5]
                scores[models[row,:gene]][4] = models[row,Symbol("snps")]
            end
        end
    return scores
end

function readOELoF(scores::Dict{String,Array{Float64,1}})
    GZip.open("../data/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz") do inf
        for line in eachline(inf)
            if startswith(line,"gene") continue end
            l = split(chomp(line),'\t')
            if (haskey(scores,l[64])) && (l[24] != "NA")
                # scores[l[64]][2] = parse(Float64,l[66]) #to pull CDE length
                scores[l[64]][5] = parse(Float64,l[24])
            end
        end
    end
    return scores
end

function readPhyloP(scores::Dict{String,Array{Float64,1}})
    open("data/phyloP_avg200kb_genes.txt") do inf
        for line in eachline(inf)
            l = split(chomp(line),"\t")
            gene = split(l[1],".")[1]
            if (haskey(scores,gene))
                scores[gene][6] = parse(Float64,l[6])
            end
        end
    end
    return scores
end

function printCorr(scores::Dict{String,Array{Float64,1}})
    #assemble into array, remove genes that are missing a score
    all_genes = collect(keys(scores))
    keep_inds = Array{Int64,1}[]
    score_array = repeat([0.0],length(all_genes),NUM_OTHER_SCORES+2*length(SCORE_DIRS)) #gene x score
    for gene in 1:length(all_genes)
        if !in(-Inf,scores[all_genes[gene]])
            keep_inds = vcat(keep_inds,[gene])
            score_array[gene,:] = scores[all_genes[gene]]
        end
    end
    println("Num of genes total: $(length(all_genes))")
    all_genes = all_genes[keep_inds]
    score_array = score_array[keep_inds,:]
    println("Num of genes with all scores: $(length(all_genes))")

    sorted_genes = all_genes[sortperm(score_array[:,1],rev=true)]
    println("$(ORDER[1]) vs $(ORDER[2])")
    rho = corspearman(score_array[:,1],score_array[:,2])
    p = pvalue(OneSampleZTest(atanh(rho), 1, length(all_genes)))
    println("Spearman Corr: $(round(rho,digits=4)), P = $(p)")
    println("$(ORDER[1]) vs $(ORDER[3])")
    rho = corspearman(score_array[:,1],score_array[:,3])
    p = pvalue(OneSampleZTest(atanh(rho), 1, length(all_genes)))
    println("Spearman Corr: $(round(rho,digits=4)), P = $(p)")
    for dir in 1:length(SCORE_DIRS)
        println("$(ORDER[1]) vs $(ORDER[4]), $(SCORE_DIRS[dir]) window:")
        rho = corspearman(score_array[:,1],score_array[:,NUM_OTHER_SCORES+dir])
        p = pvalue(OneSampleZTest(atanh(rho), 1, length(all_genes)))
        println("Spearman Corr: $(round(rho,digits=4)), P = $(p)")

        println("\n$(ORDER[1]) vs $(ORDER[5]), $(SCORE_DIRS[dir]) window:")
        rho = corspearman(score_array[:,1],score_array[:,NUM_OTHER_SCORES+dir+length(SCORE_DIRS)])
        p = pvalue(OneSampleZTest(atanh(rho), 1, length(all_genes)))
        println("Spearman Corr: $(round(rho,digits=4)), P = $(p)")

        println("\n$(ORDER[4]) vs $(ORDER[5]), $(SCORE_DIRS[dir]) window:")
        rho = corspearman(score_array[:,NUM_OTHER_SCORES+dir],score_array[:,NUM_OTHER_SCORES+dir+length(SCORE_DIRS)])
        p = pvalue(OneSampleZTest(atanh(rho), 1, length(all_genes)))
        println("Spearman Corr: $(round(rho,digits=4)), P = $(p)")

        sorted_genes = all_genes[sortperm(score_array[:,NUM_OTHER_SCORES+dir],rev=true)]
        sorted_genes = all_genes[sortperm(score_array[:,NUM_OTHER_SCORES+dir+length(SCORE_DIRS)],rev=true)]
    end
end

function main()
    scores = readQx()
    scores = readiHSnSL(scores)
    scores = readLength(scores)
    scores = readOELoF(scores)
    scores = readPhyloP(scores)
    scores = readModels(scores)

    println(scores["ENSG00000108510"])
    println(scores["ENSG00000165494"])

    # printCorr(scores)

    #assemble into array, remove genes that are missing a score
    all_genes = collect(keys(scores))
    keep_inds = Array{Int64,1}[]
    score_array = repeat([0.0],length(all_genes),NUM_OTHER_SCORES+2*length(SCORE_DIRS)) #gene x score
    for gene in 1:length(all_genes)
        if !in(-Inf,scores[all_genes[gene]])
            keep_inds = vcat(keep_inds,[gene])
            score_array[gene,:] = scores[all_genes[gene]]
        end
    end
    println("Num of genes total: $(length(all_genes))")
    all_genes = all_genes[keep_inds]
    score_array = score_array[keep_inds,:]
    println("Num of genes with all scores: $(length(all_genes))")

    s_plot = scatterplot(score_array[:,findfirst(x->x=="iHS",ORDER)],score_array[:,findfirst(x->x=="Best_Qx_P",ORDER)],
            color=:black,alpha=0.5)
    s_plot.set_ylabel("Best Qx P Value")
    s_plot.set_xlabel("iHS")
    Seaborn.savefig("ihs_qx_scatter.pdf")
    clf()

    s_plot = scatterplot(score_array[:,findfirst(x->x=="gene_length",ORDER)],score_array[:,findfirst(x->x=="Best_Qx_P",ORDER)],
            color=:black,alpha=0.5)
    s_plot.set_ylabel("Best Qx P Value")
    s_plot.set_xlabel("Gene length")
    Seaborn.savefig("gene_length_qx_scatter.pdf")
    clf()

    corrs = corspearman(score_array)
    h = Seaborn.heatmap(corrs,xticklabels=ORDER,yticklabels=ORDER,square=true,cmap=COL_MAP,vmin=-1,vmax=1,center=0,annot=true,annot_kws=Dict("size"=> 10))
    h.tick_params(axis="x",labelsize="8")
    h.tick_params(axis="y",labelsize="8")
    Seaborn.savefig("score_corr_heatmap.pdf")
    clf()
end

main()
