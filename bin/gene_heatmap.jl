# gene_heatmap.jl
# @author Laura colbran
#
# julia1.5.3

using ArgParse,GZip
using Statistics
using Seaborn


Seaborn.set(style="white", palette="muted")
COL_MAP = "coolwarm"
GROUP_SIZES  = [22,14,5,1,7]
COLOURS = ["#b93e0f","#a1c6d8","#5f3912","#feae02","#6b878f","#c6d737","#171b17","#fe8515","#b7c1c2","#9a5611"]

# ENV["GKSwstype"] = "100" #fixes a segfault bug in GR backend for Plots.jl: https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--qx","-q"
            help = "path to tissue files for Qx"
            arg_type = String
            nargs = '*'
        "--genes","-g"
            arg_type = String
            help = "Ensembl ID for gene to plot"
        "--populations","-p"
            arg_type = String
            help = "File path for sample -> population assignments."
        "--tiss_order","-t"
            help = "list of tissue names in order you want them plotted. Also used to filter."
            arg_type = String
            default = ["JTI_Adipose_Subcutaneous","JTI_Adipose_Visceral_Omentum",
                    "JTI_Artery_Aorta","JTI_Artery_Coronary","JTI_Artery_Tibial",
                    "JTI_Breast_Mammary_Tissue","JTI_Colon_Sigmoid","JTI_Colon_Transverse",
                    "JTI_Esophagus_Gastroesophageal_Junction","JTI_Esophagus_Muscularis",
                    "JTI_Heart_Atrial_Appendage","JTI_Heart_Left_Ventricle","JTI_Liver",
                    "JTI_Lung","JTI_Muscle_Skeletal","JTI_Nerve_Tibial","JTI_Ovary",
                    "JTI_Pancreas","JTI_Prostate","JTI_Stomach","JTI_Uterus","JTI_Vagina",
                    "JTI_Brain_Amygdala",
                    "JTI_Brain_Anterior_cingulate_cortex_BA24","JTI_Brain_Caudate_basal_ganglia",
                    "JTI_Brain_Cerebellar_Hemisphere","JTI_Brain_Cerebellum","JTI_Brain_Cortex",
                    "JTI_Brain_Frontal_Cortex_BA9","JTI_Brain_Hippocampus","JTI_Brain_Hypothalamus",
                    "JTI_Brain_Nucleus_accumbens_basal_ganglia","JTI_Brain_Putamen_basal_ganglia",
                    "JTI_Brain_Spinal_cord_cervical_c-1","JTI_Brain_Substantia_nigra","JTI_Pituitary",
                    "JTI_Esophagus_Mucosa","JTI_Skin_Not_Sun_Exposed_Suprapubic",
                    "JTI_Skin_Sun_Exposed_Lower_leg","JTI_Thyroid","JTI_Cells_Cultured_fibroblasts",
                    "JTI_Testis",
                    "JTI_Adrenal_Gland","JTI_Cells_EBV-transformed_lymphocytes","JTI_Kidney_Cortex",
                    "JTI_Minor_Salivary_Gland","JTI_Small_Intestine_Terminal_Ileum","JTI_Spleen",
                    "JTI_Whole_Blood"]
            nargs='*'
        "--pop_order","-o"
            help = "list of tissue names in order you want them plotted. Also used to filter."
            arg_type = String
            default = ["YRI","ASW","MSL","ESN","LWK","ACB","GWD","CEU","FIN","TSI","IBS","GBR",
                        "CHB","JPT","CHS","CDX","KHV","GIH","PJL","BEB","ITU","STU",
                        "PUR","PEL","MXL","CLM"]
            nargs='*'
        "--pop_colors","-l"
            help = "list of matplotlib colors for the population axis labels."
            arg_type = String
            nargs='*'
            default = vcat(repeat(["#946bb3"],7),vcat(repeat(["#4f7ac9"],5),
                vcat(repeat(["#d48c63"],5),vcat(repeat(["#c86f6f"],5),repeat(["#74c46e"],4)))))
        "--tiss_colors","-r"
            help = "list of matplotlib colors for the tissue axis labels."
            arg_type = String
            nargs='*'
            default = vcat(repeat([COLOURS[1]],GROUP_SIZES[1]),vcat(repeat([COLOURS[2]],GROUP_SIZES[2]),
                vcat(repeat([COLOURS[3]],GROUP_SIZES[3]),vcat(repeat([COLOURS[4]],GROUP_SIZES[4]),repeat([COLOURS[5]],GROUP_SIZES[5])))))
        "--scale","-s"
            help = "to scale by expression amount"
            action=:store_true
        "--group"
            help = "if you want to summarize heatmap by groups of tissues"
            action=:store_true
        "--expression","-e"
            help = "only if scaling. path to expression matrix to scale by. Assume GCT file"
            default = "."
            arg_type = String
    end
    return parse_args(s)
end

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

# scale alpha by expression
function calcAlphas(gene_id::String,tissues::Array{String,1},expr_path::String,mask_inds::Array{Int64,1})
    # println(mask_inds)
    meds = Float64[]
    tiss_cols = String[]
    # read median expression
    GZip.open(expr_path) do inf
        for line in eachline(inf)
            if startswith(line,"Name")
                tiss_cols = split(chomp(line),'\t')[3:end] # pull tissues to match indices
            elseif startswith(line,gene_id)
                meds = parse.(Float64,split(chomp(line),'\t')[3:end]) # pull line for medians
                break
            end
        end
    end
    # convert tissue names to match what the heatmap uses
    for i in 1:length(tiss_cols)
        tiss_cols[i] = replace(tiss_cols[i]," - "=>"_")
        tiss_cols[i] = replace(tiss_cols[i]," "=>"_")
        tiss_cols[i] = replace(tiss_cols[i],"("=>"")
        tiss_cols[i] = replace(tiss_cols[i],")"=>"")
        tiss_cols[i] = "JTI_$(tiss_cols[i])"
    end
    #for each tissue, calculate the alpha
    alphas = repeat([1.0],length(tissues))
    for (i,tiss) in enumerate(tissues)
        alphas[i] = meds[findfirst(x->x==tiss,tiss_cols)]
    end
    alphas[mask_inds] .= -100.00
    max_med = maximum(alphas)
    alphas = round.(alphas./max_med;digits=2)
    return alphas
end

function geneHeatMap(qx_paths::Array{String,1},samp_path::String,tiss_order::Array{String,1},pop_order::Array{String,1},
        gene_id::String,pop_label_colors::Array{String,1},tiss_label_colors::Array{String,1},scale::Bool,expr_file::String,
        group::Bool)
    groups = sampDict(samp_path,pop_order) # group -> sample IDs
    samples = String[]
    tissDict = Dict{String,Array{Float64,1}}() #to hold predictions
    for t in tiss_order
        tissDict[t] = [-100.0]
    end
    indices = Dict{String,Array{Int64,1}}() #group -> indices of its samples
    for pop_path in qx_paths
        tiss = join(split(splitdir(pop_path)[2],"_")[1:end-2],"_")
        println(tiss)
        GZip.open("$(pop_path)") do f
            for line in eachline(f)
                l = split(chomp(line),'\t')
                if startswith(line,"gene")
                    samples = l[2:end]
                    indices = pullInds(groups,l[2:end])
                    continue
                end
                if l[1] != gene_id continue end
                tissDict[tiss] = [parse(Float64,i) for i in l[2:end]]
            end
        end
    end
    # assemble matrix of tiss x pop in correct order
    mat = zeros(Float64,length(tiss_order),length(pop_order))
    for i in 1:length(tiss_order)
        if tissDict[tiss_order[i]] == [-100.0]
            mat[i,1:length(pop_order)] .= -100.0
        else
            for j in 1:length(pop_order)
                mat[i,j] = median(tissDict[tiss_order[i]][indices[pop_order[j]]])
            end
        end
    end
    alphas = repeat([1.0],length(tiss_order))
    if scale
        alphas = calcAlphas(gene_id,tiss_order,expr_file,findall(x->x==-100,mat[:,1])) #sorted by tiss_order
    end
    # mat = [0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         1.0 1 1 1 1.8 1.7 1.6 1 1 1 1 -2.0 1 1.1 1.4 -1.1 1.5 1.7 1.9 1 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 0 0 0 0.8 0.7 0.6 0 0 0 0 -2.0 0 0.1 0.4 -0.1 0.5 0.7 0.9 0 2.1 3.2 0.1 2.2 1.6 2.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3;
    #         0.0 3.2 1.2 2.2 0.8 0.7 0.6 -0.6 -0.5 -1.0 -1.6 -2.0 -1.2 0.1 0.4 -0.1 0.5 0.7 0.9 -0.9 2.1 3.2 0.1 0.2 0.6 0.3]
    # alphas = [1.0,0.5,0.1,0.1,0.2,0.4,0.3,0.2,0.1,0.7,0.6,0.5,0.2,0.1,0.5,0.1,0.8,0.4,0.2,0.6,0.2,
    #                 0.5,0.8,0.5,0.6,0.2,0.8,0.7,0.6,0.5,0.2,0.3,0.4,0.2,0.5,0.6,0.7,0.4,0.3,0.2,0.1,
    #                 0.7,1.0,1.0,0.5,0.3,0.2,0.1,0.4]
    # gene_id = "test"
    if group
        new_mat = zeros(Float64,length(GROUP_SIZES),length(pop_order))
        new_alphas = repeat([1.0],length(GROUP_SIZES))
        new_cols = repeat([COLOURS[1]],length(GROUP_SIZES))
        tiss_order = ["Group$(i)" for i in 1:length(GROUP_SIZES)]
        start_ind = 1
        for grp in 1:length(GROUP_SIZES)
            tmp = mat[start_ind:start_ind+GROUP_SIZES[grp]-1,:]
            new_mat[grp,:] =  [median(tmp[findall(x->x!=-100.0,tmp[:,pop]),pop]) for pop in 1:length(pop_order)] #avg prediction
            new_cols[grp] = tiss_label_colors[start_ind]
            if scale
                new_alphas[grp] = maximum(alphas[start_ind:start_ind+GROUP_SIZES[grp]-1]) #pick largest alpha for display
            end
            start_ind = start_ind+GROUP_SIZES[grp]
        end
        mat = new_mat
        alphas = new_alphas
        tiss_label_colors = new_cols
    end
    mask_inds = repeat([true],length(tiss_order),length(pop_order))
    # parameters for heatmap
    color_max = maximum(mat)
    color_min = minimum(mat[findall(x->x!=-100.0,mat)])
    cent = color_min + (color_max-color_min)/2.0
    # first, plot an empty heatmap with the right axis labels, cbar, etc.
    h = Seaborn.heatmap(mat,xticklabels=pop_order,yticklabels=tiss_order,square=true,cmap=COL_MAP,mask=mask_inds,vmin=color_min,vmax=color_max,center=cent)
    [t.set_color(i) for (i,t) in zip(pop_label_colors,h.xaxis.get_ticklabels())]
    [t.set_color(i) for (i,t) in zip(tiss_label_colors,h.yaxis.get_ticklabels())]
    h.tick_params(axis="x",labelsize="4")
    h.tick_params(axis="y",labelsize="4")
    ylabel("Tissue")
    xlabel("Population")
    # plot heatmap for each alpha
    for alph in unique(alphas)
        if alph < 0 continue end #skip alphas on masked tissues
        mask_inds = repeat([false],length(tiss_order),length(pop_order)) # clear previous mask
        mask_inds[findall(x->x==-100.0,mat)] .= true #mask missing tissues

        mask_inds[findall(x->x!=alph,alphas),:] .= true #mask tissues with a different alpha
        h = Seaborn.heatmap(mat,xticklabels=pop_order,yticklabels=tiss_order,square=true,cmap=COL_MAP,alpha=alph,
            mask=mask_inds,vmin=color_min,vmax=color_max,center=cent,cbar = false)
    end
    Seaborn.savefig("$(gene_id)_qx_heatmap.pdf")
    clf()
end

function main()
    parsed_args = parseCommandLine()
    geneHeatMap(parsed_args["qx"],parsed_args["populations"],parsed_args["tiss_order"],parsed_args["pop_order"],parsed_args["genes"],parsed_args["pop_colors"],parsed_args["tiss_colors"],parsed_args["scale"],parsed_args["expression"],parsed_args["group"])
end

main()
