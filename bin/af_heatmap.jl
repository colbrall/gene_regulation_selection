using Seaborn,DataFrames
using StatsBase:mean
COL_MAP = "coolwarm"
NUM_POPS = 26
# NUM_SNPS = 12
NUM_SNPS = 4
# NUM_SNPS = 14
# NUM_SNPS = 9
# NUM_SNPS = 10
# NUM_SNPS = 13
# GENE = "samd10"
GENE = "fads1"
# GENE = "aco2"
# GENE = "itgam"
# GENE = "sult1a2"
# GENE = "khk"

open("$(GENE)_info.txt") do inf

Seaborn.set(style="white", palette="muted")

pop_order= ["YRI","ASW","MSL","ESN","LWK","ACB","GWD","CEU","FIN","TSI","IBS","GBR",
                "GIH","PJL","BEB","ITU","STU","CHB","JPT","CHS","CDX","KHV",
                "PUR","PEL","MXL","CLM"]
pop_label_colors= vcat(repeat(["#946bb3"],7),vcat(repeat(["#4f7ac9"],5),
    vcat(repeat(["#c86f6f"],5),vcat(repeat(["#d48c63"],5),repeat(["#74c46e"],4)))))

header_pops = ["ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH",
"GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI"]

df = DataFrames.DataFrame(pop=String[],snp=String[],effect=Float64[],af=Float64[])
# for line in eachline(inf)
#     pop_count = 0
#     l = split(chomp(line),"\t")
#     for pop in l[5:end]
#         pop_count +=1
#         afs= split(pop,",")
#         af = afs[findall(x->startswith(x,l[4]),afs)]
#         push!(df,[header_pops[pop_count],l[1],parse(Float64,l[2]),parse(Float64,split(af[1],":")[2])])
#     end
# end
#
# println(first(df,6))

# new_order = zeros(Int64,nrow(df))
# tot_snps = length(unique(df[!,:snp]))
# for index in 1:nrow(df)
#     num_snp = findfirst(x->x==df[index,:snp],unique(df[!,:snp]))
#     pop_place = findfirst(x->x==df[index,:pop],pop_order)
#     new_order[tot_snps*(pop_place-1)+num_snp] = index
# end
# df = df[new_order,:]
# df = df[df[!,:af] .!= 0,:] #don't plot boxes where the SNP is absent
# max_eff = maximum(abs.(df[!,:effect]))
#
# r = Seaborn.relplot(x=df[!,:pop],y=df[!,:snp],hue=df[!,:effect],size=df[!,:af],palette=COL_MAP,
#             aspect=2.25,hue_norm=(-max_eff,max_eff), legend=true,
#             sizes = (40,390),size_norm = (minimum(df[!,:af]),maximum(df[!,:af])))
# [t.set_color(i) for (i,t) in zip(pop_label_colors,r.ax.get_xticklabels())]
# r.ax.set_xticklabels(r.ax.get_xticklabels(),rotation = 90)
#
# ylabel("rsID")
# xlabel("Population")
# ghost = r.ax.scatter([], [], c=[], vmin=-max_eff, vmax=max_eff, cmap=COL_MAP)
# r.fig.colorbar(ghost)
# Seaborn.savefig("$(GENE)_snp_heatmap.pdf")
# clf()

rsids = String[]
mafs = zeros(Float64,NUM_SNPS,NUM_POPS)
row = 1
for line in eachline(inf)
    l = split(chomp(line),"\t")
    rsids = vcat(rsids,l[1])
    # eff_dir = vcat(eff_dir,sign(parse(Float64,l[2])))
    col = 1
    for pop in l[5:end]
        afs= split(pop,",")
        # println(afs)
        af = afs[findall(x->startswith(x,l[4]),afs)]
        mafs[row,col] = parse(Float64,split(af[1],":")[2]) .* parse(Float64,l[2])
        col +=1
    end
    row+=1
end

mafs=mafs[:,[findfirst(x->x==pop,header_pops) for pop in pop_order]]
#set mean to zero for each SNP
for row in 1:size(mafs,1)
    mafs[row,:] = mafs[row,:] .- mean(mafs[row,:])
end

# snp_label_colors = [ifelse(eff == 1.0,"#b80000","#0106a5") for eff in eff_dir]

h = Seaborn.heatmap(mafs,xticklabels=pop_order,yticklabels=rsids,square=true,cmap=COL_MAP,center=0)
[t.set_color(i) for (i,t) in zip(pop_label_colors,h.xaxis.get_ticklabels())]
# [t.set_color(i) for (i,t) in zip(snp_label_colors,h.yaxis.get_ticklabels())]
h.tick_params(axis="x",labelsize="8")
h.tick_params(axis="y",labelsize="8")
ylabel("rsID")
xlabel("Population")
Seaborn.savefig("$(GENE)_snpProduct_heatmap.pdf")
clf()

end
