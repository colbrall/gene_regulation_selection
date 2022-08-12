using Seaborn
open("samd10_info.txt") do inf

Seaborn.set(style="white", palette="muted")

pop_order= ["YRI","ASW","MSL","ESN","LWK","ACB","GWD","CEU","FIN","TSI","IBS","GBR",
                "CHB","JPT","CHS","CDX","KHV","GIH","PJL","BEB","ITU","STU",
                "PUR","PEL","MXL","CLM"]
pop_label_colors= vcat(repeat(["#946bb3"],7),vcat(repeat(["#4f7ac9"],5),
    vcat(repeat(["#d48c63"],5),vcat(repeat(["#c86f6f"],5),repeat(["#74c46e"],4)))))

header_pops = ["ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH",
"GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI"]

COL_MAP = "Greys"
mafs = zeros(Float64,12,26)
rsids = String[]
row = 1
for line in eachline(inf)
    l = split(chomp(line)," ")
    rsids = vcat(rsids,l[1])
    col = 1
    for pop in l[16:end]
        afs= split(pop,",")
        af = afs[findall(x->startswith(x,l[5]),afs)]
        mafs[row,col] = parse(Float64,split(af[1],":")[2])
        col +=1
    end
    row+=1
end
display(mafs)

mafs=mafs[:,[findfirst(x->x==pop,header_pops) for pop in pop_order]]

h = Seaborn.heatmap(mafs,xticklabels=pop_order,yticklabels=rsids,square=true,cmap=COL_MAP)
[t.set_color(i) for (i,t) in zip(pop_label_colors,h.xaxis.get_ticklabels())]
h.tick_params(axis="x",labelsize="8")
h.tick_params(axis="y",labelsize="8")
ylabel("rsID")
xlabel("Population")
Seaborn.savefig("samd10_snps.pdf")
clf()

end
