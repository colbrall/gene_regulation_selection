# qx_per_gene.jl
#
# @author Laura Colbran 2021-07-16
# given database of models, pulls matching SNPs by GTEx fqcy for use in Qx calculation
# outputs file with gene and matching SNP ids
# julia 1.5

using ArgParse
using SQLite
using StatsBase
using GZip

BIN_COL = 3 #column in bin file that contains bin id
SNP_COL = :rsid #saving variable model DB uses
GENE_COL = :gene #saving variable model DB uses
DBSNP_FILE = "/project/mathilab/colbranl/gene_regulation_selection/data/jti_snp_coordinates.txt.gz"

function parseCommandLine()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--db","-d"
            help = "sqlite database with genes and snp effect sizes"
            arg_type = String
        "--snps","-s"
            help = "path to file with ascertainment SNP AFs and any other bin info. assumes snp ids are in form chr_pos_ref_alt_b38"
            arg_type = String
            default = ""
        "--num_match","-n"
            help = "number of neutral SNPs to pull to match each test SNP. Default 100"
            arg_type = Int64
            default = 100
    end
    return parse_args(s)
end

# read in database, then for each gene pull snps and effect sizes
function parseDB(path::String)
    genes = Dict{String,Array{String,1}}() # Dict{gene => [(id,weight)]}
    # open SQLite connection, read weights table
    q = "SELECT * FROM 'weights'"
    for snp in DBInterface.execute(SQLite.DB(path),q)
        if !(snp[GENE_COL] in keys(genes))
            genes[snp[GENE_COL]] = [snp[SNP_COL]]
        else
            genes[snp[GENE_COL]] = vcat(genes[snp[GENE_COL]],[snp[SNP_COL]])
        end
    end
    return genes
end

function mapIDs(snp_path)
    ids = Dict{String,Array{String,1}}()
    GZip.open(snp_path) do f
        for line in eachline(f)
            if startswith(line,"#") continue end
            l = split(chomp(line),"\t")
            if occursin("_",l[1]) continue end #skip non-standard chromosomes
            ids[l[4]] = [join([l[1],l[3],l[5],x,"b38"],"_") for x in split(l[7],",")[1:end-1]]
        end
    end
    return ids
end

function coordinateID(ids::Dict{String,Array{String,1}},rsids::Array{String,1})
    coord_ids = String[]
    for rsid in rsids
        try
            coord_ids = vcat(coord_ids,ids[rsid])
        catch e
            coord_ids = vcat(coord_ids,[""])
        end
    end
    return coord_ids
end

# identifies matched sets of SNPs based on bins
function matchSNPs(db_path::String,bin_path::String,num_to_match::Int64)
    bin_snps = Dict{String,Array{String,1}}()
    coord_ids = mapIDs(DBSNP_FILE)
    genes = parseDB(db_path)

    open(bin_path) do inf
        # search for and save all SNPs for each bin
        # N = 1
        for line in eachline(inf)
            # if N%300000 == 0
            #     println(line)
            #     break
            # end
            id = split(chomp(line),"\t")[1]
            bin = split(chomp(line),"\t")[BIN_COL]
            if in(bin,keys(bin_snps))
                bin_snps[bin] = vcat(bin_snps[bin],"$id")
            else
                bin_snps[bin] = ["$id"]
            end
            # N+=1
        end
    end

    outdir = "$(splitext(db_path)[1])/"
    if !isdir(outdir) mkdir(outdir) end
    out_path= "$(outdir)/matched_snps.txt"
    open(out_path,"w") do outf
        for gene in keys(genes)
            gene_snps = genes[gene]
            if !occursin("_",gene_snps[1]) #if IDs are not in coordinate ID form
                gene_snps = coordinateID(coord_ids,gene_snps)
            end
            # if !occursin("chr1_",gene_snps[1]) continue end
            # println(gene_snps)
            gene_string = String[]
            for bin in keys(bin_snps)
                targ_snps = intersect(gene_snps,bin_snps[bin])
                # println(targ_snps)
                options = bin_snps[bin][findall(x->!in(x,targ_snps),bin_snps[bin])] #ensure we don't include the target snps
                matches = sample(options,num_to_match*length(targ_snps),replace=false)
                for snp in matches
                    gene_string = vcat(gene_string,"$(gene)\t$snp")
                end
            end
            outs = join(gene_string,"\n")
            write(outf,"$(outs)\n")
            # exit()
        end
    end
end

function main()
    parsed_args = parseCommandLine()
    matchSNPs(parsed_args["db"],parsed_args["snps"],parsed_args["num_match"])
end

main()
