#filter_dbs.jl
#
# given list of gene-tissue pairs, puts all modle information into one DB

using SQLite, DataFrames, DBInterface
using ArgParse


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
        "--out_db","-o"
            arg_type = String
            help = "Name for the new DB"
            default = "filter.db"
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

function main()
    parsed_args = parseCommandLine()
    models_to_pull = readModels(parsed_args["models"]) #tiss->genes dictionary

    new_db = SQLite.DB(parsed_args["out_db"])
    SQLite.createtable!(new_db, "weights", Tables.Schema(("rsid", "gene", "weight", "ref_allele", "eff_allele"),
                        (String, String, Float64, String, String)))
    weight_q = "INSERT INTO weights (rsid,gene, weight, ref_allele, eff_allele) VALUES (?, ?, ?, ?, ?)"
    # SQLite.createtable!(new_db, "extra", Tables.Schema(("gene", "genename", "geneannot", "tissue", "R2", "snps", "pval", "qval"),
    #                     (String, String, String, String, Float64, Int, Float64, Float64)))
    SQLite.createtable!(new_db, "extra", Tables.Schema(("gene", "genename", "tissue", "R2", "snps", "pval", "qval"),
                        (String, String, String, Float64, Int, Float64, Float64)))
    # extra_q = "INSERT INTO extra (gene, genename, geneannot, tissue, R2, snps, pval, qval) VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
    extra_q = "INSERT INTO extra (gene, genename, tissue, R2, snps, pval, qval) VALUES (?, ?, ?, ?, ?, ?, ?)"
    gene_q = "SELECT * FROM 'extra'"

    for tiss in keys(models_to_pull)
        db = SQLite.DB("$(parsed_args["databases"])$(tiss).db")
        for gene in DBInterface.execute(db,gene_q)
            if !in(gene[:gene], models_to_pull[tiss]) continue end
            stmt = DBInterface.prepare(new_db, extra_q)
            # row = (gene=gene[:gene], genename=gene[:genename], geneannot=gene[:geneannot],tissue=tiss, R2=gene[Symbol("pred.perf.R2")], snps=gene[Symbol("n.snps.in.model")], pval=gene[Symbol("pred.perf.pval")], qval=gene[Symbol("pred.perf.qval")])
            row = (gene=gene[:gene], genename=gene[:genename], tissue=tiss, R2=gene[Symbol("pred.perf.R2")], snps=gene[Symbol("n.snps.in.model")], pval=gene[Symbol("pred.perf.pval")], qval=gene[Symbol("pred.perf.qval")])
            # DBInterface.execute(stmt, (row.gene, row.genename, row.geneannot, row.tissue, row.R2, row.snps, row.pval, row.qval))
            DBInterface.execute(stmt, (row.gene, row.genename, row.tissue, row.R2, row.snps, row.pval, row.qval))
            snp_q = "SELECT * FROM 'weights' WHERE gene='$(gene[:gene])'"
            for snp in DBInterface.execute(db,snp_q)
                stmt = DBInterface.prepare(new_db, weight_q)
                row = (rsid=snp[:rsid],gene=snp[:gene], weight=snp[:weight], ref_allele=snp[:ref_allele], eff_allele=snp[:eff_allele])
                DBInterface.execute(stmt, (row.rsid, row.gene, row.weight, row.ref_allele, row.eff_allele))
            end
        end
    end



end

main()
