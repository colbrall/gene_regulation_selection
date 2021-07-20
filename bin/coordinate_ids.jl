# coordinate_ids.jl
#
# @author Laura Colbran
# maps rsIDs in PrediXcan dbs to coordinate IDs and writes a map file
# depends on sqlite3

using GZip
using ArgParse
# using SQLite

SNP_COL = "rsid" #saving variable model DB uses
# GENE_COL = :gene #saving variable model DB uses
# REF_COL = :ref_allele
# EFF_COL = :eff_allele
# WEIGHT_COL = :weight #saving variable model DB uses

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--db","-d"
            help = "path to sqlite database(s)"
            nargs = '*'
            arg_type = String
            required=true
        "--snp_file","-s"
            help = "path to dbSNP file to parse"
            default = "/project/mathilab/colbranl/gene_regulation_selection/data/snp150.txt.gz"
        "--out_dir","-o"
            help = "directory to write ID map file to"
            default = "./"
    end
    return parse_args(s)
end

# search db for each SNP, return rsID, ref, alt
function parseDB(path::String,rsid::String,chr::String,pos::String)
    genes = ""
    # open SQLite connection, read weights table
    q = "SELECT * FROM 'weights' WHERE $SNP_COL=='$rsid'"
    # println("$(splitdir(path)[2]): $q")
    command = `sqlite3 $(path) "$(q)"`
    open(command) do inf
        for snp in eachline(inf)
        # for snp in DBInterface.execute(SQLite.DB(path),q)
            l = split(snp,"|")
            genes = "$(genes)\n$(splitdir(path)[2])\t$(l[2])\t$(rsid)\t$(chr)_$(pos)_$(l[4])_$(l[5])_b38"
            # genes = "$(genes)\n$(splitdir(path)[2])\t$(snp[GENE_COL])\t$(rsid)\t$(chr)_$(pos)_$(snp[REF_COL])_$(snp[EFF_COL])_b38"
        end
    end
    return genes
end

function coordinateID(db_paths::Array{String,1},snp_path::String,out_dir::String)
    GZip.open(snp_path) do inf
        open("$(out_dir)mapped_snps.txt","w") do outf
            write(outf,"#tiss\tgene\trsID\tcoordID")
            for line in eachline(inf)
                l = split(chomp(line),"\t")
                for db in db_paths
                    # println(db)
                    out_lines = parseDB(db,"$(l[5])","$(l[2])","$(l[4])")
                    write(outf,out_lines)
                end
            end
            write(outf,"\n")
        end
    end
end


function main()
    parsed_args = parseCommandLine()
    coordinateID(parsed_args["db"],parsed_args["snp_file"],parsed_args["out_dir"])
end

main()
