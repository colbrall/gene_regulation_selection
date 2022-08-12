# best_models.jl
# for every gene, searches databases, and outputs the tissue with the best models
# you can also give it a gene x tissue matrix to filter to just these best models
using ArgParse,SQLite,DataFrames

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--db_dir","-d"
            help = "path to directory with all DBs to search"
            arg_type = String
        "--matrix","-m"
            help = "path to gene x tissue matrix that you want to filter down to"
            arg_type = String
            default = ""
        "--output","-o"
            help = "path to save filtered matrix to"
            arg_type = String
            default = "filtered.txt"
        "--write", "-w"
            help = "if you want to save best models to file"
            action=:store_true
    end
    return parse_args(s)
end

function pullBestModels(db_dir::String,wr::Bool)
    mod_dict = Dict{String,Tuple{String,Int64,Float64}}()
    for item in readdir(db_dir)
        if !endswith(item,"db") continue end
        db_name = join(split(basename(item),".")[1:(end-1)],".")
        q = "SELECT * FROM 'extra'"
        tmp = DataFrame(DBInterface.execute(SQLite.DB("$(db_dir)$(item)"),q))
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
    #save list of best models
    if wr
        println("Saving Best Models...")
        open("best_models.txt","w") do outf
            write(outf,"gene\ttissue\tnumSNPs\tR2\n")
            for gene in keys(mod_dict)
                write(outf,"$(gene)\t$(join(mod_dict[gene],'\t'))\n")
            end
        end
    end
    return mod_dict
end

function filterBest(db_dir::String,mat_path::String,wr::Bool,output::String)
    mods = pullBestModels(db_dir,wr)
    if mat_path == "" exit() end #if you didn't give a matrix to filter, end here
    open(output,"w") do outf
        open(mat_path) do inf
            header = String[]
            for line in eachline(inf)
                if startswith(line,"gene")
                    header = split(chomp(line),"\t")
                    continue
                end
                l = split(chomp(line),"\t")
                if !haskey(mods,l[1]) continue end
                best_tiss = mods[l[1]][1]
                write(outf,"$(l[1])\t$(best_tiss)\t$(l[[t==best_tiss for t in header]][1])\n")
            end
        end
    end
end

function main()
    parsed_args = parseCommandLine()
    filterBest(parsed_args["db_dir"],parsed_args["matrix"],parsed_args["write"],parsed_args["output"])
end

main()
