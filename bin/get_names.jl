# get_names.jl
# @author Laura Colbran
#
# assumes a bed file with col 4 as gene name and col 5 as ensembl ID
#
# prints everything from input with additional name column
# julia 1.5

using ArgParse

function parseCommandLine()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--ref_file","-r"
            arg_type = String
            help = "path to file with gene names and IDs"

        "--in_file","-i"
            arg_type = String
            help = "path to file you want to convert"

        "--column","-c"
            arg_type = Int64
            help = "column in convert file with the ID"
    end
    return parse_args(s)
end

# reads reference file into a dictionary
function readDict(path::String)
    dict = Dict{String,String}()
    open(path) do f
        for line in eachline(f)
          e = split(chomp(line),"\t")
          dict[split(e[5],".")[1]] = e[4]
        end
    end
    return dict
end

#id -> name
function idToName(ref_path::String,quest_path::String,col::Int64)
    ref_dict = readDict(ref_path)
    open(quest_path) do f
        for line in eachline(f)
            if startswith(line, "gene")
                println("$(chomp(line))\tName")
                continue
            end
            l = split(chomp(line),"\t")
            println("$(join(l,"\t"))\t$(ref_dict[l[col]])")
        end
    end
end

function main()
    parsed_args = parseCommandLine()
    idToName(parsed_args["ref_file"],parsed_args["in_file"],parsed_args["column"])
end

main()
