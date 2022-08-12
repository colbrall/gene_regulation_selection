# join_selection_output.jl
#
# joins raw output from selscan and converts to bed format for intersecting

using ArgParse

function parseCommandLine()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tables","-t"
            help = "paths to results tables to join"
            arg_type = String
            nargs='*'
        "--column","-c"
            arg_type = Int64
            help = "column of stat to keep"
        "--output","-o"
            arg_type = String
            help = "path of file to save to"
    end
    return parse_args(s)
end

function joinTables(in_files::Array{String,1},col::Int64,out_path::String)
    pops = String["rsID"]
    chr = split(split(basename(in_files[1]),".")[1],"_")[2]
    snps = Dict{String,Array{String,1}}()
    num_pops = length(in_files)
    iter = 2
    for file in in_files
        pop = split(basename(file),"_")[1]
        pops = vcat(pops,"$(pop)")
        open(file) do inf
            for line in eachline(inf)
                l = split(chomp(line),"\t")
                if !haskey(snps,l[2])
                    snps[l[2]] = vcat(l[1],repeat(["NA"],num_pops))
                end
                snps[l[2]][iter] = "$(l[col])"
            end
        end
        iter += 1
    end
    open("header$(chr)","w") do outf
        write(outf,"#chr\tpos-1\tpos\t$(join(pops,"\t"))\n")
    end
    open("tmp$(chr).txt","w") do outf
        for snp in keys(snps)
            pos = parse(Int64,snp)
            write(outf,"$(chr)\t$(pos-1)\t$(pos)\t$(join(snps[snp],"\t"))\n")
        end
    end
    cmda = `sortBed -i tmp$(chr).txt`
    cmdb = `cat header$(chr) -`
    cmdc = `bgzip`
    run(pipeline(pipeline(pipeline(cmda,cmdb),cmdc),stdout=out_path))
    run(`rm header$(chr)`)
    run(`rm tmp$(chr).txt`)
end

function main()
  parsed_args = parseCommandLine()
  joinTables(parsed_args["tables"],parsed_args["column"],parsed_args["output"])
end

main()
