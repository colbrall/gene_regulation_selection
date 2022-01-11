# predixcan_search.jl
# given ensembl id(s), searches predixcan model database(s) for a model. returns SNP count and r2
# @author Laura Colbran
# julia 1.4

using SQLite
using ArgParse
using DataFrames
using DBInterface

# parses command-line arguments
function parseCommandLine()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--databases", "-d"
			nargs='*'
			help="path to model database(s)"
			arg_type = String
			required = true
		"--genes","-g"
			help="file with gene IDs to search for"
			arg_type=String
			default ="."
		"--snps","-s"
			help="file(s) with gene IDs to search for"
			arg_type=String
			nargs='*'
		"--column","-c"
			help = "1-indexed column containing gene IDs or rsIDs to search for"
			arg_type=Int64
	end
	return parse_args(s)
end

function geneSearch(dbs::Array{String,1},gene_file::String,col::Int64)
	println("gene_id\tdb\tnum_SNPs\tr2")
	open(gene_file) do f
		for line in eachline(f)
			if startswith(line, "#") continue end
			l = split(chomp(line),"\t")
			q = "SELECT * FROM 'extra' WHERE gene='$(l[col])'"
			for db in dbs
				db_name = join(split(basename(db),".")[1:(end-1)],".")
				model = DataFrame(DBInterface.execute(SQLite.DB(db),q))
				for r in 1:nrow(model)
					println("$(l[col])\t$(db_name)\t$(model[r,Symbol("n.snps.in.model")])\t$(model[r,Symbol("pred.perf.R2")])")
				end
			end
		end
	end
end

function snpSearch(dbs::Array{String,1},snp_files::Array{String,1},col::Int64)
	q = "SELECT * FROM 'weights'"
	# written = false
	for db in dbs
		tiss = split(basename(db),".")[1]
		println(tiss)
		db_name = join(split(basename(db),".")[1:(end-1)],".")
		all_models = DataFrame(DBInterface.execute(SQLite.DB(db),q))

		for chr in snp_files
			out_file = "$(splitext(basename(chr))[1])_weights.txt"
			open(chr) do inf
				open(out_file,"a") do outf
					for line in eachline(inf)
						out_line = chomp(line)
						if (startswith(line, "#"))
							# write(outf,"$(out_line)\ttissue\tgene\teff_allele\tweight\n")
							# written = true
							continue
						end
						l = split(out_line,"\t")
						model = all_models[all_models[!,:rsid].==l[col],:]
						# println(model)
						for r in 1:nrow(model)
							##output-- add db name,gene to end of line
							write(outf,"$(out_line)\t$(tiss)\t$(model[r,Symbol("gene")])\t$(model[r,Symbol("eff_allele")])\t$(model[r,Symbol("weight")])\n")
						end
					end
				end
			end
		end
	end
end

function main()
	parsed_args = parseCommandLine()
	if parsed_args["genes"] != "."
		geneSearch(parsed_args["databases"],parsed_args["genes"],parsed_args["column"])
	else
		snpSearch(parsed_args["databases"],parsed_args["snps"],parsed_args["column"])
	end
end

main()
