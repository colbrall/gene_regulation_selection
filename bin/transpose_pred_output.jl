# copied from 1kG_analyses.jl
# transposes predicted_expression.txt files (predixcan output) into gene x ind matrix
# this is the format needed for a lot of the population-level analyses already written
# julia 1.4

using DataFrames
using CSV

# reads file into DataFrame
function readDF(f_path::String)
  println("Reading file $f_path as data frame.......")
  if ispath(f_path)
    df = CSV.read(f_path,DataFrame;delim='\t')
  else
    df = "NA"
  end
  return df
end


# converts predicted_expression.txt into an Eric-style file
# takes root directory that holds output directories
function transPredExp(path::String)
  for item in readdir(path)
    if isdir(realpath("$path/$item"))
      out_f = "$(realpath(path))/$(item)_elasticNet0_0.5.full"
      if ispath("$(out_f).gz")
        println("$out_f already exists!")
        continue
      end
      prd_exp = readDF("$path$item/predicted_expression.txt")
      if prd_exp == "NA"
        println("No predicted_expression.txt file in $path$item/")
        continue
      end
      select!(prd_exp, Not(:FID))
      genes = names(prd_exp)
      prd_exp = unstack(stack(prd_exp,names(prd_exp)[2:end]),:variable,:IID,:value)
      rename!(prd_exp,:variable => :gene)
      CSV.write(out_f, prd_exp, delim='\t')
      run(`gzip $out_f`)
    end
  end
end


function main()
#### Transposing predicted_expression.txt file
    in_files = ARGS[:,1]::Array{String,1}
    transPredExp(in_files[length(in_files)])
end

main()
