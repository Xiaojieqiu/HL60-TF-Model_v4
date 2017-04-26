# Script to generate a simple ensemle -
using Distributions
using PyPlot

# load the error and parameter archives -
file_suffix_number_array = [10,13]
top_set_cutoff = 20

# read from disk -
parameter_archive = []
error_archive = Array{Float64,1}()
counter = 1
for file_suffix_index in file_suffix_number_array

  # load -
  raw_error_archive = readdlm("./error_archive.dat."*string(file_suffix_index))
  raw_parameter_archive = readdlm("./parameter_archive.dat."*string(file_suffix_index))

  # cut the lead -
  raw_error_archive = raw_error_archive[2:end]
  raw_parameter_archive = raw_parameter_archive[:,2:end]

  # sort the error archive -
  index_sort_error = sortperm(raw_error_archive )

  # which sets should we grab?
  index_of_good_sets = index_sort_error[1:top_set_cutoff]

  # initialize -
  if counter == 1

    # what is the size?
    (number_of_rows,number_of_cols) = size(raw_parameter_archive)
    parameter_archive = zeros(number_of_rows,1)
  end

  for parameter_set_index in index_of_good_sets
    parameter_set = raw_parameter_archive[:,parameter_set_index]
    parameter_archive = [parameter_archive parameter_set]
  end


  # update the counter =
  counter = counter + 1
end

# cut the leading zero col -
parameter_archive = parameter_archive[:,2:end]

# write to disk -
writedlm("parameter_ensemble.dat",parameter_archive)
