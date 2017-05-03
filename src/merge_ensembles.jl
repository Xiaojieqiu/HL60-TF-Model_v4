# Script to merge ensembles -

# load ensemble 1 -
ensemble_base_error = readdlm("error_archive.dat.21")
ensemble_base_parameters = readdlm("parameter_archive.dat.21")
ensemble_base_error = ensemble_base_error[2:end]
ensemble_base_parameters = ensemble_base_parameters[:,2:end]

# load ensemble 2 -
ensemble_new_error = readdlm("error_archive.dat.22")
ensemble_new_parameters = readdlm("parameter_archive.dat.22")
ensemble_new_error = ensemble_new_error[2:end]
ensemble_new_parameters = ensemble_new_parameters[:,2:end]

# How many new parameter sets do we want?
number_of_parameter_sets = 10
(number_of_parameters,number_of_samples) = size(ensemble_new_parameters)

# grab the number_of_parameter_sets best sets from the base -
idx_base_sort = sortperm(ensemble_base_error,rev=false)
base_parameters = ensemble_base_parameters[:,idx_base_sort[1:number_of_parameter_sets]]

idx_new_sort = sortperm(ensemble_new_error,rev=false)
new_parameters = ensemble_new_parameters[:,idx_new_sort[1:number_of_parameter_sets]]

# ok, lets create the merged ensemble -
merged_parameter_array = zeros(number_of_parameters,number_of_parameter_sets)
idx_patch = [collect(23:32) ; collect(108:113)]
number_of_patch_parameters = length(idx_patch)
for parameter_set_index = 1:number_of_parameter_sets

  base_set = base_parameters[:,parameter_set_index]
  new_set = new_parameters[:,parameter_set_index]

  # copy -
  for parameter_index = 1:number_of_parameters
    merged_parameter_array[parameter_index,parameter_set_index] = base_set[parameter_index]
  end

  # patch -
  for patch_index in idx_patch
    merged_parameter_array[patch_index,parameter_set_index] = new_set[patch_index]
  end
end

# dump to disk -
writedlm("parameter_ensemble.dat",merged_parameter_array)
