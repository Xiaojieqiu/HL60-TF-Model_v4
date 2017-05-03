# script to randomly perturb the parameters -

# parameters for sampling -
number_of_samples = 10
coefficient_of_variation = 0.10

# load the parameters -
parameter_array = readdlm("parameter_archive.dat.23")
parameter_array = parameter_array[:,2:end]

error_archive = readdlm("error_archive.dat.23")
error_archive = error_archive[2:end]

idx_error_sort = sortperm(error_archive)
mean_parameter_array = mean(parameter_array[:,idx_error_sort[1:10]],2)

# how many parameters?
number_of_parameters = length(mean_parameter_array)
ensemble_array = zeros(number_of_parameters,number_of_samples)
for sample_index = 1:number_of_samples

  # create new array -
  new_parameter_array = mean_parameter_array.*(1+coefficient_of_variation*randn(number_of_parameters))

  # copy -
  for parameter_index = 1:number_of_parameters
    ensemble_array[parameter_index,sample_index] = new_parameter_array[parameter_index]
  end

end

# dump to disk -
writedlm("randomized_parameter_ensemble.dat",ensemble_array)
