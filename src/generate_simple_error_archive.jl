# script to generate an error archive -
include("Include.jl")

# load the parameter ensemble -
par_array = readdlm("randomized_parameter_ensemble.dat")
(number_of_parameters,number_of_samples) = size(par_array)


error_archive = Array{Float64,1}()
for sample_index = 1:number_of_samples

  parameter_array = par_array[:,sample_index]
  error_value = objective_function(vec(parameter_array))

  push!(error_archive,error_value)

  @show sample_index

end

# dump -
writedlm("error_archive_random_ensemble.dat",error_archive)
