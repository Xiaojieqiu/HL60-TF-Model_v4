# include -
include("Include.jl")

# Setup the simulation timescale -
time_start = 0.0
time_stop = 120.0
time_step_size = 0.01

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the data_dictionary -
# initial_parameter_guess = Float64[]
# parameter_mapping_array = data_dictionary["parameter_name_mapping_array"]
# number_of_binding_parameters = length(data_dictionary["binding_parameter_dictionary"])
# number_of_control_parameters = length(data_dictionary["control_parameter_dictionary"])
#
# # Update data dictionary to match new parameters before calculating obj
# for index = 1:length(parameter_mapping_array)
#
#   parameter_name = parameter_mapping_array[index]
#
#   if index <= number_of_binding_parameters
#     push!(initial_parameter_guess,data_dictionary["binding_parameter_dictionary"][parameter_name])
#   elseif (index>number_of_binding_parameters && index<=(number_of_binding_parameters+number_of_control_parameters))
#     push!(initial_parameter_guess,data_dictionary["control_parameter_dictionary"][parameter_name])
#   else
#     push!(initial_parameter_guess,data_dictionary[parameter_name])
#   end
# end

# get my initial parameter guess from the previous run -
obj_array = readdlm("objective_archive.dat.8")
par_array = readdlm("parameter_archive.dat.8")
min_index = indmin(obj_array)
initial_parameter_guess = par_array[:,min_index]

# Search exposes the *run loop* method -
(objective_archive,parameter_archive) = estimate_model_parameters(objective_function,generation_function,acceptance_function,constraints_function,
  initial_parameter_guess; maximum_number_of_iterations=1000,show_trace=true)

# convert wrapper to actual array -
number_of_parameters = length(initial_parameter_guess)
parameter_array = zeros(number_of_parameters,1)
for parameter_wrapper in parameter_archive

  parameter_vector = parameter_wrapper.array
  parameter_array = [parameter_array parameter_vector]
end
parameter_array = parameter_array[:,2:end]

# write results to disk -
writedlm("objective_archive.dat.9",objective_archive)
writedlm("parameter_archive.dat.9",parameter_array)
