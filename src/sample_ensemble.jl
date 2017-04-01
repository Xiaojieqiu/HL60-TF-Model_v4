# include -
include("Include.jl")

# Setup the simulation timescale -
time_start = 0.0
time_stop = 120.0
time_step_size = 0.01

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the previous run -
obj_array = readdlm("error_archive.dat.4")
par_array = readdlm("parameter_archive.dat.4")
obj_array = obj_array[2:end]
par_array = par_array[:,2:end]

# what sort by error -
index_sort_error = sortperm(obj_array)

# Which obj values are *less than or equal* to the median error?
index_of_good_sets = index_sort_error[1:10]

# how many ensemble members do we have?
number_of_samples = length(index_of_good_sets)

# main loop -
for outer_sample_index = 1:number_of_samples

  # get the *local* sample index?
  sample_index = index_of_good_sets[outer_sample_index]

  # what parameters are we looking at?
  parameter_set = par_array[:,sample_index]

  # make a copy of the data dictionary -
  copy_data_dictionary = deepcopy(data_dictionary)

  # how many different types of parameters do we have?
  parameter_mapping_array = copy_data_dictionary["parameter_name_mapping_array"]
  number_of_binding_parameters = length(copy_data_dictionary["binding_parameter_dictionary"])
  number_of_control_parameters = length(copy_data_dictionary["control_parameter_dictionary"])

  # Update data dictionary to match new parameters before running the simulation -
  for index = 1:length(parameter_mapping_array)

    parameter_name = parameter_mapping_array[index]

    if index <= number_of_binding_parameters
      copy_data_dictionary["binding_parameter_dictionary"][parameter_name] = parameter_set[index]
    elseif (index>number_of_binding_parameters && index<=(number_of_binding_parameters+number_of_control_parameters))
      copy_data_dictionary["control_parameter_dictionary"][parameter_name] = parameter_set[index]
    else
      copy_data_dictionary[parameter_name] = parameter_set[index]
    end
  end

  # run the simulation -
  (time_array,simulation_state_array) = add_atra_simulation(time_start,time_stop,time_step_size,copy_data_dictionary);

  # dump data to disk -
  local_data = [time_array simulation_state_array];
  data_filename = "./simulation/sim_data_E1_EI_"*string(sample_index)*"_LI_"*string(outer_sample_index)*".dat"
  writedlm(data_filename,local_data);

  message_string = "Completed $(outer_sample_index) of $(number_of_samples)"
  println(message_string)
end
