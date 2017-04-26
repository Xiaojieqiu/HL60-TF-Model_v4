# include -
include("Include.jl")

# Setup the simulation timescale -
time_start = 0.0
time_stop = 48.0
time_step_size = 0.01

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the previous run -
par_array = readdlm("parameter_ensemble.dat")

# how many ensemble members do we have?
(number_of_rows,number_of_samples) = size(par_array)

# main loop -
for outer_sample_index = 1:number_of_samples

  # what parameters are we looking at?
  parameter_set = par_array[:,outer_sample_index]

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

  pairwise_ko_array = zeros(16,16)
  gene_i_index_array = collect(1:1:16)
  gene_j_index_array = collect(1:1:16)
  for gene_i_index in gene_i_index_array
    for gene_j_index in gene_j_index_array

      # run the simulation -
      distance = simulate_pairwise_gene_knockout(time_start,time_stop,time_step_size,gene_i_index,gene_j_index,copy_data_dictionary);

      # grab -
      pairwise_ko_array[gene_i_index,gene_j_index] = distance

      message_string = "Completed $(outer_sample_index) of $(number_of_samples) for (Gi,Gj) = ($(gene_i_index),$(gene_j_index))"
      println(message_string)

    end
  end

  # dump data to disk -
  data_filename = "./pairwise_gene_ko/delta_data_LI_"*string(outer_sample_index)*".dat"
  writedlm(data_filename,pairwise_ko_array);
end
