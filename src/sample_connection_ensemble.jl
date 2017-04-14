# include -
include("Include.jl")

# Setup the simulation timescale -
time_start = 0.0
time_stop = 48.0
time_step_size = 0.01

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the previous run -
obj_array = readdlm("error_archive.dat.4")
par_array = readdlm("parameter_archive.dat.4")
obj_array = obj_array[2:end]
par_array = par_array[:,2:end]

# sort by error -
index_sort_error = sortperm(obj_array)

# Take the top 10 -
index_of_good_sets = index_sort_error[1:10]

# how many ensemble members do we have?
number_of_samples = length(index_of_good_sets)

# which parameters do we want to KO -
list_of_ko_parameters = [

  "W_gene_AP1_gene_PU1" ; # 1
  "W_gene_CD11b_gene_PU1_gene_cRAF" ; # 2
  "W_gene_EGR1_gene_PU1"  ; # 3
  "W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF" ; # 4
  "W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"  ; # 5
  "W_gene_PPARg_gene_PU1" ; # 6
  "W_gene_PU1_gene_PU1" ; # 7
  "W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF" ; # 8
  "W_gene_PPARg_gene_AP1" ; # 9
  "W_gene_PU1_gene_AP1" ; # 10
];



number_of_ko_connection_parameters = length(list_of_ko_parameters)

# main loop -
for outer_sample_index = 1:1

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


  pairwise_cko_array = zeros(number_of_ko_connection_parameters,number_of_ko_connection_parameters)
  for outer_connection_index = 1:number_of_ko_connection_parameters

    # grab the outer symbol -
    outer_connection_symbol = list_of_ko_parameters[outer_connection_index]
    for inner_connection_index = 1:number_of_ko_connection_parameters

      # grab the inner symbol -
      inner_connection_symbol = list_of_ko_parameters[inner_connection_index]

      distance = simulate_pairwise_connection_knockout(time_start,time_stop,time_step_size,outer_connection_symbol,inner_connection_symbol,data_dictionary)

      # grab -
      pairwise_cko_array[outer_connection_index,inner_connection_index] = distance

      # let the user know progress -
      message_string = "Completed $(outer_sample_index) of $(number_of_samples) for (Wi,Wj) = ($(outer_connection_symbol),$(inner_connection_symbol))"
      println(message_string)

    end
  end

  # dump data to disk -
  data_filename = "./pairwise_connection_ko/delta_data_EI_"*string(sample_index)*"_LI_"*string(outer_sample_index)*".dat"
  writedlm(data_filename,pairwise_cko_array);
end
