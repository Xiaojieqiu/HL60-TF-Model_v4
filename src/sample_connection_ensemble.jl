# include -
include("Include.jl")

# Setup the simulation timescale -
time_start = 0.0
time_stop = 48.0
time_step_size = 0.01

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the previous run -
par_array = readdlm("randomized_parameter_ensemble.dat")

# how many ensemble members do we have?
(number_of_rows,number_of_samples) = size(par_array)

# which parameters do we want to KO -
list_of_ko_parameters = [
  "W_gene_AP1_RNAP"	;	# 91
  "W_gene_AP1_gene_AhR"	;	# 92
  "W_gene_AP1_gene_PU1"	;	# 93
  "W_gene_AP1_gene_PPARg"	;	# 94
  "W_gene_AhR_RNAP"	;	# 95
  "W_gene_AhR_gene_Trigger"	;	# 96
  "W_gene_CD11b_RNAP"	;	# 97
  "W_gene_CD11b_gene_PU1_gene_cRAF"	;	# 98
  "W_gene_CD14_RNAP"	;	# 99
  "W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 100
  "W_gene_CD38_RNAP"	;	# 101
  "W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 102
  "W_gene_CEBPa_RNAP"	;	# 103
  "W_gene_CEBPa_gene_Trigger"	;	# 104
  "W_gene_CEBPa_gene_PPARg"	;	# 105
  "W_gene_CEBPa_gene_CEBPa"	;	# 106
  "W_gene_CEBPa_gene_GFI1"	;	# 107
  "W_gene_E2F_RNAP"	;	# 108
  "W_gene_E2F_gene_E2F"	;	# 109
  "W_gene_E2F_gene_PPARg"	;	# 110
  "W_gene_E2F_gene_CEBPa"	;	# 111
  "W_gene_E2F_gene_GFI1"	;	# 112
  "W_gene_E2F_gene_cRAF"	;	# 113
  "W_gene_EGR1_RNAP"	;	# 114
  "W_gene_EGR1_gene_Trigger"	;	# 115
  "W_gene_EGR1_gene_PU1"	;	# 116
  "W_gene_EGR1_gene_PPARg"	;	# 117
  "W_gene_EGR1_gene_GFI1"	;	# 118
  "W_gene_GFI1_RNAP"	;	# 119
  "W_gene_GFI1_gene_CEBPa"	;	# 120
  "W_gene_GFI1_gene_EGR1"	;	# 121
  "W_gene_IRF1_RNAP"	;	# 122
  "W_gene_IRF1_gene_Trigger"	;	# 123
  "W_gene_IRF1_gene_AhR"	;	# 124
  "W_gene_IRF1_gene_PPARg"	;	# 125
  "W_gene_OCT1_RNAP"	;	# 126
  "W_gene_OCT1_gene_PPARg"	;	# 127
  "W_gene_OCT4_RNAP"	;	# 128
  "W_gene_OCT4_gene_Trigger"	;	# 129
  "W_gene_OCT4_gene_AhR"	;	# 130
  "W_gene_OCT4_gene_cRAF"	;	# 131
  "W_gene_P21_RNAP"	;	# 132
  "W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 133
  "W_gene_P21_gene_GFI1"	;	# 134
  "W_gene_P47Phox_RNAP"	;	# 135
  "W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 136
  "W_gene_P47Phox_gene_PPARg"	;	# 137
  "W_gene_PPARg_RNAP"	;	# 138
  "W_gene_PPARg_gene_Trigger"	;	# 139
  "W_gene_PPARg_gene_CEBPa"	;	# 140
  "W_gene_PPARg_gene_EGR1"	;	# 141
  "W_gene_PPARg_gene_PU1"	;	# 142
  "W_gene_PPARg_gene_AP1"	;	# 143
  "W_gene_PU1_RNAP"	;	# 144
  "W_gene_PU1_gene_Trigger"	;	# 145
  "W_gene_PU1_gene_CEBPa"	;	# 146
  "W_gene_PU1_gene_PU1"	;	# 147
  "W_gene_PU1_gene_AP1"	;	# 148
  "W_gene_PU1_gene_OCT1"	;	# 149
  "W_gene_PU1_gene_AhR"	;	# 150
  "W_gene_PU1_gene_GFI1"	;	# 151
];

number_of_ko_connection_parameters = length(list_of_ko_parameters)

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
  data_filename = "./pairwise_connection_ko/delta_data_"*string(outer_sample_index)*".dat"
  writedlm(data_filename,pairwise_cko_array);
end
