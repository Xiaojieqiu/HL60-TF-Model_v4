function simulate_pairwise_gene_knockout(time_start,time_stop,time_step,gene_i_index,gene_j_index,data_dictionary)

  # Solve the nominal -
  (TN,XN) = add_atra_simulation(time_start,time_stop,time_step_size,data_dictionary);

  # get the initial condition array -
  copy_data_dictionary = deepcopy(data_dictionary)
  initial_condition_array = copy_data_dictionary["initial_condition_array"]
  initial_condition_array[gene_i_index] = 0.0
  initial_condition_array[gene_j_index] = 0.0

  # solve the model equations -
  (TP,XP) = add_atra_simulation(time_start,time_stop,time_step_size,copy_data_dictionary);

  # compute a scale factor -
  SF = norm(XN)

  # calculate the scaled differnce -
  return (1/SF)*norm((XP-XN))
end

function simulate_pairwise_connection_knockout(time_start,time_stop,time_step,connection_symbol_i,connection_symbol_j,data_dictionary)

  # Solve the nominal -
  (TN,XN) = add_atra_simulation(time_start,time_stop,time_step_size,data_dictionary);

  # get the initial condition array -
  copy_data_dictionary = deepcopy(data_dictionary)
  control_parameter_dictionary = copy_data_dictionary["control_parameter_dictionary"]
  control_parameter_dictionary[connection_symbol_i] = 0.0
  control_parameter_dictionary[connection_symbol_j] = 0.0

  # solve the model equations -
  (TP,XP) = add_atra_simulation(time_start,time_stop,time_step_size,copy_data_dictionary);

  # compute a scale factor -
  SF = norm(XN)

  # calculate the scaled differnce -
  return (1/SF)*norm((XP-XN))
end
