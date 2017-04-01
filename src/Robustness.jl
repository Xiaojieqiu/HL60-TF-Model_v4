function simulate_pairwise_gene_knockout(time_start,time_stop,time_step,gene_i_index,gene_j_index,data_dictionary)

  # get the initial condition array -
  initial_condition_array = data_dictionary["initial_condition_array"]
  initial_condition_array[gene_i_index] = 0.0
  initial_condition_array[gene_j_index] = 0.0

  # solve the model equations -
  (T,X) = add_atra_simulation(time_start,time_stop,time_step_size,data_dictionary);

  return (T,X)
end
