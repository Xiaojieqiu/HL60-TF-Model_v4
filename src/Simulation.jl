function add_atra_simulation(time_start,time_stop,time_step,data_dictionary)

  # Run the model to steady-state -
  steady_state_array = estimate_steady_state(0.01,data_dictionary);

  # Ok, reset the TRIGGER and CRAF -
  steady_state_array[53] = 1;
  steady_state_array[54] = 1;

  #@show steady_state_array[37:54]

  # Reset the IC's -
  data_dictionary["initial_condition_array"] = steady_state_array;

  # Solve the model equations -
  (T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)

  # return my simulation time, and state -
  return (T,X)

end

function add_atra_adj_simulation(time_start,time_stop,time_step,data_dictionary)
end
