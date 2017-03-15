function add_atra_simulation(time_start,time_stop,time_step,data_dictionary)

  # First - run the model to steady-state w/no ATRA forcing -
  XSS = estimate_steady_state(0.001,data_dictionary)

  # Next, set the IC to the steady-state value -
  initial_condition_array = XSS;
  number_of_states = length(XSS)
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Phase 1: Run the model 1/4 of the final time w/o ATRA
  time_start_phase_1 = 0.0;
  time_stop_phase_1 = 1.0;

  # Solve the model equations -
  (TP1,XP1) = SolveBalances(time_start_phase_1,time_stop_phase_1,time_step_size,data_dictionary);

  # set the ic and add the cRaf and activated trigger -
  initial_condition_array = XP1[end,:];
  initial_condition_array[53] = 1.0 # activate trigger
  initial_condition_array[54] = 1.0 # activated cRaf -
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_2 = time_stop_phase_1+time_step_size
  time_stop_phase_2 = time_start_phase_2 + time_stop

  # Solve the model equations -
  (TP2,XP2) = SolveBalances(time_start_phase_2,time_stop_phase_2,time_step_size,data_dictionary);

  # Package the two phases together -
  T = [TP1 ; TP2];
  X = [XP1 ; XP2];

  # return -
  return (T,X);
end

function add_atra_adj_simulation(time_start,time_stop,time_step,data_dictionary)
end
