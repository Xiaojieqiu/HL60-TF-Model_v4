function Error(experimental_data_array,simulation_time_array,simulation_state_array,output_index,species_symbol,data_dictionary)

  # initialize -
  total_error_value = Inf

  # Interpole the simulation data to the experimental time scale -
  interpolated_simulated_trajectory = np.interp(experimental_data_array[:,1],simulation_time_array,simulation_state_array[:,output_index]);

  # Compute the error -
  # Error is the sum of a shape function, a scale function and a background_copy_number function -

  # Compute the shape term -
  scaled_data_array = alpha_scaling(experimental_data_array[:,2]);
  scaled_simulation_array = alpha_scaling(interpolated_simulated_trajectory);
  shape_term = sum((scaled_data_array - scaled_simulation_array).^2);

  # Compute the initial value term -
  background_copy_number_dictionary = data_dictionary["background_copy_number_dictionary"];
  background_copy_number = background_copy_number_dictionary[species_symbol]
  cell_volume = data_dictionary["cell_volume"];
  av_number = data_dictionary["av_number"];
  initial_concentration = background_copy_number*(1/av_number)*(1/cell_volume)*1e9; # nM
  initial_value_term = ((interpolated_simulated_trajectory[1] - initial_concentration)^2);

  # Compute the scale term -
  simulated_fold_change = (interpolated_simulated_trajectory)./(interpolated_simulated_trajectory[1]);
  scale_term = sum((simulated_fold_change - experimental_data_array[:,2]).^2);

  # @show (interpolated_simulated_trajectory[1],initial_concentration,interpolated_simulated_trajectory[end],species_symbol)

  # compute the total error -
  total_error_value = shape_term+scale_term+initial_value_term;

  # return -
  return total_error_value
end


function Error_E2F(experimental_data_array,simulation_time_array,simulation_state_array,output_index,species_symbol,data_dictionary)

  # initialize -
  total_error_value = Inf

  # Interpole the simulation data to the experimental time scale -
  interpolated_simulated_trajectory = np.interp(experimental_data_array[:,1],simulation_time_array,simulation_state_array[:,output_index]);

  # Compute G1G0 -
  simulated_G1G0 = 1./interpolated_simulated_trajectory;

  # Compute the shape term -
  scaled_data_array = alpha_scaling(experimental_data_array[:,2]);
  scaled_simulation_array = alpha_scaling(simulated_G1G0);
  shape_term = sum((scaled_data_array - scaled_simulation_array).^2);


  # Compute the initial value term -
  background_copy_number_dictionary = data_dictionary["background_copy_number_dictionary"];
  background_copy_number = background_copy_number_dictionary[species_symbol]
  cell_volume = data_dictionary["cell_volume"];
  av_number = data_dictionary["av_number"];
  initial_concentration = background_copy_number*(1/av_number)*(1/cell_volume)*1e9; # nM
  initial_value_term = ((interpolated_simulated_trajectory[1] - initial_concentration)^2);

  # Compute the fold change term -
  scaling_factor = simulated_G1G0[1];
  simulated_fold_change = (simulated_G1G0)./(scaling_factor);
  scale_term = sum((simulated_fold_change - experimental_data_array[:,2]).^2);

  # compute the total error -
  total_error_value = shape_term+scale_term+initial_value_term;

  # return -
  return total_error_value;
end

function alpha_scaling(data_array)

  # get the dimension -
  number_of_points = length(data_array)

  # what is the min index?
  min_index = indmin(data_array)
  min_value = data_array[min_index]

  # What is the max index?
  max_index = indmax(data_array)
  max_value = data_array[max_index]

  # compute the scaling -
  scaled_array = Float64[]
  for index = 1:number_of_points

    value = data_array[index]
    scaled_value = (value - min_value)/(max_value - min_value)
    push!(scaled_array,scaled_value)
  end

  return scaled_array
end
