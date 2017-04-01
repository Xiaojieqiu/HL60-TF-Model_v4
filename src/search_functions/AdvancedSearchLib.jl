# ================ GLOBAL INITIALIZATION CODE DO NOT EDIT ====================================== %
# some global parameters -
BIG = 1e10
SMALL = 1e-8

# Load the experiment specifications -
tmp_value = JSON.parsefile("./search_functions/experiments/Experiments.json")
experiment_array = tmp_value["experiment_array"]

# preload the data -
cached_data_dictionary = Dict()
for (experiment_index,experiment_dictionary) in enumerate(experiment_array)

  # Grab the data file -
  data_file_path = experiment_dictionary["data_file"]
  experiment_id = experiment_dictionary["experiment_id"]

  # Load the data -
  experimental_data_array = readdlm(data_file_path)

  # Cache the data w/experiment_id -
  cached_data_dictionary[experiment_id] = experimental_data_array
end

# preload the error function calls -
cached_error_function_dictionary = Dict()
for (experiment_index,experiment_dictionary) in enumerate(experiment_array)

  # Grab the experimental id -
  experiment_id = experiment_dictionary["experiment_id"]

  # Grab the error functions from the experiment_dictionary -
  error_function = eval(parse(experiment_dictionary["error_function"]))

  # cache these functions -
  cached_error_function_dictionary[experiment_id] = error_function
end
# ================ GLOBAL INITIALIZATION CODE DO NOT EDIT ====================================== %

function objective_function(parameter_guess::Array{Float64,1})

  # Script to solve the balance equations -
  time_start = 0.0
  time_stop = 120.0
  time_step_size = 0.01

  # Load the data dictionary -
  data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

  # how many different types of parameters do we have?
  parameter_mapping_array = data_dictionary["parameter_name_mapping_array"]
  number_of_binding_parameters = length(data_dictionary["binding_parameter_dictionary"])
  number_of_control_parameters = length(data_dictionary["control_parameter_dictionary"])

  # Update data dictionary to match new parameters before calculating obj
  for index = 1:length(parameter_mapping_array)

    parameter_name = parameter_mapping_array[index]

    if index <= number_of_binding_parameters
      data_dictionary["binding_parameter_dictionary"][parameter_name] = parameter_guess[index]
    elseif (index>number_of_binding_parameters && index<=(number_of_binding_parameters+number_of_control_parameters))
      data_dictionary["control_parameter_dictionary"][parameter_name] = parameter_guess[index]
    else
      data_dictionary[parameter_name] = parameter_guess[index]
    end
  end

  # how many objectives do we have?
  obj_array = zeros(13)

  # Call simulation routine -
  # (run the model to SS, and then set the ICs to the SS for this parameter set)
  (time_array,simulation_state_array) = add_atra_simulation(time_start,time_stop,time_step_size,data_dictionary);

  wght_array = ones(13)
  wght_array[4] = 5.4 # PU1
  wght_array[5] = 5.2 # P47
  wght_array[9] = 2.84 # IRF1
  wght_array[10] = 2.43 # EGR1


  # Call the error functions -
  # loop through the experimental dictionary, and call the appropriate error function -
  for (experiment_index,experiment_dictionary) in enumerate(experiment_array)

    # Get the error function pointer -
    error_function_pointer = experiment_dictionary["error_function"]
    output_index = parse(Int,experiment_dictionary["output_index"])
    experiment_id = experiment_dictionary["experiment_id"]
    species_symbol = experiment_dictionary["protein_symbol"];

    # Get the experimental data array -
    experimental_data_array = cached_data_dictionary[experiment_id]

    # Call the error function -
    error_function = cached_error_function_dictionary[experiment_id]
    error_value = error_function(experimental_data_array,time_array,simulation_state_array,output_index,species_symbol,data_dictionary)

    # Add the error to the objective array -
    obj_array[experiment_index] = wght_array[experiment_index]*error_value;
  end




  return sum(obj_array)
end

function sample_function(parameter_guess,kVI,NPARAMETERS,A,SIGMA)

  # Generate a new vector -
  zV = randn(NPARAMETERS);

  # Compute the new perturbed parameter vector -
  perturbed_parameter_array = parameter_guess.*(1 + SIGMA*A*zV);
  corrected_parameter_array = parameter_bounds_function(perturbed_parameter_array);

  # return -
  return (corrected_parameter_array,zV)
end


function parameter_bounds_function(parameter_guess::Array{Float64,1})

  # How many parameters do we have?
  number_of_parameters = length(parameter_guess)

  # Setup my upper bound, and lower bounds on parameters -
  bounds_array = [

    # binding parameters -
    0.0 4.0    ; # "n_gene_AP1_gene_AhR"	;	# 1
		0.0 1000.0  ; # "K_gene_AP1_gene_AhR"	;	# 2
		0.0 4.0    ; # "n_gene_AP1_gene_PU1"	;	# 3
		0.0 1000.0  ; # "K_gene_AP1_gene_PU1"	;	# 4
		0.0 4.0    ; # "n_gene_AP1_gene_PPARg"	;	# 5
		0.0 1000.0  ; # "K_gene_AP1_gene_PPARg"	;	# 6
		0.0 4.0    ; # "n_gene_AhR_gene_Trigger"	;	# 7
		0.0 1000.0  ; # "K_gene_AhR_gene_Trigger"	;	# 8
		0.0 4.0    ; # "n_gene_CD11b_gene_PU1_gene_cRAF"	;	# 9
		0.0 1000.0  ; # "K_gene_CD11b_gene_PU1_gene_cRAF"	;	# 10
		0.0 4.0    ; # "n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 11
		0.0 1000.0  ; # "K_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 12
		0.0 4.0    ; # "n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 13
		0.0 1000.0  ; # "K_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 14
		0.0 4.0    ; # "n_gene_CEBPa_gene_Trigger"	;	# 15
		0.0 1000.0  ; # "K_gene_CEBPa_gene_Trigger"	;	# 16
		0.0 4.0    ; # "n_gene_CEBPa_gene_PPARg"	;	# 17
		0.0 1000.0  ; # "K_gene_CEBPa_gene_PPARg"	;	# 18
		0.0 4.0    ; # "n_gene_CEBPa_gene_CEBPa"	;	# 19
		0.0 1000.0  ; # "K_gene_CEBPa_gene_CEBPa"	;	# 20
		0.0 4.0    ; # "n_gene_CEBPa_gene_GFI1"	;	# 21
		0.0 1000.0  ; # "K_gene_CEBPa_gene_GFI1"	;	# 22
		0.0 4.0    ; # "n_gene_E2F_gene_E2F"	;	# 23
		0.0 1000.0  ; # "K_gene_E2F_gene_E2F"	;	# 24
		0.0 4.0    ; # "n_gene_E2F_gene_PPARg"	;	# 25
		0.0 1000.0  ; # "K_gene_E2F_gene_PPARg"	;	# 26
		0.0 4.0    ; # "n_gene_E2F_gene_CEBPa"	;	# 27
		0.0 1000.0  ; # "K_gene_E2F_gene_CEBPa"	;	# 28
		0.0 4.0    ; # "n_gene_E2F_gene_GFI1"	;	# 29
		0.0 1000.0  ; # "K_gene_E2F_gene_GFI1"	;	# 30
		0.0 4.0    ; # "n_gene_E2F_gene_cRAF"	;	# 31
		0.0 1000.0  ; # "K_gene_E2F_gene_cRAF"	;	# 32
		0.0 4.0    ; # "n_gene_EGR1_gene_Trigger"	;	# 33
		0.0 1000.0  ; # "K_gene_EGR1_gene_Trigger"	;	# 34
		0.0 4.0    ; # "n_gene_EGR1_gene_PU1"	;	# 35
		0.0 1000.0  ; # "K_gene_EGR1_gene_PU1"	;	# 36
		0.0 4.0    ; # "n_gene_EGR1_gene_PPARg"	;	# 37
		0.0 1000.0  ; # "K_gene_EGR1_gene_PPARg"	;	# 38
		0.0 4.0    ; # "n_gene_EGR1_gene_GFI1"	;	# 39
		0.0 1000.0  ; # "K_gene_EGR1_gene_GFI1"	;	# 40
		0.0 4.0    ; # "n_gene_GFI1_gene_CEBPa"	;	# 41
		0.0 1000.0  ; # "K_gene_GFI1_gene_CEBPa"	;	# 42
		0.0 4.0    ; # "n_gene_GFI1_gene_EGR1"	;	# 43
		0.0 1000.0  ; # "K_gene_GFI1_gene_EGR1"	;	# 44
		0.0 4.0    ; # "n_gene_IRF1_gene_Trigger"	;	# 45
		0.0 1000.0  ; # "K_gene_IRF1_gene_Trigger"	;	# 46
		0.0 4.0    ; # "n_gene_IRF1_gene_AhR"	;	# 47
		0.0 1000.0  ; # "K_gene_IRF1_gene_AhR"	;	# 48
		0.0 4.0    ; # "n_gene_IRF1_gene_PPARg"	;	# 49
		0.0 1000.0  ; # "K_gene_IRF1_gene_PPARg"	;	# 50
		0.0 4.0    ; # "n_gene_OCT1_gene_PPARg"	;	# 51
		0.0 1000.0  ; # "K_gene_OCT1_gene_PPARg"	;	# 52
		0.0 4.0    ; # "n_gene_OCT4_gene_Trigger"	;	# 53
		0.0 1000.0  ; # "K_gene_OCT4_gene_Trigger"	;	# 54
		0.0 4.0    ; # "n_gene_OCT4_gene_AhR"	;	# 55
		0.0 1000.0  ; # "K_gene_OCT4_gene_AhR"	;	# 56
		0.0 4.0    ; # "n_gene_OCT4_gene_cRAF"	;	# 57
		0.0 1000.0  ; # "K_gene_OCT4_gene_cRAF"	;	# 58
		0.0 4.0    ; # "n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 59
		0.0 1000.0  ; # "K_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 60
		0.0 4.0    ; # "n_gene_P21_gene_GFI1"	;	# 61
		0.0 1000.0  ; # "K_gene_P21_gene_GFI1"	;	# 62
		0.0 4.0    ; # "n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 63
		0.0 1000.0  ; # "K_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 64
		0.0 4.0    ; # "n_gene_P47Phox_gene_PPARg"	;	# 65
		0.0 1000.0  ; # "K_gene_P47Phox_gene_PPARg"	;	# 66
		0.0 4.0    ; # "n_gene_PPARg_gene_Trigger"	;	# 67
		0.0 1000.0  ; # "K_gene_PPARg_gene_Trigger"	;	# 68
		0.0 4.0    ; # "n_gene_PPARg_gene_CEBPa"	;	# 69
		0.0 1000.0  ; # "K_gene_PPARg_gene_CEBPa"	;	# 70
		0.0 4.0    ; # "n_gene_PPARg_gene_EGR1"	;	# 71
		0.0 1000.0  ; # "K_gene_PPARg_gene_EGR1"	;	# 72
		0.0 4.0    ; # "n_gene_PPARg_gene_PU1"	;	# 73
		0.0 1000.0  ; # "K_gene_PPARg_gene_PU1"	;	# 74
		0.0 4.0    ; # "n_gene_PPARg_gene_AP1"	;	# 75
		0.0 1000.0  ; # "K_gene_PPARg_gene_AP1"	;	# 76
		0.0 4.0    ; # "n_gene_PU1_gene_Trigger"	;	# 77
		0.0 1000.0  ; # "K_gene_PU1_gene_Trigger"	;	# 78
		0.0 4.0    ; # "n_gene_PU1_gene_CEBPa"	;	# 79
		0.0 1000.0  ; # "K_gene_PU1_gene_CEBPa"	;	# 80
		0.0 4.0    ; # "n_gene_PU1_gene_PU1"	;	# 81
		0.0 1000.0  ; # "K_gene_PU1_gene_PU1"	;	# 82
		0.0 4.0    ; # "n_gene_PU1_gene_AP1"	;	# 83
		0.0 1000.0  ; # "K_gene_PU1_gene_AP1"	;	# 84
		0.0 4.0    ; # "n_gene_PU1_gene_OCT1"	;	# 85
		0.0 1000.0  ; # "K_gene_PU1_gene_OCT1"	;	# 86
		0.0 4.0    ; # "n_gene_PU1_gene_AhR"	;	# 87
		0.0 1000.0  ; # "K_gene_PU1_gene_AhR"	;	# 88
		0.0 4.0    ; # "n_gene_PU1_gene_GFI1"	;	# 89
		0.0 1000.0  ; # "K_gene_PU1_gene_GFI1"	;	# 90

    # weight parameters -
		0 100 ; # "W_gene_AP1_RNAP"	;	# 91
		0 100 ; # "W_gene_AP1_gene_AhR"	;	# 92
		0 100 ; # "W_gene_AP1_gene_PU1"	;	# 93
		0 100 ; # "W_gene_AP1_gene_PPARg"	;	# 94
		0 100 ; # "W_gene_AhR_RNAP"	;	# 95
		0 100 ; # "W_gene_AhR_gene_Trigger"	;	# 96
		0 100 ; # "W_gene_CD11b_RNAP"	;	# 97
		0 100 ; # "W_gene_CD11b_gene_PU1_gene_cRAF"	;	# 98
		0 0 ; # "W_gene_CD14_RNAP"	;	# 99
		0 0 ; # "W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 100
		0 100 ; # "W_gene_CD38_RNAP"	;	# 101
		0 100 ; # "W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 102
		0 100 ; # "W_gene_CEBPa_RNAP"	;	# 103
		0 100 ; # "W_gene_CEBPa_gene_Trigger"	;	# 104
		0 100 ; # "W_gene_CEBPa_gene_PPARg"	;	# 105
		0 100 ; # "W_gene_CEBPa_gene_CEBPa"	;	# 106
		0 100 ; # "W_gene_CEBPa_gene_GFI1"	;	# 107
		0 100 ; # "W_gene_E2F_RNAP"	;	# 108
		0 100 ; # "W_gene_E2F_gene_E2F"	;	# 109
		0 100 ; # "W_gene_E2F_gene_PPARg"	;	# 110
		0 100 ; # "W_gene_E2F_gene_CEBPa"	;	# 111
		0 100 ; # "W_gene_E2F_gene_GFI1"	;	# 112
		0 100 ; # "W_gene_E2F_gene_cRAF"	;	# 113
		0 100 ; # "W_gene_EGR1_RNAP"	;	# 114
		0 100 ; # "W_gene_EGR1_gene_Trigger"	;	# 115
		0 100 ; # "W_gene_EGR1_gene_PU1"	;	# 116
		0 100 ; # "W_gene_EGR1_gene_PPARg"	;	# 117
		0 100 ; # "W_gene_EGR1_gene_GFI1"	;	# 118
		0 100 ; # "W_gene_GFI1_RNAP"	;	# 119
		0 100 ; # "W_gene_GFI1_gene_CEBPa"	;	# 120
		0 100 ; # "W_gene_GFI1_gene_EGR1"	;	# 121
		0 100 ; # "W_gene_IRF1_RNAP"	;	# 122
		0 100 ; # "W_gene_IRF1_gene_Trigger"	;	# 123
		0 100 ; # "W_gene_IRF1_gene_AhR"	;	# 124
		0 100 ; # "W_gene_IRF1_gene_PPARg"	;	# 125
		0 100 ; # "W_gene_OCT1_RNAP"	;	# 126
		0 100 ; # "W_gene_OCT1_gene_PPARg"	;	# 127
		0 100 ; # "W_gene_OCT4_RNAP"	;	# 128
		0 100 ; # "W_gene_OCT4_gene_Trigger"	;	# 129
		0 100 ; # "W_gene_OCT4_gene_AhR"	;	# 130
		0 100 ; # "W_gene_OCT4_gene_cRAF"	;	# 131
		0 100 ; # "W_gene_P21_RNAP"	;	# 132
		0 100 ; # "W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 133
		0 100 ; # "W_gene_P21_gene_GFI1"	;	# 134
		0 100 ; # "W_gene_P47Phox_RNAP"	;	# 135
		0 100 ; # "W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 136
		0 100 ; # "W_gene_P47Phox_gene_PPARg"	;	# 137
		0 100 ; # "W_gene_PPARg_RNAP"	;	# 138
		0 100 ; # "W_gene_PPARg_gene_Trigger"	;	# 139
		0 100 ; # "W_gene_PPARg_gene_CEBPa"	;	# 140
		0 100 ; # "W_gene_PPARg_gene_EGR1"	;	# 141
		0 100 ; # "W_gene_PPARg_gene_PU1"	;	# 142
		0 100 ; # "W_gene_PPARg_gene_AP1"	;	# 143
		0 100 ; # "W_gene_PU1_RNAP"	;	# 144
		0 100 ; # "W_gene_PU1_gene_Trigger"	;	# 145
		0 100 ; # "W_gene_PU1_gene_CEBPa"	;	# 146
		0 100 ; # "W_gene_PU1_gene_PU1"	;	# 147
		0 100 ; # "W_gene_PU1_gene_AP1"	;	# 148
		0 100 ; # "W_gene_PU1_gene_OCT1"	;	# 149
		0 100 ; # "W_gene_PU1_gene_AhR"	;	# 150
		0 100 ; # "W_gene_PU1_gene_GFI1"	;	# 151
		0 0 ; # "W_gene_Trigger_RNAP"	;	# 152
		0 0 ; # "W_gene_cRAF_RNAP"	;	# 153
	];

  # Split into lower and upper bound arrays -
  lower_bound_array = bounds_array[:,1]
  upper_bound_array = bounds_array[:,2]

  # iterate through and fix the parameters -
  new_parameter_array = copy(parameter_guess)
  for (index,value) in enumerate(parameter_guess)

    lower_bound = lower_bound_array[index]
    upper_bound = upper_bound_array[index]

    if (value<lower_bound)
      new_parameter_array[index] = lower_bound
    elseif (value>upper_bound)
      new_parameter_array[index] = upper_bound
    end
  end

  return new_parameter_array
end

function updateCholesky(A,zV,CCOV,AVG_P_SUCC,PTHRESH)

  if (AVG_P_SUCC<PTHRESH)
		CA = sqrt(1-CCOV);
		F = 1 + ((1-CA^2)*(norm(zV))^2)/(CA^2);
		A = CA*A+CA/(norm(zV)^2)*(sqrt(F)-1)*A*zV*transpose(zV);
	else
		A = A;
	end

  return A;
end

function updateStepSize(LAMBDA_SUCC,AVG_P_SUCC,P_SUCC_TARGET,CP,D,SIGMA)

	# Compute AVG_P_SUCC -
	AVG_P_SUCC = (1-CP)*AVG_P_SUCC+CP*LAMBDA_SUCC;

	# Compute the step-size -
	F = (1/D)*(AVG_P_SUCC-(P_SUCC_TARGET/(1-P_SUCC_TARGET))*(1-AVG_P_SUCC));
	SIGMA = SIGMA*exp(F);

	if (SIGMA<= 1e-5)
		SIGMA = 1e-5;
	end

  return (SIGMA,AVG_P_SUCC);
end

function estimate_parameters(pObjectiveFunction,initial_parameter_array,data_dictionary,number_of_iterations,error_target)

  # Get the initial parameter set -
  kVI = initial_parameter_array;
  N = length(kVI);

  # Initialize the algorithm parameters -
  D = 1+N/2;
  P_SUCC_TARGET = 2/11;
  CP = 1/12;
  CC = 2/(N+2);
  CCOV = 2/(N^2+6);
  PTHRESH = 0.99;

  # Set the initial search directions - apply -
  A = eye(N,N);
  AVG_P_SUCC = P_SUCC_TARGET;

  # Setup the run -
  should_loop_continue = true
  parameter_array = kVI;
  SIGMA = 0.10;

  # Initialize -
  error_array = zeros(2)
  error_array[1] = Inf
  counter = 1;

  # How many failed steps are we allowed to take?
  number_of_failed_steps = 0;
  max_number_of_failed_steps = 10;

  # setup archive -
  parameter_archive = zeros(N,1)
  error_archive = 0;
  while (should_loop_continue == true)

    # Generate new parameter set -
  	(perturbed_parameter_array,zV) = sample_function(parameter_array,initial_parameter_array,N,A,SIGMA);

    # Evaluate the parameter array -
    objective_function_array = pObjectiveFunction(perturbed_parameter_array);
    error_array[2] = sum(objective_function_array);

    # Was the step good or bad?
    delta_error = error_array[2] - error_array[1];
  	if (delta_error<0)
  		LAMBDA_SUCC = 1;
  	else
  		LAMBDA_SUCC = 0;
  	end

  	# Update the step size -
  	(SIGMA,AVG_P_SUCC) = updateStepSize(LAMBDA_SUCC,AVG_P_SUCC,P_SUCC_TARGET,CP,D,SIGMA);

    # Ok, recompute the direction array -
    if (LAMBDA_SUCC == 1)

      # Message -
      msg = "Step improved the objective function. New error = "*string(error_array[2])*" Old error = "*string(error_array[1])*" counter = "*string(counter)
      println(msg)

  		# Keep the peturbed parameters -
  		parameter_array = perturbed_parameter_array;

  		# Update the covariance matrix -
  		A = updateCholesky(A,zV,CCOV,AVG_P_SUCC,PTHRESH);

  		# Update the error -
  		error_array[1] = error_array[2];

  		# If we have a successful step and we are stuck at the min, reset this
  		if (SIGMA<=1e-5)
  			SIGMA = 0.05;
  		end;

  		# Update the counter =
  		counter = counter + 1;

      # Write parameters to disk ...
      parameter_archive = [parameter_archive parameter_array]
      writedlm("./parameter_archive.dat.4",parameter_archive);

      error_archive = [error_archive error_array[1]]
      writedlm("./error_archive.dat.4",error_archive);

      # reset the number of failed steps -
      number_of_failed_steps = 0;

  	else

      # Message -
      msg = "No improved. New error = "*string(error_array[2])*" Old error = "*string(error_array[1])*" counter = "*string(counter)
      println(msg)

      # update the failed step counter -
  		number_of_failed_steps = number_of_failed_steps + 1;

  		# if we are beyond the max failed steps, reset the system w/a random parameter guess -
  		if (number_of_failed_steps>max_number_of_failed_steps)

  			# Reset the step length -
  			SIGMA = 0.10;

  			# Reset the direction array -
  			A = eye(N,N);
  			AVG_P_SUCC = P_SUCC_TARGET;

        # Generate new parameter set -
        (perturbed_parameter_array,zV) = sample_function(parameter_array,initial_parameter_array,N,A,SIGMA);

  			# Evaluate the error, that becomes the bext ...
        objective_function_array = pObjectiveFunction(perturbed_parameter_array);
        error_array[2] = sum(objective_function_array);

        # Reset -
        error_array[1] = error_array[2];

  			# let the user know that we are reseting -
  			msg = "Number of failed steps excedded - reset the parameters and search direction"
  			print(msg);

  			# reset the counter -
  			number_of_failed_steps = 0;
  		end

  	   # Check to see if we need to set FLAG -
  	  if (counter>=number_of_iterations || error_array[1] <= error_target)
  		    should_loop_continue = false;
  	  end
    end
  end
end
