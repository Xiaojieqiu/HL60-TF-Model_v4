# include -
include("Include.jl")

const P3_color = "gray"
const P3_shaded_color="lightgray"

# load the experimental data -
experimental_data_array = readdlm("./search_functions/data/G1G0-WTHL60.txt")

# load the data_dictionary -
time_start = 0.0
time_stop = 80.0
time_step_size = 0.01
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# load the *best* fit set -
obj_array = readdlm("error_archive.dat.4")
par_array = readdlm("parameter_archive.dat.4")
obj_array = obj_array[2:end]
par_array = par_array[:,2:end]

# What is the index of the best fit set?
index_best_fit_set = indmin(obj_array)

# what parameters are we looking at?
parameter_set = par_array[:,index_best_fit_set]

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

# Create the Vandermonde matrix -

# create a synthetic time scale -
synthethic_time_scale = collect(0:0.1:72)
number_of_basis_rows = length(synthethic_time_scale)
number_of_basis_cols = 4

# interpolate the simulation results onto the experimental time scale -
inter_e2f = np.interp(synthethic_time_scale,time_array,simulation_state_array[:,43]);
inter_p21 = np.interp(synthethic_time_scale,time_array,simulation_state_array[:,49]);


# calculate the basis function entries -
basis_array = (inter_p21)./(inter_e2f)+(synthethic_time_scale./72.0)

# build the Vandermonde matrix
vandermonde_matrix = zeros(number_of_basis_rows,number_of_basis_cols)
for row_index = 1:number_of_basis_rows
  for col_index = 1:number_of_basis_cols

    if (col_index == 1)
      vandermonde_matrix[row_index,col_index] = 1.0
    else
      vandermonde_matrix[row_index,col_index] = basis_array[row_index]^(col_index-1)
    end
  end
end

# ok, calculate the coefficients?
vhat = inv(transpose(vandermonde_matrix)*vandermonde_matrix)*transpose(vandermonde_matrix)

# ok, we need to sample to get a set of possible curves -
number_of_samples = 10
yhat_archive = zeros(length(time_array),1)
for sample_index = 1:number_of_samples

  # generate a synthetic data set from the original -
  synthetic_data_set = zeros(4,2)
  synthetic_data_set[:,1] = experimental_data_array[:,1]

  # how many time samples do we have?
  number_of_time_samples = length(experimental_data_array[:,1])
  for time_sample_index = 1:number_of_time_samples

    # create a distribution -
    d = Normal(experimental_data_array[time_sample_index,2],experimental_data_array[time_sample_index,3])
    synthetic_data_set[time_sample_index,2] = rand(d)
  end

  # create the interpolated e-data
  inter_e_data = np.interp(synthethic_time_scale,synthetic_data_set[:,1],synthetic_data_set[:,2]);
  c_array = vhat*inter_e_data

  # make plot -
  full_basis_array = (simulation_state_array[:,49]./simulation_state_array[:,43])+(time_array./72.0)
  number_of_time_steps = length(time_array)
  y_hat_data_array = Array{Float64,1}()
  term_array = zeros(number_of_basis_cols)
  for time_index = 1:number_of_time_steps

      term_array[1] = c_array[1]

    for basis_index = 2:number_of_basis_cols

      term_array[basis_index] = c_array[basis_index]*(full_basis_array[time_index])^(basis_index-1)

    end
    total_value = sum(term_array)
    push!(y_hat_data_array,total_value)
  end

  yhat_archive = [yhat_archive y_hat_data_array]

end

yhat_archive = yhat_archive[:,2:end]
mean_yhat = mean(yhat_archive,2)
std_yhat = std(yhat_archive,2)
SF = 1.96
P3_lower_bound = mean_yhat - SF*std_yhat
P3_upper_bound = mean_yhat + SF*std_yhat
plot(time_array,mean_yhat,lw=2,color=P3_color)
fill_between(time_array,vec(P3_lower_bound),vec(P3_upper_bound),color=P3_shaded_color,lw=3)
errorbar(experimental_data_array[:,1],experimental_data_array[:,2],yerr=experimental_data_array[:,3],fmt="ko")
