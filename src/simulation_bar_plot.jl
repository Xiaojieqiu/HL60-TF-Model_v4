# include -
include("Include.jl")

const simulation_color = "gray"
const measured_color="lightgray"

# What species are we going to plot?
species_name_array = [
  "AhR",
  "IRF1",
  "P47Phox",
  "OCT4",
  "PU1",
  "EGR1",
  "CEBPa",
  "GFI1"
]

# load the experimental, and simulation data arrays -
simulation_data_dictionary = Dict{String,Any}()
experimental_data_dictionary = Dict{String,Any}()
for species_name in species_name_array

  # load the simulation data -
  simulation_path = "./plot_ensemble/scaled_data_$(species_name).dat"
  tmp_array = readdlm(simulation_path)
  tmp_mean = mean(tmp_array,2)
  tmp_std = std(tmp_array,2)

  # how many time points do we have?
  number_of_timepoints = length(tmp_mean)
  simulation_data_array = zeros(number_of_timepoints,2)
  for time_index = 1:number_of_timepoints
    simulation_data_array[time_index,1] = tmp_mean[time_index]
    simulation_data_array[time_index,2] = tmp_std[time_index]
  end
  simulation_data_dictionary[species_name] = simulation_data_array

  # load the experimental data -
  experimental_path = "./search_functions/data/$(species_name)-WTHL60.txt"
  experimental_data_dictionary[species_name] = readdlm(experimental_path)

end

# make the simulation data array -
simulation_data_array = zeros(length(species_name_array),2)
local_index = 1
for species_name in species_name_array

  simulation_data_array[local_index,1] = simulation_data_dictionary[species_name][2,1]
  simulation_data_array[local_index,2] = (1/sqrt(3))*simulation_data_dictionary[species_name][2,2]
  local_index = local_index + 1
end

# make simulation bars -
error_bar_dictionary = PyDict()
error_bar_dictionary["ecolor"]="black"
error_bar_dictionary["lw"]=1
error_bar_dictionary["capsize"]=5
error_bar_dictionary["capthick"]=1
ind = np.arange(length(species_name_array))
width = 0.35
(fig, ax) = subplots()
ax[:bar](ind, simulation_data_array[:,1], width, color=simulation_color,yerr=simulation_data_array[:,2],error_kw=error_bar_dictionary)
ax[:margins](0.05)

# make the experiment_data_array
experiment_data_array = zeros(length(species_name_array),2)
local_index = 1
for species_name in species_name_array

  experiment_data_array[local_index,1] = experimental_data_dictionary[species_name][2,2]
  experiment_data_array[local_index,2] = (1/sqrt(3))*experimental_data_dictionary[species_name][2,3]
  local_index = local_index + 1
end
ax[:bar](ind+width, experiment_data_array[:,1], width, color=measured_color,yerr=experiment_data_array[:,2],error_kw=error_bar_dictionary)

# write -
savefig("./raw_figs/SimBar-Raw.pdf")
