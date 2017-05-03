# include -
include("Include.jl")

const simulation_color = "gray"
const measured_color="lightgray"

# load the CD38/CD11b data -
CD38_experimental_data = readdlm("./search_functions/data/CD38-WTHL60.txt")
CD11b_experimental_data = readdlm("./search_functions/data/CD11b-WTHL60.txt")

# load the simulation -
plot_index = 39
result_set = process_simulation_files("./simulation")
state_archive = []
time_archive = []
for (index,data_array) in enumerate(result_set)

  time_array = data_array[:,1];
  state_array = data_array[:,plot_index+1]*(1/data_array[1,plot_index+1]);

  if (index == 1)
    number_of_time_steps = length(time_array)
    state_archive = zeros(number_of_time_steps)
    time_archive = time_array
  end

  state_archive = [state_archive state_array]

end

state_archive = state_archive[:,2:end]
mean_state = mean(state_archive,2)
std_state = std(state_archive,2)
SF = 1.96
state_lower_bound = mean_state - SF*std_state
state_upper_bound = mean_state + SF*std_state


if (plot_index == 41)
  plot(time_archive,mean_state,lw=2,color="dodgerblue")
  fill_between(time_archive,vec(state_lower_bound),vec(state_upper_bound),color="lightskyblue",lw=3)
  errorbar(CD38_experimental_data[:,1],CD38_experimental_data[:,2],color="dodgerblue",yerr=CD38_experimental_data[:,3],fmt="o")
else
  plot(time_archive,mean_state,lw=2,color=simulation_color)
  fill_between(time_archive,vec(state_lower_bound),vec(state_upper_bound),color=measured_color,lw=3)
  errorbar(CD11b_experimental_data[:,1],CD11b_experimental_data[:,2],color="dimgray",yerr=CD11b_experimental_data[:,3],fmt="o")
end
