# include -
include("Include.jl")

# Setup the simulation timescale -
time_start = 0.0
time_stop = 120.0
time_step_size = 0.01

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the previous run -
# initial_parameter_guess = vec(readdlm("parameter_best_v1.dat.1"))
parameter_archive = readdlm("parameter_archive.dat.1")
objective_archive = readdlm("error_archive.dat.1")

# strip -
parameter_archive = parameter_archive[:,2:end]
objective_archive = objective_archive[:,2:end]

# find min -
index_min = indmin(objective_archive)

# grab the initial_parameter_guess -
initial_parameter_guess = parameter_archive[:,index_min]

# Setup stoping criteria -
NUMBER_OF_ITERATIONS = 100;
ERR_TARGET = 0.001;

# Setup function pointers -
pObjectiveFunction = objective_function;

# estimate parameters -
estimate_parameters(pObjectiveFunction,initial_parameter_guess,data_dictionary,NUMBER_OF_ITERATIONS,ERR_TARGET);
