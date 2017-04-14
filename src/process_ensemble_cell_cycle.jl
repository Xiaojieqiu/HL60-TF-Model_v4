include("Include.jl")

# Where are my files stored?
path_to_sim_files = "./simulation"

# Load the ensemble for P21 -
time_P21 = [0.0 24 48 72.0]
(raw_data_P21,scaled_data_P21) = process_ensemble_data(path_to_sim_files,time_P21,49);
writedlm("./plot_ensemble/scaled_data_P21.dat",scaled_data_P21);

# Load the ensemble for E2F -
time_E2F = [0.0 24 48 72.0]
(raw_data_E2F,scaled_data_E2F) = process_ensemble_data(path_to_sim_files,time_E2F,43);
writedlm("./plot_ensemble/scaled_data_E2F.dat",scaled_data_E2F);
