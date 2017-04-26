include("Include.jl")

# Where are my files stored?
path_to_sim_files = "./simulation"

# Load the ensemble for AhR -
time_AhR = [0.0 48.0]
(raw_data_AhR,scaled_data_AhR) = process_ensemble_data(path_to_sim_files,time_AhR,38);
writedlm("./plot_ensemble/scaled_data_AhR.dat",scaled_data_AhR);

# Load the ensemble for CD38 -
time_CD38 = [0.0 24.0 48.0 72.0]
(raw_data_CD38,scaled_data_CD38) = process_ensemble_data(path_to_sim_files,time_CD38,41);
writedlm("./plot_ensemble/scaled_data_CD38.dat",scaled_data_CD38);

# Load the ensemble for CD11b -
time_CD11b = [0.0 24.0 48.0 72.0]
(raw_data_CD11b,scaled_data_CD11b) = process_ensemble_data(path_to_sim_files,time_CD11b,39);
writedlm("./plot_ensemble/scaled_data_CD11b.dat",scaled_data_CD11b);

# Load the ensemble for P47Phox -
time_P47Phox = [0.0 48.0]
(raw_data_P47Phox,scaled_data_P47Phox) = process_ensemble_data(path_to_sim_files,time_P47Phox,50);
writedlm("./plot_ensemble/scaled_data_P47Phox.dat",scaled_data_P47Phox);

# Load the ensemble for EGR1 -
time_EGR1 = [0.0 48.0]
(raw_data_EGR1,scaled_data_EGR1) = process_ensemble_data(path_to_sim_files,time_EGR1,44);
writedlm("./plot_ensemble/scaled_data_EGR1.dat",scaled_data_EGR1);

# Load the ensemble for IRF1 -
time_IRF1 = [0.0 48.0]
(raw_data_IRF1,scaled_data_IRF1) = process_ensemble_data(path_to_sim_files,time_IRF1,46);
writedlm("./plot_ensemble/scaled_data_IRF1.dat",scaled_data_IRF1);

# Load the ensemble for GIF1 -
time_GIF1 = [0.0 48.0]
(raw_data_GIF1,scaled_data_GIF1) = process_ensemble_data(path_to_sim_files,time_GIF1,45);
writedlm("./plot_ensemble/scaled_data_GFI1.dat",scaled_data_GIF1);

# Load the ensemble for PU1 -
time_PU1 = [0.0 48.0]
(raw_data_PU1,scaled_data_PU1) = process_ensemble_data(path_to_sim_files,time_PU1,52);
writedlm("./plot_ensemble/scaled_data_PU1.dat",scaled_data_PU1);

# Load the ensemble for OCT4 -
time_OCT4 = [0.0 48.0]
(raw_data_OCT4,scaled_data_OCT4) = process_ensemble_data(path_to_sim_files,time_OCT4,48);
writedlm("./plot_ensemble/scaled_data_OCT4.dat",scaled_data_OCT4);

# Load the ensemble for CEBPa -
time_CEBPa = [0.0 48.0]
(raw_data_PU1,scaled_data_CEBPa) = process_ensemble_data(path_to_sim_files,time_CEBPa,42);
writedlm("./plot_ensemble/scaled_data_CEBPa.dat",scaled_data_CEBPa);
