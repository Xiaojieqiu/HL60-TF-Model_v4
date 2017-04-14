# include -
include("Include.jl")

const P3_color = "gray"
const P3_shaded_color="lightgray"

# load the P21 and E2F data -
data_P21 = readdlm("./plot_ensemble/scaled_data_P21.dat")
data_E2F = readdlm("./plot_ensemble/scaled_data_E2F.dat")


# find the mean, and std of the data sets -
mean_P21 = data_P21[:,1]
mean_E2F = data_E2F[:,2]

ind = np.arange(4)
width = 0.35
(fig, ax) = subplots()
ax[:bar](ind, mean_P21, width, color=P3_color)
ax[:bar](ind+width, mean_E2F, width, color=P3_shaded_color)
