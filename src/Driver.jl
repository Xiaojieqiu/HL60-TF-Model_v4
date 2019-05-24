# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# include -
include("Include.jl")

# Script to solve the balance equations -
time_start = 0.0
time_stop = 120.0
time_step_size = 0.01

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the previous run -
par_array = readdlm("./raw_ensemble/parameter_archive.dat.23")
parameter_set = par_array[:,end]

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
