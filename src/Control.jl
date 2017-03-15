# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
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
#
# ----------------------------------------------------------------------------------- #
# Function: Control
# Description: Calculate the transcriptional control array at time t
# Generated on: 2017-03-09T06:52:55.849
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# initialize the control - 
	control_array = zeros(18)

	# Alias the species - 
	gene_AP1 = x[1]
	gene_AhR = x[2]
	gene_CD11b = x[3]
	gene_CD14 = x[4]
	gene_CD38 = x[5]
	gene_CEBPa = x[6]
	gene_E2F = x[7]
	gene_EGR1 = x[8]
	gene_GFI1 = x[9]
	gene_IRF1 = x[10]
	gene_OCT1 = x[11]
	gene_OCT4 = x[12]
	gene_P21 = x[13]
	gene_P47Phox = x[14]
	gene_PPARg = x[15]
	gene_PU1 = x[16]
	gene_Trigger = x[17]
	gene_cRAF = x[18]
	mRNA_gene_AP1 = x[19]
	mRNA_gene_AhR = x[20]
	mRNA_gene_CD11b = x[21]
	mRNA_gene_CD14 = x[22]
	mRNA_gene_CD38 = x[23]
	mRNA_gene_CEBPa = x[24]
	mRNA_gene_E2F = x[25]
	mRNA_gene_EGR1 = x[26]
	mRNA_gene_GFI1 = x[27]
	mRNA_gene_IRF1 = x[28]
	mRNA_gene_OCT1 = x[29]
	mRNA_gene_OCT4 = x[30]
	mRNA_gene_P21 = x[31]
	mRNA_gene_P47Phox = x[32]
	mRNA_gene_PPARg = x[33]
	mRNA_gene_PU1 = x[34]
	mRNA_gene_Trigger = x[35]
	mRNA_gene_cRAF = x[36]
	protein_gene_AP1 = x[37]
	protein_gene_AhR = x[38]
	protein_gene_CD11b = x[39]
	protein_gene_CD14 = x[40]
	protein_gene_CD38 = x[41]
	protein_gene_CEBPa = x[42]
	protein_gene_E2F = x[43]
	protein_gene_EGR1 = x[44]
	protein_gene_GFI1 = x[45]
	protein_gene_IRF1 = x[46]
	protein_gene_OCT1 = x[47]
	protein_gene_OCT4 = x[48]
	protein_gene_P21 = x[49]
	protein_gene_P47Phox = x[50]
	protein_gene_PPARg = x[51]
	protein_gene_PU1 = x[52]
	protein_gene_Trigger = x[53]
	protein_gene_cRAF = x[54]

	# Alias the binding parameters - 
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_gene_AP1_gene_AhR = binding_parameter_dictionary["n_gene_AP1_gene_AhR"]
	K_gene_AP1_gene_AhR = binding_parameter_dictionary["K_gene_AP1_gene_AhR"]
	n_gene_AP1_gene_PU1 = binding_parameter_dictionary["n_gene_AP1_gene_PU1"]
	K_gene_AP1_gene_PU1 = binding_parameter_dictionary["K_gene_AP1_gene_PU1"]
	n_gene_AP1_gene_PPARg = binding_parameter_dictionary["n_gene_AP1_gene_PPARg"]
	K_gene_AP1_gene_PPARg = binding_parameter_dictionary["K_gene_AP1_gene_PPARg"]
	n_gene_AhR_gene_Trigger = binding_parameter_dictionary["n_gene_AhR_gene_Trigger"]
	K_gene_AhR_gene_Trigger = binding_parameter_dictionary["K_gene_AhR_gene_Trigger"]
	n_gene_CD11b_gene_PU1_gene_cRAF = binding_parameter_dictionary["n_gene_CD11b_gene_PU1_gene_cRAF"]
	K_gene_CD11b_gene_PU1_gene_cRAF = binding_parameter_dictionary["K_gene_CD11b_gene_PU1_gene_cRAF"]
	n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF = binding_parameter_dictionary["n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"]
	K_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF = binding_parameter_dictionary["K_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"]
	n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF = binding_parameter_dictionary["n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"]
	K_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF = binding_parameter_dictionary["K_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"]
	n_gene_CEBPa_gene_Trigger = binding_parameter_dictionary["n_gene_CEBPa_gene_Trigger"]
	K_gene_CEBPa_gene_Trigger = binding_parameter_dictionary["K_gene_CEBPa_gene_Trigger"]
	n_gene_CEBPa_gene_PPARg = binding_parameter_dictionary["n_gene_CEBPa_gene_PPARg"]
	K_gene_CEBPa_gene_PPARg = binding_parameter_dictionary["K_gene_CEBPa_gene_PPARg"]
	n_gene_CEBPa_gene_CEBPa = binding_parameter_dictionary["n_gene_CEBPa_gene_CEBPa"]
	K_gene_CEBPa_gene_CEBPa = binding_parameter_dictionary["K_gene_CEBPa_gene_CEBPa"]
	n_gene_CEBPa_gene_GFI1 = binding_parameter_dictionary["n_gene_CEBPa_gene_GFI1"]
	K_gene_CEBPa_gene_GFI1 = binding_parameter_dictionary["K_gene_CEBPa_gene_GFI1"]
	n_gene_E2F_gene_E2F = binding_parameter_dictionary["n_gene_E2F_gene_E2F"]
	K_gene_E2F_gene_E2F = binding_parameter_dictionary["K_gene_E2F_gene_E2F"]
	n_gene_E2F_gene_PPARg = binding_parameter_dictionary["n_gene_E2F_gene_PPARg"]
	K_gene_E2F_gene_PPARg = binding_parameter_dictionary["K_gene_E2F_gene_PPARg"]
	n_gene_E2F_gene_CEBPa = binding_parameter_dictionary["n_gene_E2F_gene_CEBPa"]
	K_gene_E2F_gene_CEBPa = binding_parameter_dictionary["K_gene_E2F_gene_CEBPa"]
	n_gene_E2F_gene_GFI1 = binding_parameter_dictionary["n_gene_E2F_gene_GFI1"]
	K_gene_E2F_gene_GFI1 = binding_parameter_dictionary["K_gene_E2F_gene_GFI1"]
	n_gene_E2F_gene_cRAF = binding_parameter_dictionary["n_gene_E2F_gene_cRAF"]
	K_gene_E2F_gene_cRAF = binding_parameter_dictionary["K_gene_E2F_gene_cRAF"]
	n_gene_EGR1_gene_Trigger = binding_parameter_dictionary["n_gene_EGR1_gene_Trigger"]
	K_gene_EGR1_gene_Trigger = binding_parameter_dictionary["K_gene_EGR1_gene_Trigger"]
	n_gene_EGR1_gene_PU1 = binding_parameter_dictionary["n_gene_EGR1_gene_PU1"]
	K_gene_EGR1_gene_PU1 = binding_parameter_dictionary["K_gene_EGR1_gene_PU1"]
	n_gene_EGR1_gene_PPARg = binding_parameter_dictionary["n_gene_EGR1_gene_PPARg"]
	K_gene_EGR1_gene_PPARg = binding_parameter_dictionary["K_gene_EGR1_gene_PPARg"]
	n_gene_EGR1_gene_GFI1 = binding_parameter_dictionary["n_gene_EGR1_gene_GFI1"]
	K_gene_EGR1_gene_GFI1 = binding_parameter_dictionary["K_gene_EGR1_gene_GFI1"]
	n_gene_GFI1_gene_CEBPa = binding_parameter_dictionary["n_gene_GFI1_gene_CEBPa"]
	K_gene_GFI1_gene_CEBPa = binding_parameter_dictionary["K_gene_GFI1_gene_CEBPa"]
	n_gene_GFI1_gene_EGR1 = binding_parameter_dictionary["n_gene_GFI1_gene_EGR1"]
	K_gene_GFI1_gene_EGR1 = binding_parameter_dictionary["K_gene_GFI1_gene_EGR1"]
	n_gene_IRF1_gene_Trigger = binding_parameter_dictionary["n_gene_IRF1_gene_Trigger"]
	K_gene_IRF1_gene_Trigger = binding_parameter_dictionary["K_gene_IRF1_gene_Trigger"]
	n_gene_IRF1_gene_AhR = binding_parameter_dictionary["n_gene_IRF1_gene_AhR"]
	K_gene_IRF1_gene_AhR = binding_parameter_dictionary["K_gene_IRF1_gene_AhR"]
	n_gene_IRF1_gene_PPARg = binding_parameter_dictionary["n_gene_IRF1_gene_PPARg"]
	K_gene_IRF1_gene_PPARg = binding_parameter_dictionary["K_gene_IRF1_gene_PPARg"]
	n_gene_OCT1_gene_PPARg = binding_parameter_dictionary["n_gene_OCT1_gene_PPARg"]
	K_gene_OCT1_gene_PPARg = binding_parameter_dictionary["K_gene_OCT1_gene_PPARg"]
	n_gene_OCT4_gene_Trigger = binding_parameter_dictionary["n_gene_OCT4_gene_Trigger"]
	K_gene_OCT4_gene_Trigger = binding_parameter_dictionary["K_gene_OCT4_gene_Trigger"]
	n_gene_OCT4_gene_AhR = binding_parameter_dictionary["n_gene_OCT4_gene_AhR"]
	K_gene_OCT4_gene_AhR = binding_parameter_dictionary["K_gene_OCT4_gene_AhR"]
	n_gene_OCT4_gene_cRAF = binding_parameter_dictionary["n_gene_OCT4_gene_cRAF"]
	K_gene_OCT4_gene_cRAF = binding_parameter_dictionary["K_gene_OCT4_gene_cRAF"]
	n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF = binding_parameter_dictionary["n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"]
	K_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF = binding_parameter_dictionary["K_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"]
	n_gene_P21_gene_GFI1 = binding_parameter_dictionary["n_gene_P21_gene_GFI1"]
	K_gene_P21_gene_GFI1 = binding_parameter_dictionary["K_gene_P21_gene_GFI1"]
	n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF = binding_parameter_dictionary["n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"]
	K_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF = binding_parameter_dictionary["K_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"]
	n_gene_P47Phox_gene_PPARg = binding_parameter_dictionary["n_gene_P47Phox_gene_PPARg"]
	K_gene_P47Phox_gene_PPARg = binding_parameter_dictionary["K_gene_P47Phox_gene_PPARg"]
	n_gene_PPARg_gene_Trigger = binding_parameter_dictionary["n_gene_PPARg_gene_Trigger"]
	K_gene_PPARg_gene_Trigger = binding_parameter_dictionary["K_gene_PPARg_gene_Trigger"]
	n_gene_PPARg_gene_CEBPa = binding_parameter_dictionary["n_gene_PPARg_gene_CEBPa"]
	K_gene_PPARg_gene_CEBPa = binding_parameter_dictionary["K_gene_PPARg_gene_CEBPa"]
	n_gene_PPARg_gene_EGR1 = binding_parameter_dictionary["n_gene_PPARg_gene_EGR1"]
	K_gene_PPARg_gene_EGR1 = binding_parameter_dictionary["K_gene_PPARg_gene_EGR1"]
	n_gene_PPARg_gene_PU1 = binding_parameter_dictionary["n_gene_PPARg_gene_PU1"]
	K_gene_PPARg_gene_PU1 = binding_parameter_dictionary["K_gene_PPARg_gene_PU1"]
	n_gene_PPARg_gene_AP1 = binding_parameter_dictionary["n_gene_PPARg_gene_AP1"]
	K_gene_PPARg_gene_AP1 = binding_parameter_dictionary["K_gene_PPARg_gene_AP1"]
	n_gene_PU1_gene_Trigger = binding_parameter_dictionary["n_gene_PU1_gene_Trigger"]
	K_gene_PU1_gene_Trigger = binding_parameter_dictionary["K_gene_PU1_gene_Trigger"]
	n_gene_PU1_gene_CEBPa = binding_parameter_dictionary["n_gene_PU1_gene_CEBPa"]
	K_gene_PU1_gene_CEBPa = binding_parameter_dictionary["K_gene_PU1_gene_CEBPa"]
	n_gene_PU1_gene_PU1 = binding_parameter_dictionary["n_gene_PU1_gene_PU1"]
	K_gene_PU1_gene_PU1 = binding_parameter_dictionary["K_gene_PU1_gene_PU1"]
	n_gene_PU1_gene_AP1 = binding_parameter_dictionary["n_gene_PU1_gene_AP1"]
	K_gene_PU1_gene_AP1 = binding_parameter_dictionary["K_gene_PU1_gene_AP1"]
	n_gene_PU1_gene_OCT1 = binding_parameter_dictionary["n_gene_PU1_gene_OCT1"]
	K_gene_PU1_gene_OCT1 = binding_parameter_dictionary["K_gene_PU1_gene_OCT1"]
	n_gene_PU1_gene_AhR = binding_parameter_dictionary["n_gene_PU1_gene_AhR"]
	K_gene_PU1_gene_AhR = binding_parameter_dictionary["K_gene_PU1_gene_AhR"]
	n_gene_PU1_gene_GFI1 = binding_parameter_dictionary["n_gene_PU1_gene_GFI1"]
	K_gene_PU1_gene_GFI1 = binding_parameter_dictionary["K_gene_PU1_gene_GFI1"]

	# Alias the control function parameters - 
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_gene_AP1_RNAP = control_parameter_dictionary["W_gene_AP1_RNAP"]
	W_gene_AP1_gene_AhR = control_parameter_dictionary["W_gene_AP1_gene_AhR"]
	W_gene_AP1_gene_PU1 = control_parameter_dictionary["W_gene_AP1_gene_PU1"]
	W_gene_AP1_gene_PPARg = control_parameter_dictionary["W_gene_AP1_gene_PPARg"]
	W_gene_AhR_RNAP = control_parameter_dictionary["W_gene_AhR_RNAP"]
	W_gene_AhR_gene_Trigger = control_parameter_dictionary["W_gene_AhR_gene_Trigger"]
	W_gene_CD11b_RNAP = control_parameter_dictionary["W_gene_CD11b_RNAP"]
	W_gene_CD11b_gene_PU1_gene_cRAF = control_parameter_dictionary["W_gene_CD11b_gene_PU1_gene_cRAF"]
	W_gene_CD14_RNAP = control_parameter_dictionary["W_gene_CD14_RNAP"]
	W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF = control_parameter_dictionary["W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"]
	W_gene_CD38_RNAP = control_parameter_dictionary["W_gene_CD38_RNAP"]
	W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF = control_parameter_dictionary["W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"]
	W_gene_CEBPa_RNAP = control_parameter_dictionary["W_gene_CEBPa_RNAP"]
	W_gene_CEBPa_gene_Trigger = control_parameter_dictionary["W_gene_CEBPa_gene_Trigger"]
	W_gene_CEBPa_gene_PPARg = control_parameter_dictionary["W_gene_CEBPa_gene_PPARg"]
	W_gene_CEBPa_gene_CEBPa = control_parameter_dictionary["W_gene_CEBPa_gene_CEBPa"]
	W_gene_CEBPa_gene_GFI1 = control_parameter_dictionary["W_gene_CEBPa_gene_GFI1"]
	W_gene_E2F_RNAP = control_parameter_dictionary["W_gene_E2F_RNAP"]
	W_gene_E2F_gene_E2F = control_parameter_dictionary["W_gene_E2F_gene_E2F"]
	W_gene_E2F_gene_PPARg = control_parameter_dictionary["W_gene_E2F_gene_PPARg"]
	W_gene_E2F_gene_CEBPa = control_parameter_dictionary["W_gene_E2F_gene_CEBPa"]
	W_gene_E2F_gene_GFI1 = control_parameter_dictionary["W_gene_E2F_gene_GFI1"]
	W_gene_E2F_gene_cRAF = control_parameter_dictionary["W_gene_E2F_gene_cRAF"]
	W_gene_EGR1_RNAP = control_parameter_dictionary["W_gene_EGR1_RNAP"]
	W_gene_EGR1_gene_Trigger = control_parameter_dictionary["W_gene_EGR1_gene_Trigger"]
	W_gene_EGR1_gene_PU1 = control_parameter_dictionary["W_gene_EGR1_gene_PU1"]
	W_gene_EGR1_gene_PPARg = control_parameter_dictionary["W_gene_EGR1_gene_PPARg"]
	W_gene_EGR1_gene_GFI1 = control_parameter_dictionary["W_gene_EGR1_gene_GFI1"]
	W_gene_GFI1_RNAP = control_parameter_dictionary["W_gene_GFI1_RNAP"]
	W_gene_GFI1_gene_CEBPa = control_parameter_dictionary["W_gene_GFI1_gene_CEBPa"]
	W_gene_GFI1_gene_EGR1 = control_parameter_dictionary["W_gene_GFI1_gene_EGR1"]
	W_gene_IRF1_RNAP = control_parameter_dictionary["W_gene_IRF1_RNAP"]
	W_gene_IRF1_gene_Trigger = control_parameter_dictionary["W_gene_IRF1_gene_Trigger"]
	W_gene_IRF1_gene_AhR = control_parameter_dictionary["W_gene_IRF1_gene_AhR"]
	W_gene_IRF1_gene_PPARg = control_parameter_dictionary["W_gene_IRF1_gene_PPARg"]
	W_gene_OCT1_RNAP = control_parameter_dictionary["W_gene_OCT1_RNAP"]
	W_gene_OCT1_gene_PPARg = control_parameter_dictionary["W_gene_OCT1_gene_PPARg"]
	W_gene_OCT4_RNAP = control_parameter_dictionary["W_gene_OCT4_RNAP"]
	W_gene_OCT4_gene_Trigger = control_parameter_dictionary["W_gene_OCT4_gene_Trigger"]
	W_gene_OCT4_gene_AhR = control_parameter_dictionary["W_gene_OCT4_gene_AhR"]
	W_gene_OCT4_gene_cRAF = control_parameter_dictionary["W_gene_OCT4_gene_cRAF"]
	W_gene_P21_RNAP = control_parameter_dictionary["W_gene_P21_RNAP"]
	W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF = control_parameter_dictionary["W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"]
	W_gene_P21_gene_GFI1 = control_parameter_dictionary["W_gene_P21_gene_GFI1"]
	W_gene_P47Phox_RNAP = control_parameter_dictionary["W_gene_P47Phox_RNAP"]
	W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF = control_parameter_dictionary["W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"]
	W_gene_P47Phox_gene_PPARg = control_parameter_dictionary["W_gene_P47Phox_gene_PPARg"]
	W_gene_PPARg_RNAP = control_parameter_dictionary["W_gene_PPARg_RNAP"]
	W_gene_PPARg_gene_Trigger = control_parameter_dictionary["W_gene_PPARg_gene_Trigger"]
	W_gene_PPARg_gene_CEBPa = control_parameter_dictionary["W_gene_PPARg_gene_CEBPa"]
	W_gene_PPARg_gene_EGR1 = control_parameter_dictionary["W_gene_PPARg_gene_EGR1"]
	W_gene_PPARg_gene_PU1 = control_parameter_dictionary["W_gene_PPARg_gene_PU1"]
	W_gene_PPARg_gene_AP1 = control_parameter_dictionary["W_gene_PPARg_gene_AP1"]
	W_gene_PU1_RNAP = control_parameter_dictionary["W_gene_PU1_RNAP"]
	W_gene_PU1_gene_Trigger = control_parameter_dictionary["W_gene_PU1_gene_Trigger"]
	W_gene_PU1_gene_CEBPa = control_parameter_dictionary["W_gene_PU1_gene_CEBPa"]
	W_gene_PU1_gene_PU1 = control_parameter_dictionary["W_gene_PU1_gene_PU1"]
	W_gene_PU1_gene_AP1 = control_parameter_dictionary["W_gene_PU1_gene_AP1"]
	W_gene_PU1_gene_OCT1 = control_parameter_dictionary["W_gene_PU1_gene_OCT1"]
	W_gene_PU1_gene_AhR = control_parameter_dictionary["W_gene_PU1_gene_AhR"]
	W_gene_PU1_gene_GFI1 = control_parameter_dictionary["W_gene_PU1_gene_GFI1"]
	W_gene_Trigger_RNAP = control_parameter_dictionary["W_gene_Trigger_RNAP"]
	W_gene_cRAF_RNAP = control_parameter_dictionary["W_gene_cRAF_RNAP"]

	# Transfer function target:gene_AP1 actor:gene_AhR
	actor_set_gene_AP1_gene_AhR = [
		protein_gene_AhR
	]
	actor = prod(actor_set_gene_AP1_gene_AhR)
	b_gene_AP1_gene_AhR = (actor^(n_gene_AP1_gene_AhR))/(K_gene_AP1_gene_AhR^(n_gene_AP1_gene_AhR)+actor^(n_gene_AP1_gene_AhR))

	# Transfer function target:gene_AP1 actor:gene_PU1
	actor_set_gene_AP1_gene_PU1 = [
		protein_gene_PU1
	]
	actor = prod(actor_set_gene_AP1_gene_PU1)
	b_gene_AP1_gene_PU1 = (actor^(n_gene_AP1_gene_PU1))/(K_gene_AP1_gene_PU1^(n_gene_AP1_gene_PU1)+actor^(n_gene_AP1_gene_PU1))

	# Transfer function target:gene_AP1 actor:gene_PPARg
	actor_set_gene_AP1_gene_PPARg = [
		protein_gene_PPARg
	]
	actor = prod(actor_set_gene_AP1_gene_PPARg)
	b_gene_AP1_gene_PPARg = (actor^(n_gene_AP1_gene_PPARg))/(K_gene_AP1_gene_PPARg^(n_gene_AP1_gene_PPARg)+actor^(n_gene_AP1_gene_PPARg))

	# Control function for gene_AP1 - 
	control_array[1] = (W_gene_AP1_RNAP+W_gene_AP1_gene_AhR*b_gene_AP1_gene_AhR+W_gene_AP1_gene_PU1*b_gene_AP1_gene_PU1)/(1+W_gene_AP1_RNAP+W_gene_AP1_gene_AhR*b_gene_AP1_gene_AhR+W_gene_AP1_gene_PU1*b_gene_AP1_gene_PU1+W_gene_AP1_gene_PPARg*b_gene_AP1_gene_PPARg)

	# Transfer function target:gene_AhR actor:gene_Trigger
	actor_set_gene_AhR_gene_Trigger = [
		protein_gene_Trigger
	]
	actor = prod(actor_set_gene_AhR_gene_Trigger)
	b_gene_AhR_gene_Trigger = (actor^(n_gene_AhR_gene_Trigger))/(K_gene_AhR_gene_Trigger^(n_gene_AhR_gene_Trigger)+actor^(n_gene_AhR_gene_Trigger))

	# Control function for gene_AhR - 
	control_array[2] = (W_gene_AhR_RNAP+W_gene_AhR_gene_Trigger*b_gene_AhR_gene_Trigger)/(1+W_gene_AhR_RNAP+W_gene_AhR_gene_Trigger*b_gene_AhR_gene_Trigger)

	# Transfer function target:gene_CD11b actor:gene_PU1_gene_cRAF
	actor_set_gene_CD11b_gene_PU1_gene_cRAF = [
		protein_gene_PU1
		protein_gene_cRAF
	]
	actor = prod(actor_set_gene_CD11b_gene_PU1_gene_cRAF)
	b_gene_CD11b_gene_PU1_gene_cRAF = (actor^(n_gene_CD11b_gene_PU1_gene_cRAF))/(K_gene_CD11b_gene_PU1_gene_cRAF^(n_gene_CD11b_gene_PU1_gene_cRAF)+actor^(n_gene_CD11b_gene_PU1_gene_cRAF))

	# Control function for gene_CD11b - 
	control_array[3] = (W_gene_CD11b_RNAP+W_gene_CD11b_gene_PU1_gene_cRAF*b_gene_CD11b_gene_PU1_gene_cRAF)/(1+W_gene_CD11b_RNAP+W_gene_CD11b_gene_PU1_gene_cRAF*b_gene_CD11b_gene_PU1_gene_cRAF)

	# Transfer function target:gene_CD14 actor:gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF
	actor_set_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF = [
		protein_gene_PPARg
		protein_gene_CEBPa
		protein_gene_EGR1
		protein_gene_cRAF
	]
	actor = prod(actor_set_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF)
	b_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF = (actor^(n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF))/(K_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF^(n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF)+actor^(n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF))

	# Control function for gene_CD14 - 
	control_array[4] = (W_gene_CD14_RNAP+W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF*b_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF)/(1+W_gene_CD14_RNAP+W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF*b_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF)

	# Transfer function target:gene_CD38 actor:gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF
	actor_set_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF = [
		protein_gene_IRF1
		protein_gene_PPARg
		protein_gene_Trigger
		protein_gene_cRAF
	]
	actor = prod(actor_set_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF)
	b_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF = (actor^(n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF))/(K_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF^(n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF)+actor^(n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF))

	# Control function for gene_CD38 - 
	control_array[5] = (W_gene_CD38_RNAP+W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF*b_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF)/(1+W_gene_CD38_RNAP+W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF*b_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF)

	# Transfer function target:gene_CEBPa actor:gene_Trigger
	actor_set_gene_CEBPa_gene_Trigger = [
		protein_gene_Trigger
	]
	actor = prod(actor_set_gene_CEBPa_gene_Trigger)
	b_gene_CEBPa_gene_Trigger = (actor^(n_gene_CEBPa_gene_Trigger))/(K_gene_CEBPa_gene_Trigger^(n_gene_CEBPa_gene_Trigger)+actor^(n_gene_CEBPa_gene_Trigger))

	# Transfer function target:gene_CEBPa actor:gene_PPARg
	actor_set_gene_CEBPa_gene_PPARg = [
		protein_gene_PPARg
	]
	actor = prod(actor_set_gene_CEBPa_gene_PPARg)
	b_gene_CEBPa_gene_PPARg = (actor^(n_gene_CEBPa_gene_PPARg))/(K_gene_CEBPa_gene_PPARg^(n_gene_CEBPa_gene_PPARg)+actor^(n_gene_CEBPa_gene_PPARg))

	# Transfer function target:gene_CEBPa actor:gene_CEBPa
	actor_set_gene_CEBPa_gene_CEBPa = [
		protein_gene_CEBPa
	]
	actor = prod(actor_set_gene_CEBPa_gene_CEBPa)
	b_gene_CEBPa_gene_CEBPa = (actor^(n_gene_CEBPa_gene_CEBPa))/(K_gene_CEBPa_gene_CEBPa^(n_gene_CEBPa_gene_CEBPa)+actor^(n_gene_CEBPa_gene_CEBPa))

	# Transfer function target:gene_CEBPa actor:gene_GFI1
	actor_set_gene_CEBPa_gene_GFI1 = [
		protein_gene_GFI1
	]
	actor = prod(actor_set_gene_CEBPa_gene_GFI1)
	b_gene_CEBPa_gene_GFI1 = (actor^(n_gene_CEBPa_gene_GFI1))/(K_gene_CEBPa_gene_GFI1^(n_gene_CEBPa_gene_GFI1)+actor^(n_gene_CEBPa_gene_GFI1))

	# Control function for gene_CEBPa - 
	control_array[6] = (W_gene_CEBPa_RNAP+W_gene_CEBPa_gene_Trigger*b_gene_CEBPa_gene_Trigger+W_gene_CEBPa_gene_PPARg*b_gene_CEBPa_gene_PPARg+W_gene_CEBPa_gene_CEBPa*b_gene_CEBPa_gene_CEBPa)/(1+W_gene_CEBPa_RNAP+W_gene_CEBPa_gene_Trigger*b_gene_CEBPa_gene_Trigger+W_gene_CEBPa_gene_PPARg*b_gene_CEBPa_gene_PPARg+W_gene_CEBPa_gene_CEBPa*b_gene_CEBPa_gene_CEBPa+W_gene_CEBPa_gene_GFI1*b_gene_CEBPa_gene_GFI1)

	# Transfer function target:gene_E2F actor:gene_E2F
	actor_set_gene_E2F_gene_E2F = [
		protein_gene_E2F
	]
	actor = prod(actor_set_gene_E2F_gene_E2F)
	b_gene_E2F_gene_E2F = (actor^(n_gene_E2F_gene_E2F))/(K_gene_E2F_gene_E2F^(n_gene_E2F_gene_E2F)+actor^(n_gene_E2F_gene_E2F))

	# Transfer function target:gene_E2F actor:gene_PPARg
	actor_set_gene_E2F_gene_PPARg = [
		protein_gene_PPARg
	]
	actor = prod(actor_set_gene_E2F_gene_PPARg)
	b_gene_E2F_gene_PPARg = (actor^(n_gene_E2F_gene_PPARg))/(K_gene_E2F_gene_PPARg^(n_gene_E2F_gene_PPARg)+actor^(n_gene_E2F_gene_PPARg))

	# Transfer function target:gene_E2F actor:gene_CEBPa
	actor_set_gene_E2F_gene_CEBPa = [
		protein_gene_CEBPa
	]
	actor = prod(actor_set_gene_E2F_gene_CEBPa)
	b_gene_E2F_gene_CEBPa = (actor^(n_gene_E2F_gene_CEBPa))/(K_gene_E2F_gene_CEBPa^(n_gene_E2F_gene_CEBPa)+actor^(n_gene_E2F_gene_CEBPa))

	# Transfer function target:gene_E2F actor:gene_GFI1
	actor_set_gene_E2F_gene_GFI1 = [
		protein_gene_GFI1
	]
	actor = prod(actor_set_gene_E2F_gene_GFI1)
	b_gene_E2F_gene_GFI1 = (actor^(n_gene_E2F_gene_GFI1))/(K_gene_E2F_gene_GFI1^(n_gene_E2F_gene_GFI1)+actor^(n_gene_E2F_gene_GFI1))

	# Transfer function target:gene_E2F actor:gene_cRAF
	actor_set_gene_E2F_gene_cRAF = [
		protein_gene_cRAF
	]
	actor = prod(actor_set_gene_E2F_gene_cRAF)
	b_gene_E2F_gene_cRAF = (actor^(n_gene_E2F_gene_cRAF))/(K_gene_E2F_gene_cRAF^(n_gene_E2F_gene_cRAF)+actor^(n_gene_E2F_gene_cRAF))

	# Control function for gene_E2F - 
	control_array[7] = (W_gene_E2F_RNAP+W_gene_E2F_gene_E2F*b_gene_E2F_gene_E2F)/(1+W_gene_E2F_RNAP+W_gene_E2F_gene_E2F*b_gene_E2F_gene_E2F+W_gene_E2F_gene_PPARg*b_gene_E2F_gene_PPARg+W_gene_E2F_gene_CEBPa*b_gene_E2F_gene_CEBPa+W_gene_E2F_gene_GFI1*b_gene_E2F_gene_GFI1+W_gene_E2F_gene_cRAF*b_gene_E2F_gene_cRAF)

	# Transfer function target:gene_EGR1 actor:gene_Trigger
	actor_set_gene_EGR1_gene_Trigger = [
		protein_gene_Trigger
	]
	actor = prod(actor_set_gene_EGR1_gene_Trigger)
	b_gene_EGR1_gene_Trigger = (actor^(n_gene_EGR1_gene_Trigger))/(K_gene_EGR1_gene_Trigger^(n_gene_EGR1_gene_Trigger)+actor^(n_gene_EGR1_gene_Trigger))

	# Transfer function target:gene_EGR1 actor:gene_PU1
	actor_set_gene_EGR1_gene_PU1 = [
		protein_gene_PU1
	]
	actor = prod(actor_set_gene_EGR1_gene_PU1)
	b_gene_EGR1_gene_PU1 = (actor^(n_gene_EGR1_gene_PU1))/(K_gene_EGR1_gene_PU1^(n_gene_EGR1_gene_PU1)+actor^(n_gene_EGR1_gene_PU1))

	# Transfer function target:gene_EGR1 actor:gene_PPARg
	actor_set_gene_EGR1_gene_PPARg = [
		protein_gene_PPARg
	]
	actor = prod(actor_set_gene_EGR1_gene_PPARg)
	b_gene_EGR1_gene_PPARg = (actor^(n_gene_EGR1_gene_PPARg))/(K_gene_EGR1_gene_PPARg^(n_gene_EGR1_gene_PPARg)+actor^(n_gene_EGR1_gene_PPARg))

	# Transfer function target:gene_EGR1 actor:gene_GFI1
	actor_set_gene_EGR1_gene_GFI1 = [
		protein_gene_GFI1
	]
	actor = prod(actor_set_gene_EGR1_gene_GFI1)
	b_gene_EGR1_gene_GFI1 = (actor^(n_gene_EGR1_gene_GFI1))/(K_gene_EGR1_gene_GFI1^(n_gene_EGR1_gene_GFI1)+actor^(n_gene_EGR1_gene_GFI1))

	# Control function for gene_EGR1 - 
	control_array[8] = (W_gene_EGR1_RNAP+W_gene_EGR1_gene_Trigger*b_gene_EGR1_gene_Trigger+W_gene_EGR1_gene_PU1*b_gene_EGR1_gene_PU1)/(1+W_gene_EGR1_RNAP+W_gene_EGR1_gene_Trigger*b_gene_EGR1_gene_Trigger+W_gene_EGR1_gene_PU1*b_gene_EGR1_gene_PU1+W_gene_EGR1_gene_PPARg*b_gene_EGR1_gene_PPARg+W_gene_EGR1_gene_GFI1*b_gene_EGR1_gene_GFI1)

	# Transfer function target:gene_GFI1 actor:gene_CEBPa
	actor_set_gene_GFI1_gene_CEBPa = [
		protein_gene_CEBPa
	]
	actor = prod(actor_set_gene_GFI1_gene_CEBPa)
	b_gene_GFI1_gene_CEBPa = (actor^(n_gene_GFI1_gene_CEBPa))/(K_gene_GFI1_gene_CEBPa^(n_gene_GFI1_gene_CEBPa)+actor^(n_gene_GFI1_gene_CEBPa))

	# Transfer function target:gene_GFI1 actor:gene_EGR1
	actor_set_gene_GFI1_gene_EGR1 = [
		protein_gene_EGR1
	]
	actor = prod(actor_set_gene_GFI1_gene_EGR1)
	b_gene_GFI1_gene_EGR1 = (actor^(n_gene_GFI1_gene_EGR1))/(K_gene_GFI1_gene_EGR1^(n_gene_GFI1_gene_EGR1)+actor^(n_gene_GFI1_gene_EGR1))

	# Control function for gene_GFI1 - 
	control_array[9] = (W_gene_GFI1_RNAP+W_gene_GFI1_gene_CEBPa*b_gene_GFI1_gene_CEBPa)/(1+W_gene_GFI1_RNAP+W_gene_GFI1_gene_CEBPa*b_gene_GFI1_gene_CEBPa+W_gene_GFI1_gene_EGR1*b_gene_GFI1_gene_EGR1)

	# Transfer function target:gene_IRF1 actor:gene_Trigger
	actor_set_gene_IRF1_gene_Trigger = [
		protein_gene_Trigger
	]
	actor = prod(actor_set_gene_IRF1_gene_Trigger)
	b_gene_IRF1_gene_Trigger = (actor^(n_gene_IRF1_gene_Trigger))/(K_gene_IRF1_gene_Trigger^(n_gene_IRF1_gene_Trigger)+actor^(n_gene_IRF1_gene_Trigger))

	# Transfer function target:gene_IRF1 actor:gene_AhR
	actor_set_gene_IRF1_gene_AhR = [
		protein_gene_AhR
	]
	actor = prod(actor_set_gene_IRF1_gene_AhR)
	b_gene_IRF1_gene_AhR = (actor^(n_gene_IRF1_gene_AhR))/(K_gene_IRF1_gene_AhR^(n_gene_IRF1_gene_AhR)+actor^(n_gene_IRF1_gene_AhR))

	# Transfer function target:gene_IRF1 actor:gene_PPARg
	actor_set_gene_IRF1_gene_PPARg = [
		protein_gene_PPARg
	]
	actor = prod(actor_set_gene_IRF1_gene_PPARg)
	b_gene_IRF1_gene_PPARg = (actor^(n_gene_IRF1_gene_PPARg))/(K_gene_IRF1_gene_PPARg^(n_gene_IRF1_gene_PPARg)+actor^(n_gene_IRF1_gene_PPARg))

	# Control function for gene_IRF1 - 
	control_array[10] = (W_gene_IRF1_RNAP+W_gene_IRF1_gene_Trigger*b_gene_IRF1_gene_Trigger+W_gene_IRF1_gene_AhR*b_gene_IRF1_gene_AhR+W_gene_IRF1_gene_PPARg*b_gene_IRF1_gene_PPARg)/(1+W_gene_IRF1_RNAP+W_gene_IRF1_gene_Trigger*b_gene_IRF1_gene_Trigger+W_gene_IRF1_gene_AhR*b_gene_IRF1_gene_AhR+W_gene_IRF1_gene_PPARg*b_gene_IRF1_gene_PPARg)

	# Transfer function target:gene_OCT1 actor:gene_PPARg
	actor_set_gene_OCT1_gene_PPARg = [
		protein_gene_PPARg
	]
	actor = prod(actor_set_gene_OCT1_gene_PPARg)
	b_gene_OCT1_gene_PPARg = (actor^(n_gene_OCT1_gene_PPARg))/(K_gene_OCT1_gene_PPARg^(n_gene_OCT1_gene_PPARg)+actor^(n_gene_OCT1_gene_PPARg))

	# Control function for gene_OCT1 - 
	control_array[11] = (W_gene_OCT1_RNAP+W_gene_OCT1_gene_PPARg*b_gene_OCT1_gene_PPARg)/(1+W_gene_OCT1_RNAP+W_gene_OCT1_gene_PPARg*b_gene_OCT1_gene_PPARg)

	# Transfer function target:gene_OCT4 actor:gene_Trigger
	actor_set_gene_OCT4_gene_Trigger = [
		protein_gene_Trigger
	]
	actor = prod(actor_set_gene_OCT4_gene_Trigger)
	b_gene_OCT4_gene_Trigger = (actor^(n_gene_OCT4_gene_Trigger))/(K_gene_OCT4_gene_Trigger^(n_gene_OCT4_gene_Trigger)+actor^(n_gene_OCT4_gene_Trigger))

	# Transfer function target:gene_OCT4 actor:gene_AhR
	actor_set_gene_OCT4_gene_AhR = [
		protein_gene_AhR
	]
	actor = prod(actor_set_gene_OCT4_gene_AhR)
	b_gene_OCT4_gene_AhR = (actor^(n_gene_OCT4_gene_AhR))/(K_gene_OCT4_gene_AhR^(n_gene_OCT4_gene_AhR)+actor^(n_gene_OCT4_gene_AhR))

	# Transfer function target:gene_OCT4 actor:gene_cRAF
	actor_set_gene_OCT4_gene_cRAF = [
		protein_gene_cRAF
	]
	actor = prod(actor_set_gene_OCT4_gene_cRAF)
	b_gene_OCT4_gene_cRAF = (actor^(n_gene_OCT4_gene_cRAF))/(K_gene_OCT4_gene_cRAF^(n_gene_OCT4_gene_cRAF)+actor^(n_gene_OCT4_gene_cRAF))

	# Control function for gene_OCT4 - 
	control_array[12] = (W_gene_OCT4_RNAP)/(1+W_gene_OCT4_RNAP+W_gene_OCT4_gene_Trigger*b_gene_OCT4_gene_Trigger+W_gene_OCT4_gene_AhR*b_gene_OCT4_gene_AhR+W_gene_OCT4_gene_cRAF*b_gene_OCT4_gene_cRAF)

	# Transfer function target:gene_P21 actor:gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF
	actor_set_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF = [
		protein_gene_Trigger
		protein_gene_AP1
		protein_gene_PPARg
		protein_gene_PU1
		protein_gene_IRF1
		protein_gene_CEBPa
		protein_gene_cRAF
	]
	actor = prod(actor_set_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF)
	b_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF = (actor^(n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF))/(K_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF^(n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF)+actor^(n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF))

	# Transfer function target:gene_P21 actor:gene_GFI1
	actor_set_gene_P21_gene_GFI1 = [
		protein_gene_GFI1
	]
	actor = prod(actor_set_gene_P21_gene_GFI1)
	b_gene_P21_gene_GFI1 = (actor^(n_gene_P21_gene_GFI1))/(K_gene_P21_gene_GFI1^(n_gene_P21_gene_GFI1)+actor^(n_gene_P21_gene_GFI1))

	# Control function for gene_P21 - 
	control_array[13] = (W_gene_P21_RNAP+W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF*b_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF)/(1+W_gene_P21_RNAP+W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF*b_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF+W_gene_P21_gene_GFI1*b_gene_P21_gene_GFI1)

	# Transfer function target:gene_P47Phox actor:gene_PU1_gene_CEBPa_gene_cRAF
	actor_set_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF = [
		protein_gene_PU1
		protein_gene_CEBPa
		protein_gene_cRAF
	]
	actor = prod(actor_set_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF)
	b_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF = (actor^(n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF))/(K_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF^(n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF)+actor^(n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF))

	# Transfer function target:gene_P47Phox actor:gene_PPARg
	actor_set_gene_P47Phox_gene_PPARg = [
		protein_gene_PPARg
	]
	actor = prod(actor_set_gene_P47Phox_gene_PPARg)
	b_gene_P47Phox_gene_PPARg = (actor^(n_gene_P47Phox_gene_PPARg))/(K_gene_P47Phox_gene_PPARg^(n_gene_P47Phox_gene_PPARg)+actor^(n_gene_P47Phox_gene_PPARg))

	# Control function for gene_P47Phox - 
	control_array[14] = (W_gene_P47Phox_RNAP+W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF*b_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF)/(1+W_gene_P47Phox_RNAP+W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF*b_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF+W_gene_P47Phox_gene_PPARg*b_gene_P47Phox_gene_PPARg)

	# Transfer function target:gene_PPARg actor:gene_Trigger
	actor_set_gene_PPARg_gene_Trigger = [
		protein_gene_Trigger
	]
	actor = prod(actor_set_gene_PPARg_gene_Trigger)
	b_gene_PPARg_gene_Trigger = (actor^(n_gene_PPARg_gene_Trigger))/(K_gene_PPARg_gene_Trigger^(n_gene_PPARg_gene_Trigger)+actor^(n_gene_PPARg_gene_Trigger))

	# Transfer function target:gene_PPARg actor:gene_CEBPa
	actor_set_gene_PPARg_gene_CEBPa = [
		protein_gene_CEBPa
	]
	actor = prod(actor_set_gene_PPARg_gene_CEBPa)
	b_gene_PPARg_gene_CEBPa = (actor^(n_gene_PPARg_gene_CEBPa))/(K_gene_PPARg_gene_CEBPa^(n_gene_PPARg_gene_CEBPa)+actor^(n_gene_PPARg_gene_CEBPa))

	# Transfer function target:gene_PPARg actor:gene_EGR1
	actor_set_gene_PPARg_gene_EGR1 = [
		protein_gene_EGR1
	]
	actor = prod(actor_set_gene_PPARg_gene_EGR1)
	b_gene_PPARg_gene_EGR1 = (actor^(n_gene_PPARg_gene_EGR1))/(K_gene_PPARg_gene_EGR1^(n_gene_PPARg_gene_EGR1)+actor^(n_gene_PPARg_gene_EGR1))

	# Transfer function target:gene_PPARg actor:gene_PU1
	actor_set_gene_PPARg_gene_PU1 = [
		protein_gene_PU1
	]
	actor = prod(actor_set_gene_PPARg_gene_PU1)
	b_gene_PPARg_gene_PU1 = (actor^(n_gene_PPARg_gene_PU1))/(K_gene_PPARg_gene_PU1^(n_gene_PPARg_gene_PU1)+actor^(n_gene_PPARg_gene_PU1))

	# Transfer function target:gene_PPARg actor:gene_AP1
	actor_set_gene_PPARg_gene_AP1 = [
		protein_gene_AP1
	]
	actor = prod(actor_set_gene_PPARg_gene_AP1)
	b_gene_PPARg_gene_AP1 = (actor^(n_gene_PPARg_gene_AP1))/(K_gene_PPARg_gene_AP1^(n_gene_PPARg_gene_AP1)+actor^(n_gene_PPARg_gene_AP1))

	# Control function for gene_PPARg - 
	control_array[15] = (W_gene_PPARg_RNAP+W_gene_PPARg_gene_Trigger*b_gene_PPARg_gene_Trigger+W_gene_PPARg_gene_CEBPa*b_gene_PPARg_gene_CEBPa+W_gene_PPARg_gene_EGR1*b_gene_PPARg_gene_EGR1)/(1+W_gene_PPARg_RNAP+W_gene_PPARg_gene_Trigger*b_gene_PPARg_gene_Trigger+W_gene_PPARg_gene_CEBPa*b_gene_PPARg_gene_CEBPa+W_gene_PPARg_gene_EGR1*b_gene_PPARg_gene_EGR1+W_gene_PPARg_gene_PU1*b_gene_PPARg_gene_PU1+W_gene_PPARg_gene_AP1*b_gene_PPARg_gene_AP1)

	# Transfer function target:gene_PU1 actor:gene_Trigger
	actor_set_gene_PU1_gene_Trigger = [
		protein_gene_Trigger
	]
	actor = prod(actor_set_gene_PU1_gene_Trigger)
	b_gene_PU1_gene_Trigger = (actor^(n_gene_PU1_gene_Trigger))/(K_gene_PU1_gene_Trigger^(n_gene_PU1_gene_Trigger)+actor^(n_gene_PU1_gene_Trigger))

	# Transfer function target:gene_PU1 actor:gene_CEBPa
	actor_set_gene_PU1_gene_CEBPa = [
		protein_gene_CEBPa
	]
	actor = prod(actor_set_gene_PU1_gene_CEBPa)
	b_gene_PU1_gene_CEBPa = (actor^(n_gene_PU1_gene_CEBPa))/(K_gene_PU1_gene_CEBPa^(n_gene_PU1_gene_CEBPa)+actor^(n_gene_PU1_gene_CEBPa))

	# Transfer function target:gene_PU1 actor:gene_PU1
	actor_set_gene_PU1_gene_PU1 = [
		protein_gene_PU1
	]
	actor = prod(actor_set_gene_PU1_gene_PU1)
	b_gene_PU1_gene_PU1 = (actor^(n_gene_PU1_gene_PU1))/(K_gene_PU1_gene_PU1^(n_gene_PU1_gene_PU1)+actor^(n_gene_PU1_gene_PU1))

	# Transfer function target:gene_PU1 actor:gene_AP1
	actor_set_gene_PU1_gene_AP1 = [
		protein_gene_AP1
	]
	actor = prod(actor_set_gene_PU1_gene_AP1)
	b_gene_PU1_gene_AP1 = (actor^(n_gene_PU1_gene_AP1))/(K_gene_PU1_gene_AP1^(n_gene_PU1_gene_AP1)+actor^(n_gene_PU1_gene_AP1))

	# Transfer function target:gene_PU1 actor:gene_OCT1
	actor_set_gene_PU1_gene_OCT1 = [
		protein_gene_OCT1
	]
	actor = prod(actor_set_gene_PU1_gene_OCT1)
	b_gene_PU1_gene_OCT1 = (actor^(n_gene_PU1_gene_OCT1))/(K_gene_PU1_gene_OCT1^(n_gene_PU1_gene_OCT1)+actor^(n_gene_PU1_gene_OCT1))

	# Transfer function target:gene_PU1 actor:gene_AhR
	actor_set_gene_PU1_gene_AhR = [
		protein_gene_AhR
	]
	actor = prod(actor_set_gene_PU1_gene_AhR)
	b_gene_PU1_gene_AhR = (actor^(n_gene_PU1_gene_AhR))/(K_gene_PU1_gene_AhR^(n_gene_PU1_gene_AhR)+actor^(n_gene_PU1_gene_AhR))

	# Transfer function target:gene_PU1 actor:gene_GFI1
	actor_set_gene_PU1_gene_GFI1 = [
		protein_gene_GFI1
	]
	actor = prod(actor_set_gene_PU1_gene_GFI1)
	b_gene_PU1_gene_GFI1 = (actor^(n_gene_PU1_gene_GFI1))/(K_gene_PU1_gene_GFI1^(n_gene_PU1_gene_GFI1)+actor^(n_gene_PU1_gene_GFI1))

	# Control function for gene_PU1 - 
	control_array[16] = (W_gene_PU1_RNAP+W_gene_PU1_gene_Trigger*b_gene_PU1_gene_Trigger+W_gene_PU1_gene_CEBPa*b_gene_PU1_gene_CEBPa+W_gene_PU1_gene_PU1*b_gene_PU1_gene_PU1+W_gene_PU1_gene_AP1*b_gene_PU1_gene_AP1+W_gene_PU1_gene_OCT1*b_gene_PU1_gene_OCT1)/(1+W_gene_PU1_RNAP+W_gene_PU1_gene_Trigger*b_gene_PU1_gene_Trigger+W_gene_PU1_gene_CEBPa*b_gene_PU1_gene_CEBPa+W_gene_PU1_gene_PU1*b_gene_PU1_gene_PU1+W_gene_PU1_gene_AP1*b_gene_PU1_gene_AP1+W_gene_PU1_gene_OCT1*b_gene_PU1_gene_OCT1+W_gene_PU1_gene_AhR*b_gene_PU1_gene_AhR+W_gene_PU1_gene_GFI1*b_gene_PU1_gene_GFI1)

	# Control function for gene_Trigger - 
	control_array[17] = (W_gene_Trigger_RNAP)/(1+W_gene_Trigger_RNAP)

	# Control function for gene_cRAF - 
	control_array[18] = (W_gene_cRAF_RNAP)/(1+W_gene_cRAF_RNAP)

	# return - 
	return control_array
end
