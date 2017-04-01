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
# Function: calculate_transcription_rates
# Description: Calculate the transcriptional rate array at time t
# Generated on: 2017-03-09T06:52:55.843
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters
#
# Output arguments:
# transcription_rate_array::Array{Float64,1} => Transcriptional rate array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

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

	# Initialize the transcription rate -
	transcription_rate_array = zeros(18)
	KSAT = data_dictionary["saturation_constant_transcription"]
	kcat_transcription = data_dictionary["kcat_transcription"]
	rnapII_concentration = data_dictionary["rnapII_concentration"]
	average_transcript_length = data_dictionary["average_transcript_length"]
	gene_coding_length_array = data_dictionary["gene_coding_length_array"]

	# Populate the transcription rate array -
	# Gene: gene_AP1
	gene_length = gene_coding_length_array[1]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[1] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_AP1)/(KSAT+gene_AP1))

	# Gene: gene_AhR
	gene_length = gene_coding_length_array[2]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[2] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_AhR)/(KSAT+gene_AhR))

	# Gene: gene_CD11b
	gene_length = gene_coding_length_array[3]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[3] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_CD11b)/(KSAT+gene_CD11b))

	# Gene: gene_CD14
	gene_length = gene_coding_length_array[4]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[4] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_CD14)/(KSAT+gene_CD14))

	# Gene: gene_CD38
	gene_length = gene_coding_length_array[5]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[5] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_CD38)/(KSAT+gene_CD38))

	# Gene: gene_CEBPa
	gene_length = gene_coding_length_array[6]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[6] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_CEBPa)/(KSAT+gene_CEBPa))

	# Gene: gene_E2F
	gene_length = gene_coding_length_array[7]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[7] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_E2F)/(KSAT+gene_E2F))

	# Gene: gene_EGR1
	gene_length = gene_coding_length_array[8]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[8] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_EGR1)/(KSAT+gene_EGR1))

	# Gene: gene_GFI1
	gene_length = gene_coding_length_array[9]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[9] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_GFI1)/(KSAT+gene_GFI1))

	# Gene: gene_IRF1
	gene_length = gene_coding_length_array[10]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[10] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_IRF1)/(KSAT+gene_IRF1))

	# Gene: gene_OCT1
	gene_length = gene_coding_length_array[11]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[11] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_OCT1)/(KSAT+gene_OCT1))

	# Gene: gene_OCT4
	gene_length = gene_coding_length_array[12]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[12] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_OCT4)/(KSAT+gene_OCT4))

	# Gene: gene_P21
	gene_length = gene_coding_length_array[13]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[13] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_P21)/(KSAT+gene_P21))

	# Gene: gene_P47Phox
	gene_length = gene_coding_length_array[14]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[14] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_P47Phox)/(KSAT+gene_P47Phox))

	# Gene: gene_PPARg
	gene_length = gene_coding_length_array[15]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[15] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_PPARg)/(KSAT+gene_PPARg))

	# Gene: gene_PU1
	gene_length = gene_coding_length_array[16]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[16] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_PU1)/(KSAT+gene_PU1))

	# Gene: gene_Trigger
	gene_length = gene_coding_length_array[17]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[17] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_Trigger)/(KSAT+gene_Trigger))

	# Gene: gene_cRAF
	gene_length = gene_coding_length_array[18]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[18] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_cRAF)/(KSAT+gene_cRAF))


	# return transcription_rate_array -
	return transcription_rate_array
end

# ----------------------------------------------------------------------------------- #
# Function: calculate_background_transcription_rates
# Description: Calculate the leak transcriptional rate array at time t
# Generated on: 2017-03-09T06:52:55.843
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters
#
# Output arguments:
# background_transcription_rate_array::Array{Float64,1} => Background transcriptional rate array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_background_transcription_rates(t::Float64,x::Array{Float64,1},transcription_rate_array::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})
	return zeros(length(x))
end


# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_rates
# Description: Calculate the translation rate array at time t
# Generated on: 2017-03-09T06:52:55.843
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters
#
# Output arguments:
# translation_rate_array::Array{Float64,1} => Translation rate array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species -
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

	# Initialize the translation rate -
	translation_rate_array = zeros(18)
	KSAT = data_dictionary["saturation_constant_translation"]
	kcat_translation = data_dictionary["kcat_translation"]
	ribosome_concentration = data_dictionary["ribosome_concentration"]
	average_protein_length = data_dictionary["average_protein_length"]
	protein_coding_length_array = data_dictionary["protein_coding_length_array"]

	# Populate the translation rate array -
	# Transcript: mRNA_gene_AP1
	protein_length = protein_coding_length_array[1]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[1] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_AP1)/(KSAT+mRNA_gene_AP1))

	# Transcript: mRNA_gene_AhR
	protein_length = protein_coding_length_array[2]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[2] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_AhR)/(KSAT+mRNA_gene_AhR))

	# Transcript: mRNA_gene_CD11b
	protein_length = protein_coding_length_array[3]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[3] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_CD11b)/(KSAT+mRNA_gene_CD11b))

	# Transcript: mRNA_gene_CD14
	protein_length = protein_coding_length_array[4]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[4] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_CD14)/(KSAT+mRNA_gene_CD14))

	# Transcript: mRNA_gene_CD38
	protein_length = protein_coding_length_array[5]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[5] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_CD38)/(KSAT+mRNA_gene_CD38))

	# Transcript: mRNA_gene_CEBPa
	protein_length = protein_coding_length_array[6]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[6] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_CEBPa)/(KSAT+mRNA_gene_CEBPa))

	# Transcript: mRNA_gene_E2F
	protein_length = protein_coding_length_array[7]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[7] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_E2F)/(KSAT+mRNA_gene_E2F))

	# Transcript: mRNA_gene_EGR1
	protein_length = protein_coding_length_array[8]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[8] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_EGR1)/(KSAT+mRNA_gene_EGR1))

	# Transcript: mRNA_gene_GFI1
	protein_length = protein_coding_length_array[9]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[9] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_GFI1)/(KSAT+mRNA_gene_GFI1))

	# Transcript: mRNA_gene_IRF1
	protein_length = protein_coding_length_array[10]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[10] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_IRF1)/(KSAT+mRNA_gene_IRF1))

	# Transcript: mRNA_gene_OCT1
	protein_length = protein_coding_length_array[11]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[11] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_OCT1)/(KSAT+mRNA_gene_OCT1))

	# Transcript: mRNA_gene_OCT4
	protein_length = protein_coding_length_array[12]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[12] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_OCT4)/(KSAT+mRNA_gene_OCT4))

	# Transcript: mRNA_gene_P21
	protein_length = protein_coding_length_array[13]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[13] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_P21)/(KSAT+mRNA_gene_P21))

	# Transcript: mRNA_gene_P47Phox
	protein_length = protein_coding_length_array[14]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[14] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_P47Phox)/(KSAT+mRNA_gene_P47Phox))

	# Transcript: mRNA_gene_PPARg
	protein_length = protein_coding_length_array[15]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[15] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_PPARg)/(KSAT+mRNA_gene_PPARg))

	# Transcript: mRNA_gene_PU1
	protein_length = protein_coding_length_array[16]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[16] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_PU1)/(KSAT+mRNA_gene_PU1))

	# Transcript: mRNA_gene_Trigger
	protein_length = protein_coding_length_array[17]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[17] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_Trigger)/(KSAT+mRNA_gene_Trigger))

	# Transcript: mRNA_gene_cRAF
	protein_length = protein_coding_length_array[18]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[18] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_cRAF)/(KSAT+mRNA_gene_cRAF))


	# return translation array -
	return translation_rate_array
end

# ----------------------------------------------------------------------------------- #
# Function: calculate_mRNA_degradation_rates
# Description: Calculate the mRNA degradation rate array at time t
# Generated on: 2017-03-09T06:52:55.845
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters
#
# Output arguments:
# mRNA_degradation_rate_array::Array{Float64,1} => mRNA degradation rate array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_mRNA_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species -
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

	# Initialize the degrdation array -
	degradation_rate_array = zeros(18)
	mRNA_degrdation_constant = data_dictionary["degradation_constant_mRNA"]

	# get my corrections -
	# beta_mRNA_gene_AP1 = data_dictionary["beta_mRNA_gene_AP1"]
	# beta_mRNA_gene_AhR = data_dictionary["beta_mRNA_gene_AhR"]
	# beta_mRNA_gene_CD11b = data_dictionary["beta_mRNA_gene_CD11b"]
	# beta_mRNA_gene_CD14 = data_dictionary["beta_mRNA_gene_CD14"]
	# beta_mRNA_gene_CD38 = data_dictionary["beta_mRNA_gene_CD38"]
	# beta_mRNA_gene_CEBPa = data_dictionary["beta_mRNA_gene_CEBPa"]
	# beta_mRNA_gene_E2F = data_dictionary["beta_mRNA_gene_E2F"]
	# beta_mRNA_gene_EGR1 = data_dictionary["beta_mRNA_gene_EGR1"]
	# beta_mRNA_gene_GFI1 = data_dictionary["beta_mRNA_gene_GFI1"]
	# beta_mRNA_gene_IRF1 = data_dictionary["beta_mRNA_gene_IRF1"]
	# beta_mRNA_gene_OCT1 = data_dictionary["beta_mRNA_gene_OCT1"]
	# beta_mRNA_gene_OCT4 = data_dictionary["beta_mRNA_gene_OCT4"]
	# beta_mRNA_gene_P21 = data_dictionary["beta_mRNA_gene_P21"]
	# beta_mRNA_gene_P47Phox = data_dictionary["beta_mRNA_gene_P47Phox"]
	# beta_mRNA_gene_PPARg = data_dictionary["beta_mRNA_gene_PPARg"]
	# beta_mRNA_gene_PU1 = data_dictionary["beta_mRNA_gene_PU1"]
	# beta_mRNA_gene_Trigger = data_dictionary["beta_mRNA_gene_Trigger"]
	# beta_mRNA_gene_cRAF = data_dictionary["beta_mRNA_gene_cRAF"]

	# Calculate the degradation_rate_array -
	degradation_rate_array[1] = (mRNA_degrdation_constant)*mRNA_gene_AP1
	degradation_rate_array[2] = (mRNA_degrdation_constant)*mRNA_gene_AhR
	degradation_rate_array[3] = (mRNA_degrdation_constant)*mRNA_gene_CD11b
	degradation_rate_array[4] = (mRNA_degrdation_constant)*mRNA_gene_CD14
	degradation_rate_array[5] = (mRNA_degrdation_constant)*mRNA_gene_CD38
	degradation_rate_array[6] = (mRNA_degrdation_constant)*mRNA_gene_CEBPa
	degradation_rate_array[7] = (mRNA_degrdation_constant)*mRNA_gene_E2F
	degradation_rate_array[8] = (mRNA_degrdation_constant)*mRNA_gene_EGR1
	degradation_rate_array[9] = (mRNA_degrdation_constant)*mRNA_gene_GFI1
	degradation_rate_array[10] = (mRNA_degrdation_constant)*mRNA_gene_IRF1
	degradation_rate_array[11] = (mRNA_degrdation_constant)*mRNA_gene_OCT1
	degradation_rate_array[12] = (mRNA_degrdation_constant)*mRNA_gene_OCT4
	degradation_rate_array[13] = (mRNA_degrdation_constant)*mRNA_gene_P21
	degradation_rate_array[14] = (mRNA_degrdation_constant)*mRNA_gene_P47Phox
	degradation_rate_array[15] = (mRNA_degrdation_constant)*mRNA_gene_PPARg
	degradation_rate_array[16] = (mRNA_degrdation_constant)*mRNA_gene_PU1
	degradation_rate_array[17] = (mRNA_degrdation_constant)*mRNA_gene_Trigger
	degradation_rate_array[18] = (mRNA_degrdation_constant)*mRNA_gene_cRAF

	# return the degrdation rate array -
	return degradation_rate_array
end

# ----------------------------------------------------------------------------------- #
# Function: calculate_protein_degradation_rates
# Description: Calculate the protein degradation rate array at time t
# Generated on: 2017-03-09T06:52:55.847
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters
#
# Output arguments:
# protein_degradation_rate_array::Array{Float64,1} => protein degradation rate array (number_of_proteins x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_protein_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species -
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

	# Initialize the degrdation array -
	degradation_rate_array = zeros(18)
	protein_degrdation_constant = data_dictionary["degradation_constant_protein"]

	# Calculate the degradation_rate_array -
	degradation_rate_array[1] = (protein_degrdation_constant)*protein_gene_AP1
	degradation_rate_array[2] = (protein_degrdation_constant)*protein_gene_AhR
	degradation_rate_array[3] = (protein_degrdation_constant)*protein_gene_CD11b
	degradation_rate_array[4] = (protein_degrdation_constant)*protein_gene_CD14
	degradation_rate_array[5] = (protein_degrdation_constant)*protein_gene_CD38
	degradation_rate_array[6] = (protein_degrdation_constant)*protein_gene_CEBPa
	degradation_rate_array[7] = (protein_degrdation_constant)*protein_gene_E2F
	degradation_rate_array[8] = (protein_degrdation_constant)*protein_gene_EGR1
	degradation_rate_array[9] = (protein_degrdation_constant)*protein_gene_GFI1
	degradation_rate_array[10] = (protein_degrdation_constant)*protein_gene_IRF1
	degradation_rate_array[11] = (protein_degrdation_constant)*protein_gene_OCT1
	degradation_rate_array[12] = (protein_degrdation_constant)*protein_gene_OCT4
	degradation_rate_array[13] = (protein_degrdation_constant)*protein_gene_P21
	degradation_rate_array[14] = (protein_degrdation_constant)*protein_gene_P47Phox
	degradation_rate_array[15] = (protein_degrdation_constant)*protein_gene_PPARg
	degradation_rate_array[16] = (protein_degrdation_constant)*protein_gene_PU1
	degradation_rate_array[17] = (protein_degrdation_constant)*protein_gene_Trigger
	degradation_rate_array[18] = (protein_degrdation_constant)*protein_gene_cRAF

	# return the degrdation rate array -
	return degradation_rate_array
end
