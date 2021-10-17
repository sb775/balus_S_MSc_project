The folder contains:

1) The python script 'task4_automate_mutlist_mCSM-PPI2_v1.1'. This script generates the mutation lists for all 42 IGHV from Table 1 for mCSM-PPI2.

2) An R script called 'MSc_project_mutation_profiling_200921v2.R'. This contains: i) a function to plot all the stacked bar plots in Appendix, Figure B, ii) all statistical tests and iii)code to create all other plots in Figures 3-6.

3) The 3 files called 'PyIgClassify_and_final_list','ThreeLetterAACode' and 'Final_list.csv' are the input to run the python script. The folder called 'pdbseq&germline', also required is in the 'data_downloads' folder one level up.

4) An output folder called 'mscsm_ppi2_mutlists_auto' contains data output by mCSM-PPI2.

5) An output folder called 'mutabind2_mutlists' contains data output by MutaBind2.

6) An output folder called 'SAAMBE3D_mutlists' contains data output by SAAMBE3D.

7) An output folder called 'Webserver_1-69_results_compiled' contains:

	* Files which end in '_comp_results.csv' pertain to the combined results created from output from the 3 modellers for the 11 PDBs.  
	* Stacked barplots are saved as.bmp files according to the Physio-chemical property plotted.
	* File called 'df_all.csv' is an output file of a compiled dataframe of all 134 variants identified in this study.
	* File called 'final_list_of_interest.csv' is the output containing all 134 variants with only the columns of interest included.
	* File called 'final_list_top_21.xlsx' contains the top 21 variants of interest.
	* File called 'final_list_of_interest_sd_centile.xls' has the the multiples of sd and percentiles added on. 
