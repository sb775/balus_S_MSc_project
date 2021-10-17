### Aim:
## 1) Create pretty MSA


#if (!requireNamespace("BiocManager", quietly=TRUE))
 #  install.packages("BiocManager") 
#BiocManager::install("msa", force=TRUE)
library(msa)

## Example file from msa package
#inp_file = 'C:/Users/b/Desktop/Suchu/MSc_projects/Code/fasta/IGHV1-46_out.fasta'
#mySequences <- readAAStringSet(inp_file)
#ighv1_46 <- readAAStringSet(inp_file)

## Clustal W with default paramters
#myFirstAlignment <- msa(ighv1_46, "ClustalW")

##Create PDF file (help menu)
##"C:/Users/b/Documents/MSc_project_R"
#tmpFile <- tempfile(pattern="msa",tmpdir="C:/Users/b/Documents/MSc_project_R",fileext=".pdf")

#msaPrettyPrint(myFirstAlignment, file=tmpFile, output="pdf",showNames="left", 
 #              showLogo="top",showConsensus = "bottom",logoColors = "rasmol",
  #             verbose=TRUE,askForOverwrite=FALSE)


#### Above in a loop

#set the working directory from which the files will be read from
setwd('C:/Users/b/Desktop/Suchu/MSc_projects/Code/msa_pretty_aln')

#create a list of the files from target directory
input_path <- 'C:/Users/b/Desktop/Suchu/MSc_projects/Code/fasta/'
output_path <- 'C:/Users/b/Desktop/Suchu/MSc_projects/Code/msa_pretty_aln/'
file_list <- list.files(path = input_path, '.fasta')

for (i in 1:length(file_list)){
  inputs <- (file_list[i])
  filename <- paste0(input_path, inputs)
  allseqs <- readAAStringSet(filename)
  alignments <- msa(allseqs) ##default ClustalW
  output_name <- paste0(output_path,substr(inputs,1,nchar(inputs)-5),'pdf') 
  ##nchar counts length of string -5 to take of fasta
  
all_pdf <- msaPrettyPrint(alignments, file=output_name,output="pdf",showNames="left",
                          shadingMode="identical",shadingColors= "grays",showLogo="top",showConsensus = c("bottom"),
                          logoColors = "rasmol",showLogoScale=c("leftright"),verbose=TRUE,
                          askForOverwrite=FALSE)
}


