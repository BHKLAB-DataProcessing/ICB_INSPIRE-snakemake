library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")

#############################################################################
#############################################################################

load(file.path(input_dir, "INSPIRE_RNASEQ_TPM.RData"))

#############################################################################
#############################################################################
## Get Clinical data

clin = read.table( file.path(input_dir, "CLIN.txt") , sep="\t" , header=TRUE , stringsAsFactors = FALSE , dec=",")

clin = cbind( clin , "PD-1/PD-L1" , NA , NA , NA , NA , NA , NA )
colnames(clin) = c( "patient" , "age" , "primary" , "histo" , "recist" , "pfs" ,"os" , "t.pfs" , "t.os" , "drug_type" , "stage" , "sex" , "response" , "response.other.info" , "dna" , "rna" )
rownames(clin) = clin$patient

clin$pfs = ifelse( clin$pfs %in% "" , 0 , 1 )
clin$os = ifelse( clin$os %in% "" , 0 , 1 )

clin$response = Get_Response( data=clin )
clin$rna = "tpm"
clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

#############################################################################
#############################################################################

patient = intersect( colnames(expr) , rownames(clin) )
clin = clin[ patient , ]
expr =  expr[ , patient ]

# format patient ids by replacing any "-" with "_".
patient <- str_replace_all(patient, "-", "_")
rownames(clin) <- str_replace_all(rownames(clin), "-", "_")
clin$patient <- str_replace_all(clin$patient, "-", "_")
colnames(expr) <- str_replace_all(colnames(expr), "-", "_")

snv = as.data.frame( fread( file.path(output_dir, "SNV.csv") , sep=";" , stringsAsFactors=FALSE , header=TRUE , fill=TRUE ))
snv = snv[ snv$Sample %in% patient , ]
snv_patient = unique( sort( snv$Sample ) )

case = cbind( patient , 0 , 0 , 1 )
colnames(case ) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = case[,1]
case[ snv_patient , "snv" ] = 1

write.table( case , file = file.path(output_dir, "cased_sequenced.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( clin , file = file.path(output_dir, "CLIN.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( expr , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )


