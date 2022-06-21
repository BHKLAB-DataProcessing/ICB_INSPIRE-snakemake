library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

get_Info = function( mut ){
	mutType = sapply( mut , function(x){ unlist( strsplit( x , "|" , fixed = TRUE ) )[2] } )
	gene = sapply( mut , function(x){ unlist( strsplit( x , "|" , fixed = TRUE ) )[4] } )
	names( mutType ) = names( gene ) = NA
	
	info = cbind( mutType , gene )

	info
}

get_Sample = function( mut ){
	mutType = sapply( mut , function(x){ unlist( strsplit( x , "|" , fixed = TRUE ) )[2] } )
	gene = sapply( mut , function(x){ unlist( strsplit( x , "|" , fixed = TRUE ) )[4] } )
	names( mutType ) = names( gene ) = NA
	
	info = cbind( mutType , gene )

	info
}

################################################
################################################

snv = as.data.frame( fread( file.path(input_dir, "INSPIRE.combinedVariants.snv.released.112020.renamed.vcf.gz") , sep="\t" , stringsAsFactors=FALSE , skip = 490 ))

strelka = snv[ , grep( "STRELKA" , colnames(snv) ) ]
sample = table(sapply(colnames(strelka) , function(x){unlist(strsplit(x,".",fixed=T))[1]}) )
strelka = strelka[ , !colnames(strelka) %in% paste( names(sample[ sample == 2 ] ) , "STRELKA" , sep = "" ) ]

annot = snv[,9]

variant = c( "^AU$" , "^CU$" , "^GU$" , "^TU$")
names(variant) = c( "A" , "C" , "G" , "T" )

snv_data = NULL
for( i in 1:nrow(strelka) ){

	print( paste( i , "in" , nrow(strelka) , '(' , round( i/nrow(strelka) * 100 ) , '%)' ) )

	mu = strelka[ i , ]

	dp = grep( "^DP$", unlist( strsplit( annot[i] , ":" , fixed = TRUE )) )
	alt = grep( variant[snv[i,5]], unlist( strsplit( annot[i] , ":" , fixed = TRUE )) )
	
	dp = as.numeric( as.character( sapply( mu , function(x){ unlist( strsplit( x , ":" , fixed = TRUE ) )[dp] } ) ) )
	alt = as.numeric( as.character( sapply( mu , function(x){ unlist( strsplit(  unlist( strsplit( x , ":" , fixed = TRUE ) )[alt] , "," , fixed = "TRUE" ))[1] } ) ) )

	vaf = alt/dp *100

	mu = mu [ !is.na(vaf) & dp >= 20 & vaf > 5 ]

	if( length( mu ) > 0 ){
		snv_data = rbind( snv_data , 
							cbind( snv[ i , c( "#CHROM" , "POS" , "REF" , "ALT" ) ] , get_Info( mut = snv[ i , 8 ] ) , "SNV" , sapply( names(mu) , function(x){ unlist( strsplit( x , "." , fixed = TRUE ))[1] } ) ) )
	}
}

################################################
################################################


indel = as.data.frame( fread( file.path(input_dir, "INSPIRE.combinedVariants.indel.released.112020.renamed.vcf.gz") , sep="\t" , stringsAsFactors=FALSE , skip = 490 ))

strelka = indel[ , grep( "STRELKA" , colnames(indel) ) ]
sample = table(sapply(colnames(strelka) , function(x){unlist(strsplit(x,".",fixed=T))[1]}) )
strelka = strelka[ , !colnames(strelka) %in% paste( names(sample[ sample == 2 ] ) , "STRELKA" , sep = "" ) ]

indel_data = NULL
for( i in 1:nrow(strelka) ){

	print( paste( i , "in" , nrow(strelka) , '(' , round( i/nrow(strelka) * 100 ) , '%)' ) )

	mu = strelka[ i , ]
	
	ref = as.numeric( as.character( sapply( mu , function(x){ unlist( strsplit(  unlist( strsplit( x , ":" , fixed = TRUE ) )[5] , "," , fixed = "TRUE" ))[1] } ) ) )
	alt = as.numeric( as.character( sapply( mu , function(x){ unlist( strsplit(  unlist( strsplit( x , ":" , fixed = TRUE ) )[2] , "," , fixed = "TRUE" ))[1] } ) ) )

	vaf = alt/(ref+alt) *100

	mu = mu [ !is.na(vaf) & (ref+alt) >= 20 & vaf > 5 ]

	if( length( mu ) > 0 ){
		indel_data = rbind( indel_data , 
							cbind( indel[ i , c( "#CHROM" , "POS" , "REF" , "ALT" ) ] , get_Info( mut = indel[ i , 8 ] ) , "INDEL" , sapply( names(mu) , function(x){ unlist( strsplit( x , "." , fixed = TRUE ))[1] } ) ) )
	}
}

################################################
################################################
colnames( snv_data ) = colnames( indel_data ) = c( "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "Gene" , "MutType" , "Sample" )
data = rbind( snv_data , indel_data )

data$Sample = sapply( data$Sample , function(x){ paste( unlist(strsplit( x , "-" , fixed=TRUE))[1] ,  unlist(strsplit( x , "-" , fixed=TRUE))[2] ,  unlist(strsplit( x , "-" , fixed=TRUE))[3]  , sep = "-" ) })

data = unique(data)


data = data[ , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]

data$Effect[ data$Effect %in% "downstream_gene_variant" ] = "3'Flank" 
data$Effect[ data$Effect %in% "upstream_gene_variant" ] = "5'Flank" 
data$Effect[ data$Effect %in% c( "3_prime_UTR_variant" , "3_prime_UTR_variant&NMD_transcript_variant" ) ] = "3'UTR" 
data$Effect[ data$Effect %in% c( "5_prime_UTR_premature_start_codon_gain_variant" , "5_prime_UTR_variant" , "5_prime_UTR_variant&NMD_transcript_variant" ) ] = "5'UTR" 
data$Effect[ data$Effect %in% c( "synonymous_variant" , "initiator_codon_variant" , "synonymous_variant&NMD_transcript_variant" ) ] = "Silent" 
data$Effect[ data$Effect %in% c( "missense_variant&NMD_transcript_variant" , "missense_variant&splice_region_variant&NMD_transcript_variant" , "missense_variant" , "missense_variant&splice_region_variant" ) ] = "Missense_Mutation" 
data$Effect[ data$Effect %in% c( "intergenic_region" , "intron_variant" , "non_coding_transcript_exon_variant" , "mature_miRNA_variant" , "intergenic_variant" , "intron_variant&NMD_transcript_variant" , "intron_variant&non_coding_transcript_variant" , "regulatory_region_variant" ) ] = "Intron" 
data$Effect[ data$Effect %in% c( "splice_acceptor_variant" , "splice_acceptor_variant&NMD_transcript_variant" , "splice_acceptor_variant&non_coding_transcript_variant" , "splice_donor_variant" , "splice_donor_variant&coding_sequence_variant" , "splice_donor_variant&coding_sequence_variant&intron_variant" , "splice_donor_variant&NMD_transcript_variant" , "splice_donor_variant&non_coding_transcript_variant" , "splice_region_variant&3_prime_UTR_variant" , "splice_region_variant&3_prime_UTR_variant&NMD_transcript_variant" , "splice_region_variant&5_prime_UTR_variant" , "splice_region_variant&intron_variant&NMD_transcript_variant" , "splice_region_variant&intron_variant&non_coding_transcript_variant" , "splice_region_variant&stop_retained_variant" , "splice_region_variant&synonymous_variant&NMD_transcript_variant" , "splice_acceptor_variant&intron_variant" , "splice_donor_variant&intron_variant" , "splice_region_variant" , "splice_region_variant&intron_variant" , "splice_region_variant&non_coding_transcript_exon_variant" , "splice_region_variant&synonymous_variant" ) ] = "Splice_Site" 
data$Effect[ data$Effect %in% c( "missense_variant" , "missense_variant&splice_region_variant" ) ] = "Missense_Mutation" 

data$Effect[ data$Effect %in% c( "frameshift_variant&NMD_transcript_variant" , "frameshift_variant" , "frameshift_variant&splice_region_variant" , "frameshift_variant&start_lost" , "frameshift_variant&stop_gained" ) ] = "Frame_Shift_Del" 
data$Effect[ data$Effect %in% c( "inframe_deletion" , "inframe_deletion&splice_region_variant" ) ] = "In_Frame_Del" 
data$Effect[ data$Effect %in% c( "inframe_insertion" , "conservative_inframe_insertion" , "disruptive_inframe_insertion" , "disruptive_inframe_insertion&splice_region_variant" ) ] = "In_Frame_Ins" 

data$Effect[ data$Effect %in% c( "start_lost&NMD_transcript_variant" , "stop_gained&NMD_transcript_variant" , "stop_lost" , "stop_retained_variant" , "stop_gained" , "stop_gained&conservative_inframe_insertion" , "stop_gained&splice_region_variant" , "start_lost" , "start_retained_variant" ) ] = "Nonsense_Mutation" 

data$Sample <- str_replace_all(data$Sample, "-", "_")
write.table( data , file=file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

