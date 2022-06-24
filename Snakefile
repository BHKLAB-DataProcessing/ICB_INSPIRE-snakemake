from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]
data_source  = "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_INSPIRE-data/main/"
large_data_source = "https://github.com/BHKLAB-Pachyderm/ICB_INSPIRE-data/blob/main/"

rule get_MultiAssayExp:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    resources:
        mem_mb=4000
    shell:
        """
        Rscript -e \
        '
        load(paste0("{prefix}", "annotation/Gencode.v40.annotation.RData"))
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "INSPIRE", input_dir = paste0("{prefix}", "processed")), 
            "{prefix}{filename}"
        );
        '
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true -O {prefix}annotation/Gencode.v40.annotation.RData 
        """

rule format_snv:
    output:
        S3.remote(prefix + "processed/SNV.csv")
    input:
        S3.remote(prefix + "download/INSPIRE.combinedVariants.snv.released.112020.renamed.vcf.gz"),
        S3.remote(prefix + "download/INSPIRE.combinedVariants.indel.released.112020.renamed.vcf.gz")
    resources:
        mem_mb=5000
    shell:
        """
        Rscript scripts/Format_SNV.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_data:
    output:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    input:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/gene-expression-matrix-TPM-final.tsv"),
        S3.remote(prefix + "processed/SNV.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_Data.R \
        {prefix}download \
        {prefix}processed \
        """

rule download_data:
    output:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/gene-expression-matrix-TPM-final.tsv"),
        S3.remote(prefix + "download/INSPIRE.combinedVariants.snv.released.112020.renamed.vcf.gz"),
        S3.remote(prefix + "download/INSPIRE.combinedVariants.indel.released.112020.renamed.vcf.gz")
    resources:
        mem_mb=2000
    shell:
        """
        wget {data_source}CLIN.txt -O {prefix}download/CLIN.txt
        wget {large_data_source}gene-expression-matrix-TPM-final.tsv?raw=true -O {prefix}download/gene-expression-matrix-TPM-final.tsv
        wget {large_data_source}INSPIRE.combinedVariants.snv.released.112020.renamed.vcf.gz?raw=true -O {prefix}download/INSPIRE.combinedVariants.snv.released.112020.renamed.vcf.gz
        wget {large_data_source}INSPIRE.combinedVariants.indel.released.112020.renamed.vcf.gz?raw=true -O {prefix}download/INSPIRE.combinedVariants.indel.released.112020.renamed.vcf.gz
        """ 