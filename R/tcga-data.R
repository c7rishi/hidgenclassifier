#' TCGA somatic mutation annotation data
#' from 10 common cancer sites
#'
#' A subset of the publicly available NIH's the Cancer Genome Atlas (TCGA)
#' whole exome sequencing data containing non-synonymous SNV mutations.
#' More specifically,
#' the subset of the TCGA data with
#' \code{Variant_Type == "SNP"} and  tumors from the following 10 cancer sites:
#' BREAST, COLORECTAL, ESOPHAGEAL, KIDNEY, LIVER,
#' LUNG, OVARIAN, PANCREATIC, PROSTATE, and SKIN
#' @format A dataset with 1,395,665 rows and 14 columns:
#'
#' \describe{
#' \item{patient_id}{the patient (tumor) label. Obtained by extracting
#' the first 16 characters of the column'Tumor_Sample_Barcode' in
#' the original TCGA data.}
#' \item{Hugo_Symbol}{the gene label}
#' \item{Chromosome}{Chromosome label}
#' \item{Variant}{the variant label. Obtained by concatenating
#' the columns labeled
#' 'Hugo_Symbol', "HGVSp_Short', 'Chromosome', 'Start_Position', 'Tumor_Seq_Allele1',
#' and 'Tumor_Seq_Allele2' in the original TCGA data}
#' \item{Start_Position}{Start Position of alteration on the chromosome}
#' \item{End_Position}{End Position of alteration on the chromosome}
#' \item{Tumor_Seq_Allele1}{.}
#' \item{Tumor_Seq_Allele2}{.}
#' \item{Reference_Allele}{.}
#' \item{HGVSp_Short}{Protein code of alteration}
#' \item{Variant_Type}{Type of alteration}
#' \item{Variant_Classification}{.}
#' \item{CANCER_HISTOLOGY}{Histological subtype for each tumor}
#' \item{CANCER_SITE}{Cancer site of origin of each tumor}
#' }
#'
#' @source
#' https://github.com/cBioPortal/datahub
#'

"tcga"
