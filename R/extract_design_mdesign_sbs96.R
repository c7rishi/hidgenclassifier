#' Extract Design-Metadesign (XU) matrix for SBS-96
#' meta-feature categories from a maf file
#'
#' @inheritParams extract_design
#'
#' @param maf mutation annotation file --
#' a data frame-like object with at least six columns for:
#' Chromosome, Start_position, End_position,
#' Reference_Allele, Tumor_Seq_Allele2 and Sample ID.
#' NOTE: uniqueness of rows of maf is assumed.
#' @param return_vranges_object logical. Should the intermediate
#' \code{'VRanges'} object from \code{{SomaticSignatures}} cataloging
#' SBS-96 contexts for all mutations in \code{maf} be returned as an
#' attribute named \code{VRanges} of the output? Defaults to FALSE.
#'
#'
#' @details
#' This function calls uses functions
#' VariantAnnotation::VRanges(), IRanges::IRanges(),
#' SomaticSignatures::ucsc() and SomaticSignatures::mutationContext(),
#' together with the BSgenome.Hsapiens.UCSC.hg19 database to obtain
#' the 96 single base substitution contexts for each mutation in the
#' \code{maf} file in the form of a 'VRanges' object from {SomaticSignatures}.
#' Then using SomaticSignatures::motifMatrix a 96 x n_tumors is calculated,
#' which is subsequently transposed and converted into a sparseMatrix object
#' in the form of a n_patient x 96 design matrix
#'
#' @note
#'
#' The bioconductor packages {SomaticSignatures} and {VariantAnnotation} are not
#' automatically installed with hidgenclassifier. Please install them separately.
#'
#' Using any {SomaticSignatures} function triggers the loading of
#' {proxy} and {GGally} packages, which overwrites defaults S3 methods
#' of a few functions from {ggplot2} and {registry}.
#'
#' @return
#' An n_tumor x 96 sparse dgCMatrix, with (i, j)th entry providing the total
#' number of variants in tumor i associated with j-th SBS-96 category.
#'
#' @examples
#' data("impact")
#' sbs96_mdesign <- extract_design_mdesign_sbs96(
#'   maf = impact,
#'   chromosome_col = "Chromosome",
#'   start_position_col = "Start_Position",
#'   end_position_col = "End_Position",
#'   ref_col = "Reference_Allele",
#'   alt_col = "Tumor_Seq_Allele2",
#'   sample_id_col = "patient_id"
#' )
#' dim(sbs96_mdesign)
#'
#' @author Ronglai Shen, Saptarshi Chakraborty
#' @export

extract_design_mdesign_sbs96 <- function(
  maf,
  chromosome_col = "Chromosome",
  start_position_col = "Start_Position",
  end_position_col = "End_Position",
  ref_col = "Reference_Allele",
  alt_col = "Tumor_Seq_Allele2",
  sample_id_col = "sample",
  return_vranges_obj = FALSE,
  ...
) {

  reqd_pkgs <- c("VariantAnnotation", "SomaticSignatures")

  for (pkg in reqd_pkgs) {
    stopifnot(requireNamespace(pkg))
  }

  vranges_obj <- VariantAnnotation::VRanges(
    seqnames = paste0("chr", maf[[chromosome_col]]),
    ranges = IRanges::IRanges(
      maf[[start_position_col]],
      maf[[end_position_col]]
    ),
    ref = maf[[ref_col]],
    alt = maf[[alt_col]],
    sampleNames = maf[[sample_id_col]]
  ) %>%
    SomaticSignatures::ucsc() %>%
    SomaticSignatures::mutationContext(
      BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
      unify = TRUE
    )

  #the sca_mm is a 96 by n mutation signature matrix you add to meta-regression
  sca_mm <- vranges_obj %>%
    SomaticSignatures::motifMatrix(
    group = "sampleNames",
    normalize = F
  ) %>%
    t() %>%
    as.matrix() %>%
    Matrix(sparse = TRUE)

  if (return_vranges_obj) {
    attr(sca_mm, "VRanges") <- vranges_obj
  }

  sca_mm
}

