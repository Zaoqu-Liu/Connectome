#' FANTOM5 Ligand-Receptor Database (Human)
#'
#' A dataset containing ligand-receptor pairs curated from the FANTOM5 project,
#' converted to human gene symbols.
#'
#' @format A data frame with 2557 rows and 17 columns:
#' \describe{
#'   \item{Pair.Name}{Unique identifier for the ligand-receptor pair}
#'   \item{Ligand.ApprovedSymbol}{HGNC approved symbol for the ligand}
#'   \item{Ligand.Name}{Full name of the ligand}
#'   \item{Receptor.ApprovedSymbol}{HGNC approved symbol for the receptor}
#'   \item{Receptor.Name}{Full name of the receptor}
#'   \item{DLRP}{DLRP database reference}
#'   \item{HPMR}{HPMR database reference}
#'   \item{IUPHAR}{IUPHAR database reference}
#'   \item{HPRD}{HPRD database reference}
#'   \item{STRING.binding}{STRING binding evidence}
#'   \item{STRING.experiment}{STRING experimental evidence}
#'   \item{HPMR.Ligand}{HPMR ligand identifier}
#'   \item{HPMR.Receptor}{HPMR receptor identifier}
#'   \item{PMID.Manual}{PubMed ID for manual curation}
#'   \item{Pair.Source}{Source classification: "known" or "novel"}
#'   \item{Pair.Evidence}{Evidence level: "literature supported" or "putative"}
#'   \item{mode}{Signaling mode/category of the interaction}
#' }
#' @source Ramilowski, J.A. et al. A draft network of ligand-receptor-mediated 
#'   multicellular signalling in human. Nat Commun 6, 7866 (2015).
#'   \doi{10.1038/ncomms8866}
#' @seealso \code{\link{ncomms8866_mouse}}, \code{\link{ncomms8866_rat}}, 
#'   \code{\link{ncomms8866_pig}}
"ncomms8866_human"

#' FANTOM5 Ligand-Receptor Database (Mouse)
#'
#' A dataset containing ligand-receptor pairs from the FANTOM5 project,
#' converted to mouse gene symbols using ortholog mapping.
#'
#' @format A data frame with ligand-receptor pairs. See \code{\link{ncomms8866_human}}
#'   for column descriptions.
#' @source Ramilowski, J.A. et al. Nat Commun 6, 7866 (2015). \doi{10.1038/ncomms8866}
#' @seealso \code{\link{ncomms8866_human}}
"ncomms8866_mouse"

#' FANTOM5 Ligand-Receptor Database (Rat)
#'
#' A dataset containing ligand-receptor pairs from the FANTOM5 project,
#' converted to rat gene symbols using ortholog mapping.
#'
#' @format A data frame with ligand-receptor pairs. See \code{\link{ncomms8866_human}}
#'   for column descriptions.
#' @source Ramilowski, J.A. et al. Nat Commun 6, 7866 (2015). \doi{10.1038/ncomms8866}
#' @seealso \code{\link{ncomms8866_human}}
"ncomms8866_rat"

#' FANTOM5 Ligand-Receptor Database (Pig)
#'
#' A dataset containing ligand-receptor pairs from the FANTOM5 project,
#' converted to pig gene symbols using ortholog mapping.
#'
#' @format A data frame with ligand-receptor pairs. See \code{\link{ncomms8866_human}}
#'   for column descriptions.
#' @source Ramilowski, J.A. et al. Nat Commun 6, 7866 (2015). \doi{10.1038/ncomms8866}
#' @seealso \code{\link{ncomms8866_human}}
"ncomms8866_pig"

#' FANTOM5 Ligand-Receptor Database (Original)
#'
#' The original FANTOM5 ligand-receptor database before species conversion.
#'
#' @format A data frame with the original ligand-receptor pairs from FANTOM5.
#' @source Ramilowski, J.A. et al. Nat Commun 6, 7866 (2015). \doi{10.1038/ncomms8866}
"ncomms8866"

#' Gene Symbol Translator
#'
#' A lookup table for translating gene symbols between species.
#'
#' @format A data frame containing ortholog mappings between species.
#' @keywords internal
"translator"
