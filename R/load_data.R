read_litter <- function(file, sheet, range) {
  res <- readxl::read_excel(file, sheet = sheet, range = range,
                     col_types = c(rep("text", 4), "date", "numeric", "text", "text", "text"),
                     col_names = c("litter_id", "id",	"sex",	"genotype",	"date_birth", "age_weeks",	"date_dissection",	"dissection",	"fate"))

  cols_to_fill <- c("litter_id","date_dissection", "dissection",	"fate")
  for(c in cols_to_fill) {
    res[[c]][is.na(res[[c]])] <- res[[c]][1]
  }

  res
}

read_all_litters <- function(file, sheet, ranges) {
  res <- list()
  for(range in ranges) {
    res[[range]] <- read_litter(file, sheet, range)
  }

  do.call(rbind, res)
}

copy_merged_columns_name_repair <- function(names) {
  last_name_loc <- 1
  new_names <- rep(NA_character_, length(names))
  new_names[1] <- names[1]
  for(i in 2:length(names)) {
    if(names[i] != "") {
      if(i < length(names) && names[i + 1] == "") {
        new_names[i] <- paste0(names[i], "___rep", 1)
      } else {
        new_names[i] <- names[i]
      }
      last_name_loc <- i
    } else {
      rep_id <- i - last_name_loc + 1
      new_names[i] <- paste0(names[last_name_loc], "___rep", rep_id)
    }
  }
  vctrs::vec_as_names(new_names, repair = "unique")
}

primer_to_short <- function(x) {
  gsub(pattern = "^(MR[0-9]+/[0-9]+).*$", replacement = "\\1", x = x)
}

#' @export
#' @importFrom tidyselect starts_with
read_qpcr <- function(file, sheet, range) {
  qpcr <- readxl::read_excel(file, sheet = sheet, range = range, .name_repair = copy_merged_columns_name_repair)
  qpcr <- dplyr::filter(qpcr, !is.na(`animal no.`))
  last_run <- NULL
  for(i in 1:nrow(qpcr)) {
    if(!is.na(qpcr[["qPCR run"]][i])) {
      last_run <- qpcr[["qPCR run"]][i]
    } else {
      qpcr[["qPCR run"]][i] <- last_run
    }
  }
  qpcr_long <- tidyr::pivot_longer(qpcr,
                      c(starts_with("MR"), starts_with("Ubb"), starts_with("Gapdh"), starts_with("Actb")),
                      names_to = c("primer", "replicate"), names_sep = "___", values_to = "ct")

  qpcr_long <- dplyr::filter(qpcr_long, !is.na(ct), ct > 0)
  qpcr_long <- dplyr::mutate(qpcr_long, primer_short = primer_to_short(primer))

  qpcr_long
}
