library(argparse)
library(VariantAnnotation)

# Create the command-line argument parser
parser <- ArgumentParser(description = "Process VCF file and extract non-zero entries.")

# Add the input file argument
parser$add_argument("-f", "--file", dest = "file", help = "VCF file path")

# Parse the command-line arguments
args <- parser$parse_args()

# Read VCF file
fl <- args$file#system.file(, package="VariantAnnotation")
# get file directory
dir <- dirname(fl)
file_pre <- sub("\\..*$", "", tools::file_path_sans_ext(basename(fl)))
replace_sub_str <- function(start_ind,end_ind, replace_sub_seq, original_seq) {

  new_seq <- paste0(substr(original_seq, 1, start_ind - 1),
                    replace_sub_seq,
                    substr(original_seq, end_ind + 1, nchar(original_seq)))

  return(new_seq)
}

get_annotations <- function(ann) {
  mut_names <- names(ann)
  # create an empty dataframe with name, type, mut, and transcript
  df <- data.frame(SNV = character(),
                   mutation_type = character(),
                   mut = character(),
                   transcript = character(),
                   gene = character(),
                   mutation_number = character(),
                   stringsAsFactors = FALSE)

  for (a in seq_along(ann)) {
    ann_a <- ann[[a]]
    name <- mut_names[a]
    name_split <- strsplit(name, ";")[[1]]

    # retrieving the wild type
    start_indices <- sapply(name_split, function(x) {
      parts <- unlist(strsplit(x, "_"))
      return(as.numeric(parts[2]))
    })
    end_indices <- sapply(name_split, function(x) {
      parts <- unlist(strsplit(x, "_"))
      return(as.numeric(parts[2])+nchar(parts[3]))
    })
    wildtypes <- sapply(name_split, function(x) {
      parts <- unlist(strsplit(x, "_"))
      return(parts[3])
    })
    mutations <- sapply(name_split, function(x) {
      parts <- unlist(strsplit(x, "_"))
      return(parts[4])
    })

    min_start <- min(start_indices)
    max_end <- max(end_indices)

    wt_seq <- strrep("X", max_end-min_start)
    for (mut_ind in order(start_indices)){
      s <- start_indices[mut_ind] - min_start + 1
      e <- end_indices[mut_ind] - min_start
      seq <- wildtypes[mut_ind]
      wt_seq <- replace_sub_str(s,e, seq, wt_seq)

    }
    mutated_seqs <- c()
    for (mut_ind in order(start_indices)){
      s <- start_indices[mut_ind] - min_start + 1
      e <- end_indices[mut_ind] - min_start
      seq <- mutations[mut_ind]
      mutated_seqs <- c(mutated_seqs, replace_sub_str(s,e, seq, wt_seq))

    }
    for (var in ann_a) {
      split_str <- strsplit(var, "\\|")[[1]]
      mut_bp <- split_str[1]
      mutation_type <- split_str[2]
      mutation_number = which(mut_bp == mutated_seqs)
      transcript <- split_str[7]
      gene <- split_str[5]
      mut <- split_str[11]

      # add a row to the dataframe
      df <- rbind(df, data.frame(SNV = name,
                                 mutation_type = mutation_type,
                                 mut = mut,
                                 transcript = transcript,
                                 gene = gene,
                                 mutation_number = mutation_number,
                                 stringsAsFactors = FALSE))

  }
  }
  return(df)}


vcffile <- open(VcfFile(fl,index=paste(fl, "tbi", sep="."), yieldSize = 1000))
vcffile2 <- open(VcfFile(fl,index=paste(fl, "tbi", sep="."), yieldSize = 1000))

first_iteration <- TRUE
GT <- NULL

repeat {
  GT_i <- readGeno(vcffile, "GT")
  ANN_i <- readInfo(vcffile2, "ANN")


  if (first_iteration) {
    GT <- GT_i
    ANN <- ANN_i
    first_iteration <- FALSE
  } else {
    GT <- rbind(GT, GT_i)
    ANN <- append(ANN, ANN_i)
      }

  if (nrow(GT_i) == 0)
    break

  ## Do work on the chunk

}

non_zero_indices <- which(GT != "0/0", arr.ind = TRUE)
snvs <- row.names(GT)[non_zero_indices[, 1]]
individuals <- colnames(GT)[non_zero_indices[, 2]]
mut <- GT[non_zero_indices]
mutation <- data.frame(Mut = mut, Individual = individuals, SNV = snvs)
split_string <- function(x) {
  unlist(strsplit(x, "\\|"))
}

ann_2 <- get_annotations(ANN)


merged_df <- merge(mutation, ann_2, by = "SNV")

mut_file <- file.path(dir, paste0("mutation",".csv"))

write.csv(merged_df, mut_file, row.names = FALSE)


