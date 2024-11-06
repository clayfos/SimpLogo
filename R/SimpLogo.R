#' Convert protein sequence alignments to pseudo-1D colored representations.
#'
#' This function takes one or more sequence alignment file and converts it to a SimpLogo format for plotting.
#'
#' @param seq.dir Directory containing one or more protein sequence alignments. If multiple alignments are present, they must be of equal lengths (i.e., a master alignment split into multiple groups).
#' @param pattern File extension for sequence alignment files (e.g., "*fa").
#' @param group.names Character vector containing names for sequence groups (1 per alignment file; make sure it's in the same order as your sequence files). If NULL, names are taken from filenames.
#' @param lineage.names Character vector containing lineage assignments for each sequence group. If NULL, assignments are taken from the sequence group.names themselves. If specifying, you must provide assignments for every included sequence group (1 per group).
#' @param res.colors Character vector containing color hexcodes for the physicochemical amino acid residue groups. If NULL, will be default (RECOMMENDED; see GitHub README). If specifying, you must provide colors for all 8 residue groupings. Order for vector is as follows: H/K/R, D/E, A/G, Q/N/S/T, M/V/L/I, F/Y/W, P, and C. Keep in mind that a blending of certain colors can often be difficult to differentiate from a solid assignment.
#' @return Returns a list containing a formatted SimpLogo data table ready for plotting and the color scheme used.
#' @export

SimpLogo <- function(seq.dir, pattern = "*fa", group.names = NULL, lineage.names = NULL, res.colors = NULL){

  ## change number of sig figs (default is 7)
  options(digits=8)

  ## check if path exists
  if (file.exists(seq.dir)) {
    ## list file names for the individual fasta files
    seq.files <- list.files(path=seq.dir, pattern=pattern, all.files=TRUE, full.names=TRUE, include.dirs = FALSE)
  } else {
    stop("Could not find the specified path!")
  }


  ## convert to AAString object
  seq.list <- list()
  for (k in 1:length(seq.files)) {
    seq.list[[k]] <- ape::as.character.AAbin(ape::read.FASTA(seq.files[[k]], type="AA"))
  }

  ## stop if you haven't specified enough group names when manually labeling seq groups
  if (!is.null(group.names) & length(seq.list) != length(group.names)){
    stop("If manually specifying sequence group names, you must provide an assignment for every distinct group!")
  }

  ## stop if you haven't specified a lineage assignment for every seq group
  if (!is.null(lineage.names) & length(seq.list) != length(lineage.names)){
    stop("If manually specifying sequence group names, you must provide an assignment for every distinct group!")
  }

  ## check that all seq.files are of same length (same number of columns/positions)
  match.length <- length(seq.list[[1]][[1]]) ## get length of first sequence

  for (i in 1:length(seq.list)){
    if (length(seq.list[[i]][[1]]) != match.length){
      stop("Alignment widths are not all equal! Check your FASTA files.")
    }
  }

  ## number of sequence groups (# of files)
  num.archs <- length(seq.list)

  ## add architectural class names if desired
  if (is.null(group.names)){
    ## take name from filenames
    group.names <- c()
    group.names <- list.files(path=seq.dir, pattern=pattern, all.files=TRUE, full.names=FALSE, include.dirs = FALSE)
    names(seq.list) <- group.names
  } else {
    ## manually specify names (doublecheck your order!)
    names(seq.list) <- group.names
  }

  ## creating a position count matrix (set prob to false) ##
  ## make combined list of residue frequencies (grouping amino acids based on physicochemical properties)
  motif.list <- vector(mode = "list", length=num.archs)
  names(motif.list) <- names(seq.list)
  combined.motif.list <- vector(mode = "list", length = num.archs)
  names(combined.motif.list) <- names(seq.list)
  residue.type.list <- vector(mode = "list", length=9) ## 9 distinct categories of side chains
  names(residue.type.list) <- c("H/K/R","D/E","A/G","Q/N/S/T","M/V/L/I","F/Y/W","P","C","-")

  ## assign colors
  if (is.null(res.colors)){ ## use defaults
  motif.color.list <- c("H/K/R" = "#4363d8", "D/E" = "#e6194b", "A/G" = "#2F4F4F", "Q/N/S/T" = "#ffe119",
                        "M/V/L/I" = "#3cb44b", "F/Y/W" = "#9a6324", "P" = "purple", "C" = "#46f0f0", "-" = NA
  )
  }
  if (!is.null(res.colors)){ ## use custom colors
    if (length(res.colors)==length(residue.type.list)) {
      motif.color.list <- c(res.colors[1], res.colors[2], res.colors[3], res.colors[4],
                            res.colors[5], res.colors[6], res.colors[7], res.colors[8], res.colors[9])
    } else {
      stop("If manually specifying residue group colors, you must provide a color for every distinct group!")
    }
  }

  ## convert to RGB, create array with with the frequencies
  ##  9 different matrices, one for each residue group
  rgb.values.list <- vector(mode = "list", length = num.archs)
  names(rgb.values.list) <- names(seq.list)

  align.width <- length(seq.list[[1]][[1]])
  num.motif.colors <- length(motif.color.list)
  rgb.values.with.freqs.array <- array(dim=c(align.width,4,num.motif.colors),
                                       dimnames = list(c(seq(1,align.width,by=1)),
                                                       c("R","G","B","freq"),
                                                       c(names(residue.type.list))))
  for (i in 1:num.archs) {
    rgb.values.list[[i]] <- rgb.values.with.freqs.array
  }

  ## convert RGB to XYZ
  xyz.values.list <- vector(mode = "list", length = num.archs)
  names(xyz.values.list) <- names(seq.list)

  xyz.values.with.freqs.array <- array(dim=c(align.width,4,num.motif.colors),
                                       dimnames = list(c(seq(1,align.width,by=1)),
                                                       c("R","G","B","freq"),
                                                       c(names(residue.type.list))))
  for (i in 1:num.archs) {
    xyz.values.list[[i]] <- xyz.values.with.freqs.array
  }

  ## list of res types to combine for colSums operation
  restypes.for.summing <- vector(mode = "list", length = num.motif.colors)
  names(restypes.for.summing) <- names(residue.type.list)
  restypes.for.summing[[1]] <- c("H","K","R")
  restypes.for.summing[[2]] <- c("D","E")
  restypes.for.summing[[3]] <- c("A","G")
  restypes.for.summing[[4]] <- c("Q","N","S","T")
  restypes.for.summing[[5]] <- c("M","V","L","I")
  restypes.for.summing[[6]] <- c("F","Y","W")
  restypes.for.summing[[7]] <- c("P")
  restypes.for.summing[[8]] <- c("C")
  restypes.for.summing[[9]] <- c("-")

  for (i in 1:num.archs) {

    ## first convert to a counts matrix
    tmp.consensus <- Biostrings::consensusMatrix(substr(do.call(rbind, seq.list[[i]]), start=1, stop=align.width), as.prob = FALSE)
    ## check that all 20 amino acids + gap are present
    ## if not, add in empty rows
    if (nrow(tmp.consensus) != 21){
      full.aa.list <- unlist(restypes.for.summing)
      names(full.aa.list) <- full.aa.list
      missing.aa <- full.aa.list[!(full.aa.list %in% row.names(tmp.consensus))]
      missing.rows <- matrix(0, nrow=length(missing.aa), ncol = ncol(tmp.consensus))
      row.names(missing.rows) <- missing.aa
      tmp.consensus <- rbind(tmp.consensus, missing.rows)
      ## set row order to be consistent with other lists (and consensusMatrix default order)
      correct.row.order <- c("-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
      tmp.consensus <- tmp.consensus[correct.row.order,]
    }
    ## then convert to frequency matrix
    tmp.motif <- motifStack::pcm2pfm(tmp.consensus)
    motif.name <- paste0(names(seq.list[i]))
    ## define new position frequency matrix
    motif.list[[i]] <- new("pfm", mat=tmp.motif, name=motif.name,
                           color=motifStack::colorset(alphabet="AA",colorScheme="chemistry"))

    ## add residue types
    combined.motif.list[[i]] <- residue.type.list

    ## add residue types by physicochemical properties
    for (f in 1:(num.motif.colors)) {
      if (length(restypes.for.summing[[f]]) > 1){
        combined.motif.list[[i]][[f]] <- colSums(motif.list[[i]]@mat[restypes.for.summing[[f]],])
      } else { ## can't sum individual characters (last 3 restypes, P, C and gaps)
        combined.motif.list[[i]][[f]] <- motif.list[[i]]@mat[restypes.for.summing[[f]],]
      }
    }

    tmp.rgb.aperm.array <- aperm(xyz.values.with.freqs.array)
    for (h in 1:nrow(tmp.rgb.aperm.array[,1:3,])){
      ## convert the color (common name) to its RGB values
      tmp.rgb.aperm.array[h,1:3,] <- t(col2rgb(motif.color.list[row.names(tmp.rgb.aperm.array)[h]], alpha=F))
      ## add in residue type freq values
      tmp.rgb.aperm.array[h,4,] <- unlist(combined.motif.list[[i]][row.names(tmp.rgb.aperm.array)[h]]) ## unlist the named list
    }

    ## revert
    tmp.rgb.unaperm.array <- aperm(tmp.rgb.aperm.array)

    rgb.values.list[[i]] <- tmp.rgb.unaperm.array

    tmp.xyz.aperm.array <- aperm(xyz.values.with.freqs.array)
    for (k in 1:nrow(tmp.xyz.aperm.array[,1:3,])){
      ## convert RGB to XYZ using rgb_to_xyz_array.R
      tmp.xyz.aperm.array[k,1:3,] <- apply(as.matrix(rgb.values.list[[i]][,1:3,row.names(tmp.xyz.aperm.array)[k]]), 1, rgb_to_xyz_array)
      ## add freq values
      tmp.xyz.aperm.array[k,4,] <- unlist(combined.motif.list[[i]][row.names(tmp.xyz.aperm.array)[k]])
    }

    ## revert
    tmp.xyz.unaperm.array <- aperm(tmp.xyz.aperm.array)

    xyz.values.list[[i]] <- tmp.xyz.unaperm.array

  }

  ## get color list for converting
  top.color.list.df <- data.frame(motif.color.list)
  top.color.list.df$top.type <- row.names(top.color.list.df)
  colnames(top.color.list.df) <- c("top.color","top.type")

  secondary.color.list.df <- data.frame(motif.color.list)
  secondary.color.list.df$secondary.type <- row.names(secondary.color.list.df)
  colnames(secondary.color.list.df) <- c("secondary.color","secondary.type")


  ## define final xyz coordinates for each seq group
  final.xyz.list <- vector(mode = "list", length = num.archs)
  names(final.xyz.list) <- names(seq.list)

  ## xyz color coordinates for each position in the alignment (all residue types combined)
  final.xyz.matrix <- matrix(nrow=align.width, ncol=8, byrow=T)
  colnames(final.xyz.matrix) <- c("X","Y","Z","gap.pattern","gap.freq","top.type",
                                  "secondary.type","info.content") ## include IC in bits

  ## convert to hex for ggplotting
  final.hex.list <- vector(mode = "list", length = num.archs)
  names(final.hex.list) <- names(seq.list)

  for (i in 1:num.archs) {
    ## pre-allocate final.xyz.matrix
    final.xyz.matrix <- matrix(0, nrow=align.width, ncol=8)

    ## pre-calculate information content for the sequence alignment using ggseqlogo functions
    ## convert seq.list[[i]] object into an appropriate list of strings of equal length
    string.list <- lapply(seq.list[[i]], function(x) paste(x, collapse = ""))
    string.unlist <- unlist(string.list)
    tmp.pfm <- ggseqlogo:::makePFM(string.unlist)
    tmp.ic <- attr(tmp.pfm, "bits")

    colnames(final.xyz.matrix) <- c("X","Y","Z","gap.pattern","gap.freq","top.type","secondary.type","info.content")

    for (j in 1:align.width) {

      ## get xyz values for each position j for each group i
      xyz.values.for.summing <- xyz.values.list[[i]][j,1:4,]

      ## calculate top two freq characters (excluding gaps in last column)
      freq.col <- xyz.values.for.summing[4, 1:(num.motif.colors - 1)]
      top.indices <- order(freq.col, decreasing = TRUE)[1:2]
      top.freq.type <- colnames(xyz.values.for.summing)[top.indices[1]]
      secondary.freq.type <- colnames(xyz.values.for.summing)[top.indices[2]]

      ## create matrix of frequencies (for multiplication) ## 3 rows for X, Y and Z, 9 cols for 9 residue types
      freq.matrix <- matrix(xyz.values.for.summing[4, ], nrow=3, ncol=9, byrow = TRUE)
      ## compensate for positions with high gap (rest of the residues will be black)
      compensated.freq.matrix <- freq.matrix + (freq.matrix[, 6] / 8)
      ## multiply X, Y and Z coordinate by frequency
      multiplied.matrix <- xyz.values.for.summing[1:3, ] * compensated.freq.matrix
      multiplied.matrix[,"-"] <- 0       ## replace gap values zeroes
      ## Sum X, Y, and Z values
      summed.matrix <- rowSums(multiplied.matrix)
      ## Store the results
      final.xyz.matrix[j, 1:3] <- summed.matrix
      final.xyz.matrix[j, 4] <- "b"
      final.xyz.matrix[j, 5] <- xyz.values.for.summing["freq", "-"]
      final.xyz.matrix[j, 6] <- top.freq.type
      final.xyz.matrix[j, 7] <- secondary.freq.type
      ## add in IC
      final.xyz.matrix[j, 8] <- tmp.ic[j]

    }
    final.xyz.list[[i]] <- final.xyz.matrix
    tmp.mat <- apply(final.xyz.list[[i]][,1:3], 2, as.numeric)
    final.hex.list[[i]] <- data.frame(schemr::xyz_to_hex(tmp.mat, transformation = "sRGB", linear_func = NULL))
    final.hex.list[[i]]$position <- seq(1,align.width, by=1)
    final.hex.list[[i]]$gap.freq <- as.numeric(final.xyz.list[[i]][,"gap.freq"])
    final.hex.list[[i]]$info.content <- as.numeric(final.xyz.list[[i]][,"info.content"])
    final.hex.list[[i]]$top.type <- final.xyz.list[[i]][,"top.type"]
    final.hex.list[[i]]$secondary.type <- final.xyz.list[[i]][,"secondary.type"]
    colnames(final.hex.list[[i]]) <- c("color", "position", "gap.freq", "info.content", "top.type","secondary.type")
    ## now add in "ideal" colors for the top.type and secondary.type
    final.hex.list[[i]] <- plyr::join(final.hex.list[[i]], top.color.list.df, by="top.type", type="left")
    final.hex.list[[i]] <- plyr::join(final.hex.list[[i]], secondary.color.list.df, by="secondary.type", type="left")
  }

  ## add arch type as a id and reshape
  for (i in 1:length(final.hex.list)){
    final.hex.list[[i]]$arch <- names(final.hex.list)[i]
  }

  if (is.null(lineage.names)){
    for (i in 1:length(final.hex.list)){
      final.hex.list[[i]]$lineage <- final.hex.list[[i]]$arch
    }
  } else {
    for (i in 1:length(final.hex.list)){
      final.hex.list[[i]]$lineage <- lineage.names[i]
      final.hex.list[[i]]$lineage <- factor(final.hex.list[[i]]$lineage)
    }

  }
  ## bind together with arch as id
  final.bound <- dplyr::bind_rows(final.hex.list, .id = "arch")
  final.bound$residue <- as.character(final.bound$position)
  ## add in dummy factor for facetting IC plot
  final.bound$dummy <- "IC (bits)"
  simplogo.list <- list("results.table" = final.bound, "color.scheme" = motif.color.list)
  return(simplogo.list)

}
