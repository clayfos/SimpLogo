#' Plot the SimpLogo representation using ggplot2.
#'
#' This function is the plotting engine for SimpLogo. It takes the output of the SimpLogo() function and creates the graphical representation using ggplot2.
#'
#' @param simplogo Formatted output of SimpLogo() function.
#' @param plot.ic Boolean whether to display position-wise information content.
#' @param position.start Number of first position residue (if using reference protein that doesn't start at position 1).
#' @return Returns a ggplot2 object containing the SimpLogo representation.
#' @export

SimpLogoPlot <- function(simplogo, plot.ic = TRUE, position.start = NULL){
  ## primary plotting window
  primary.logo <- ggplot2::ggplot(data = data.frame(simplogo), ggplot2::aes(x=as.numeric(residue), y=arch)) +
    ggplot2::geom_tile(## bottom color strip
      position = ggplot2::position_nudge(y = ifelse(simplogo$gap.freq < 0.5, -0.31*(1-simplogo$gap.freq), -0.31*(1-simplogo$gap.freq)-0.05)),
      #position = ggplot2::position_nudge(y = (0.25*(1-simplogo$gap.freq))),
      height=0.1, width=0.9, size=0.1,  color="black",
      fill=simplogo$secondary.color) +
    ggplot2::geom_tile(## top color strip
      position = ggplot2::position_nudge(y = ifelse(simplogo$gap.freq < 0.5, 0.31*(1-simplogo$gap.freq), 0.31*(1-simplogo$gap.freq)+0.05)),
      height=0.1, width=0.9, size=0.1, color="black",
      fill=simplogo$top.color) +
    ggplot2::geom_tile(## add main rectangle (set heigh max to 0.5 x % non-gap characters)
      height=(0.5*(1-simplogo$gap.freq)), width=0.9, size=0.15,
      fill=simplogo$color, color="black") +
    ## artificially construct a legend with idealized (100%) residue type colors
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#4363d8"), size=0, color="black", pch=21) + ## H/K/R
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#e6194b"), size=0, color="black", pch=21) + ## D/E
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#2F4F4F"), size=0, color="black", pch=21) + ## A/G
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#ffe119"), size=0, color="black", pch=21) + ## Q/N/S/T
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#3cb44b"), size=0, color="black", pch=21) + ## M/V/L/I
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#9a6324"), size=0, color="black", pch=21) + ## F/Y/W
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#ffffff"), size=0, color="black", pch=21) + ## P
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill="#46f0f0"), size=0, color="black", pch=21) + ## C
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=NA), size=0, color="black", pch=21) + ## - or GAP
    simplogo_main_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=5))) + ## increase the size of the legend points, without increasing size of the fake geom_points
    ggplot2::labs(x="Residue position", y="") +
    #ggplot2::scale_y_discrete(position = "left") +
    ggplot2::scale_y_discrete(position = "left", limits=rev(levels(simplogo$arch))) +
    #ggplot2::coord_cartesian(clip="off", xlim=c(1,max(simplogo$position))) +
    ggplot2::scale_fill_identity(guide = "legend",
                                 name = paste0("Residue Type (100%)"),
                                 labels = c("H/K/R","D/E","A/G","Q/N/S/T","M/V/L/I","F/Y/W","P","C","-"),
                                 breaks = c("H/K/R" = "#4363d8", "D/E" = "#e6194b", "A/G" = "#2F4F4F", "Q/N/S/T" = "#ffe119",
                                            "M/V/L/I" = "#3cb44b", "F/Y/W" = "#9a6324", "P" = "#ffffff", "C" = "#46f0f0", "-" = NA)
    ) + ## use this to force it to create a legend with the colors listed
    if (!is.null(position.start)) {
      ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo$position), by=5), labels=seq(position.start,position.start+max(simplogo$position)-1, by=5), expand = c(0, 0)
      )
    } else {
      ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo$position), by=5), labels=seq(1,max(simplogo$position), by=5), expand = c(0, 0)
      )
    }


  if (plot.ic == TRUE){
    logo.ic <- ggplot2::ggplot(data = data.frame(simplogo), ggplot2::aes(x=as.numeric(residue), y=info.content)) +
      ggplot2::geom_line(size=0.8, ggplot2::aes(color=lineage, group=interaction(arch,lineage))) + ## use interaction so we get all replicates (no distinguishing) colored by lineage
      simplogo_ic_theme() +
      ggplot2::labs(x="Residue Position", y="IC (bits)", color = "Protein lineage") +
      #ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      if (!is.null(position.start)) {
        ## add axis on top
        ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo$position), by=5), position = "top", labels=seq(position.start,position.start+max(simplogo$position)-1, by=5), expand = c(0, 0)
        )
      } else {
        ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo$position), by=5), position = "top", labels=seq(1,max(simplogo$position), by=5), expand = c(0, 0)
        )
      }
  } else {
    logo.ic <- NULL
    if (!is.null(position.start)) {
      ## overwrite original axis
      primary.logo <- primary.logo + ggplot2::scale_x_continuous(limits=NULL, sec.axis = dup_axis(), breaks=seq(1,max(simplogo$position), by=5), labels=seq(position.start,position.start+max(simplogo$position)-1, by=5), expand = c(0, 0)
      )
    } else {
      ## overwrite original axis
      primary.logo <- primary.logo + ggplot2::scale_x_continuous(limits=NULL, sec.axis = dup_axis(), breaks=seq(1,max(simplogo$position), by=5), labels=seq(1,max(simplogo$position), by=5), expand = c(0, 0)
      )
    }
  }



  final.plot <- cowplot::plot_grid(logo.ic, primary.logo,
                                   ncol = 1, align = 'v',
                                   axis = "lr", rel_heights = c(0.5, 1)
  )

  plot.list <- list("final.plot" = final.plot, "primary.plot" = primary.logo, "ic.plot" = logo.ic)
  return(plot.list)
}
