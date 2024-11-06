#' Plot the SimpLogo representation using ggplot2.
#'
#' This function is the plotting engine for SimpLogo. It takes the output of the SimpLogo() function and creates the graphical representation using ggplot2.
#'
#' @param simplogo Formatted list output of SimpLogo() function.
#' @param plot.ic Boolean whether to display position-wise information content.
#' @param position.start Number of first position residue (if using reference protein that doesn't start at position 1).
#' @return Returns a ggplot2 object containing the SimpLogo representation.
#' @export

SimpLogoPlot <- function(simplogo, plot.ic = TRUE, position.start = NULL){
  ## primary plotting window
  simplogo.datatable <- simplogo[[1]]
  color.scheme <- simplogo[[2]]
  primary.logo <- ggplot2::ggplot(data = data.frame(simplogo.datatable), ggplot2::aes(x=as.numeric(residue))) +

  facet_wrap(~arch, ncol=1, strip.position="left") +

  ## make top color strip (same coordinates as main rect, just shift up a bit)
  ggplot2::geom_rect(## top color strip
    ymax=0+(0.5*(1-simplogo.datatable$gap.freq))/2+0.05,
    ymin=0, ## start at 0
    xmax=as.numeric(simplogo.datatable$residue)+0.47,
    xmin=as.numeric(simplogo.datatable$residue)-0.47,
    size=0.2,
    fill=simplogo.datatable$top.color,
    color="black") +

  ## make bottom color strip (same coordinates as main rect, just shift down a bit)
  ggplot2::geom_rect(## top color strip
    ymax=0, ## start at 0
    ymin=0-(0.5*(1-simplogo.datatable$gap.freq))/2-0.05,
    xmax=as.numeric(simplogo.datatable$residue)+0.47,
    xmin=as.numeric(simplogo.datatable$residue)-0.47,
    size=0.2,
    fill=simplogo.datatable$secondary.color,
    color="black") +

  ## make middle main rectangle (set height max to 0.5 * % non-gap characters)
  ggplot2::geom_rect(## top color strip
    ymax=0+(0.5*(1-simplogo.datatable$gap.freq))/2,
    ymin=0-(0.5*(1-simplogo.datatable$gap.freq))/2,
    xmax=as.numeric(simplogo.datatable$residue)+0.47,
    xmin=as.numeric(simplogo.datatable$residue)-0.47,
    size=0.2,
    fill=simplogo.datatable$color,
    color="black") +

  # ylim(-0.3,0.3) +

    ## artificially construct a legend with idealized (100%) residue type colors
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[1]), size=0, color="black", pch=21) + ## H/K/R
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[2]), size=0, color="black", pch=21) + ## D/E
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[3]), size=0, color="black", pch=21) + ## A/G
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[4]), size=0, color="black", pch=21) + ## Q/N/S/T
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[5]), size=0, color="black", pch=21) + ## M/V/L/I
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[6]), size=0, color="black", pch=21) + ## F/Y/W
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[7]), size=0, color="black", pch=21) + ## P
    ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=color.scheme[8]), size=0, color="black", pch=21) + ## C
    # ggplot2::geom_point(x=0, y=0, ggplot2::aes(fill=NA), size=0, color="black", pch=21) + ## - or GAP
    simplogo_main_theme() +
    theme(strip.placement = "outside",
          strip.switch.pad.grid = unit(0.2, "in"),
          legend.background = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=4))) + ## increase the size of the legend points, without increasing size of the fake geom_points
    ggplot2::labs(x="Residue position", y="") +
    ggplot2::scale_y_continuous(limits=c(-0.3,0.3), breaks=NULL, labels=NULL) +
    #ggplot2::scale_y_discrete(position = "left") +
    # ggplot2::scale_y_discrete(position = "left", limits=rev(levels(simplogo.datatable$arch))) +
    #ggplot2::coord_cartesian(clip="off", xlim=c(1,max(simplogo.datatable$position))) +
    ggplot2::scale_fill_identity(guide = "legend",
                                 name = paste0("Residue Type (100%) "),
                                 labels = c("H/K/R","D/E","A/G","Q/N/S/T","M/V/L/I","F/Y/W","P","C"),
                                 breaks = c("H/K/R" = color.scheme[1], "D/E" = color.scheme[2], "A/G" = color.scheme[3], "Q/N/S/T" = color.scheme[4],
                                            "M/V/L/I" = color.scheme[5], "F/Y/W" = color.scheme[6], "P" = color.scheme[7], "C" = color.scheme[8])
    ) + ## use this to force it to create a legend with the colors listed
    if (!is.null(position.start)) {
      ggplot2::scale_x_continuous(limits=c(0,max(simplogo.datatable$position)+1), breaks=seq(1,max(simplogo.datatable$position), by=5), labels=seq(position.start,position.start+max(simplogo.datatable$position)-1, by=5), expand = c(0, 0)
      )
    } else {
      ggplot2::scale_x_continuous(limits=c(0,max(simplogo.datatable$position)+1), breaks=seq(1,max(simplogo.datatable$position), by=5), labels=seq(1,max(simplogo.datatable$position), by=5), expand = c(0, 0)
      )
    }


  if (plot.ic == TRUE){
    logo.ic <- ggplot2::ggplot(data = data.frame(simplogo.datatable), ggplot2::aes(x=as.numeric(residue), y=info.content)) +
      ggplot2::geom_line(size=0.8, ggplot2::aes(color=lineage, group=interaction(arch,lineage))) + ## use interaction so we get all replicates (no distinguishing) colored by lineage
      simplogo_ic_theme() +
      ggplot2::labs(x="Residue Position", y="IC (bits)", color = "Protein lineage") +
      #ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      if (!is.null(position.start)) {
        ## add axis on top
        ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo.datatable$position), by=5), position = "top", labels=seq(position.start,position.start+max(simplogo.datatable$position)-1, by=5), expand = c(0, 0)
        )
      } else {
        ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo.datatable$position), by=5), position = "top", labels=seq(1,max(simplogo.datatable$position), by=5), expand = c(0, 0)
        )
      }
  } else {
    logo.ic <- NULL
    if (!is.null(position.start)) {
      ## overwrite original axis
      primary.logo <- primary.logo + ggplot2::scale_x_continuous(limits=NULL, sec.axis = dup_axis(), breaks=seq(1,max(simplogo.datatable$position), by=5), labels=seq(position.start,position.start+max(simplogo.datatable$position)-1, by=5), expand = c(0, 0)
      )
    } else {
      ## overwrite original axis
      primary.logo <- primary.logo + ggplot2::scale_x_continuous(limits=NULL, sec.axis = dup_axis(), breaks=seq(1,max(simplogo.datatable$position), by=5), labels=seq(1,max(simplogo.datatable$position), by=5), expand = c(0, 0)
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
