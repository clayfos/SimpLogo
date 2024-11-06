#' Plot the SimpLogo representation using ggplot2.
#'
#' This function is the plotting engine for SimpLogo. It takes the output of the SimpLogo() function and creates the graphical representation using ggplot2.
#'
#' @param simplogo Formatted list output of SimpLogo() function.
#' @param plot.ic Boolean whether to display position-wise information content.
#' @param position.start Number of first position residue (if using reference protein that doesn't start at position 1).
#' @return Returns a ggplot2 object containing the SimpLogo representation.
#' @export

SimpLogoPlot <- function(simplogo, plot.ic = TRUE, position.start = NULL, group.labels = NULL){
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
    aes(fill=top.color)) +

  ## make top color strip bottom border (black)
  ggplot2::geom_rect(## top color strip
    ymax=0+(0.5*(1-simplogo.datatable$gap.freq))/2+0.0075,
    ymin=0, ## start at 0
    xmax=as.numeric(simplogo.datatable$residue)+0.47,
    xmin=as.numeric(simplogo.datatable$residue)-0.47,
    size=0.2,
    fill="white") +

  ## make bottom color strip (same coordinates as main rect, just shift down a bit)
  ggplot2::geom_rect(## top color strip
    ymax=0, ## start at 0
    ymin=0-(0.5*(1-simplogo.datatable$gap.freq))/2-0.05,
    xmax=as.numeric(simplogo.datatable$residue)+0.47,
    xmin=as.numeric(simplogo.datatable$residue)-0.47,
    size=0.2,
    aes(fill=secondary.color)) +

  ## make bottom color strip border (black)
  ggplot2::geom_rect(## top color strip
    ymax=0, ## start at 0
    ymin=0-(0.5*(1-simplogo.datatable$gap.freq))/2-0.0075,
    xmax=as.numeric(simplogo.datatable$residue)+0.47,
    xmin=as.numeric(simplogo.datatable$residue)-0.47,
    size=0.2,
    fill="white") +

  ## make middle main rectangle (set height max to 0.5 * % non-gap characters)
  ggplot2::geom_rect(## top color strip
    ymax=0+(0.5*(1-simplogo.datatable$gap.freq))/2,
    ymin=0-(0.5*(1-simplogo.datatable$gap.freq))/2,
    xmax=as.numeric(simplogo.datatable$residue)+0.47,
    xmin=as.numeric(simplogo.datatable$residue)-0.47,
    size=0.2,
    aes(fill=color)) +

  # ylim(-0.3,0.3) +

    simplogo_main_theme() +
    theme(strip.placement = "outside",
          strip.switch.pad.grid = unit(0.2, "in"),
          legend.background = element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=9),
          axis.title.x = element_text(size=9),
          axis.ticks.y=element_blank(),
          strip.background =element_rect(fill="#595959", color="#000000", size=0.2),
          strip.text = element_text(color="white", size=9),
          legend.text = element_text(size=9),
          legend.title=element_text(size=9),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          legend.spacing.x = unit(0.01, 'cm'),
          legend.key.size = unit(0.35, 'cm')) +
    ggplot2::labs(x="Residue position", y="") +
    ggplot2::scale_y_continuous(limits=c(-0.3,0.3), breaks=NULL, labels=NULL) +
    ggplot2::scale_fill_identity(guide = "legend",
                                 name = paste0("Residue Type (100%)   "),
                                 labels = c("H/K/R","D/E","A/G","Q/N/S/T","M/V/L/I","F/Y/W","P","C"),
                                 breaks = c("H/K/R" = color.scheme[1],
                                            "D/E" = color.scheme[2],
                                            "A/G" = color.scheme[3],
                                            "Q/N/S/T" = color.scheme[4],
                                            "M/V/L/I" = color.scheme[5],
                                            "F/Y/W" = color.scheme[6],
                                            "P" = color.scheme[7],
                                            "C" = color.scheme[8])) + ## use this to force it to create a legend with the colors listed
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=1),
                                                 byrow=TRUE, nrow=1)) + ## increase the size of the legend points, without increasing size of the fake geom_points
    if (!is.null(position.start)) {
      ggplot2::scale_x_continuous(limits=c(0,max(simplogo.datatable$position)+1), breaks=seq(1,max(simplogo.datatable$position), by=5), labels=seq(position.start,position.start+max(simplogo.datatable$position)-1, by=5), expand = c(0, 0), sec.axis=sec_axis(~., breaks=seq(1,max(simplogo.datatable$position), by=5)))
    } else {
      ggplot2::scale_x_continuous(limits=c(0,max(simplogo.datatable$position)+1), breaks=seq(1,max(simplogo.datatable$position), by=5), labels=seq(1,max(simplogo.datatable$position), by=5), expand = c(0, 0), sec.axis=sec_axis(~., breaks=seq(1,max(simplogo.datatable$position), by=5)))
    }


  if (plot.ic == TRUE){
    logo.ic <- ggplot2::ggplot(data = data.frame(simplogo.datatable), ggplot2::aes(x=as.numeric(residue), y=info.content)) +
      ggplot2::geom_line(size=0.6, ggplot2::aes(color=lineage, group=interaction(arch,lineage))) + ## use interaction so we get all replicates (no distinguishing) colored by lineage
      simplogo_ic_theme() +
      facet_wrap(~dummy, ncol=1, strip.position="left") +
      ggplot2::labs(x="Residue Position", y="", color = "Protein lineage") +
      #ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::scale_y_continuous(position = "right", sec.axis = dup_axis(labels = NULL)) +
      theme(axis.title.y.right = element_blank(),
            axis.ticks.y.left = element_blank(),
            axis.title.y = element_text(size=9, color="white"),
            axis.text.y.right = element_text(size=8),
            # strip.switch.pad.grid = unit(0.2, "in"),
            strip.background =element_rect(fill="#595959", color="#000000", size=0.2),
            strip.text = element_text(color="white", size=9),
            axis.title.x = element_blank(),
            legend.position = "top") +
      if (!is.null(position.start)) {
        ## add axis on top
        ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo.datatable$position), by=5), position = "top", labels=seq(position.start,position.start+max(simplogo.datatable$position)-1, by=5), expand = c(0, 0))
      } else {
        ggplot2::scale_x_continuous(limits=NULL, breaks=seq(1,max(simplogo.datatable$position), by=5), position = "top", labels=seq(1,max(simplogo.datatable$position), by=5), expand = c(0, 0))
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

  ## if too many sequence groups, the facet strip labels get cut off
  ## if more than 4 distinct groups, shrink the facet label text slightly
  if (length(unique(simplogo.datatable$arch)) > 4){
    primary.logo <- primary.logo + theme(strip.text.y.left = element_text(color="white", size=7.5))
    logo.ic <- logo.ic + theme(strip.text.y.left = element_text(color="white", size=7.5)) ## shrink IC plot label font to match
  }

  # final.plot.title <- ggdraw() +
  #                     draw_label(
  #                       "PsR for each group of life",
  #                       fontface = 'plain',
  #                       size=9,
  #                       x = 0,
  #                       hjust = 0 ) +
  #                     theme(
  #                       # add margin on the left of the drawing canvas,
  #                       # so title is aligned with left edge of first plot
  #                       plot.margin = margin(l=-4))

  final.plot <- cowplot::plot_grid(logo.ic,
                                   NULL,
                                   primary.logo,
                                   ncol = 1, align = 'hv',
                                   axis = "lr", rel_heights = c(0.35, 0.0225, 1.4),
                                   labels = c("A","","B")) ## don't forget blanks for NULL elements

  plot.list <- list("final.plot" = final.plot, "primary.plot" = primary.logo, "ic.plot" = logo.ic)
  return(plot.list)
}
