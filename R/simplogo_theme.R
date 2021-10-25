#' define the simplogo_main_theme and simplogo_ic_theme functions for ggplot ##
#' @export

simplogo_main_theme <- function(){
  font <- "Arial"   #assign font family up front

  `%+replace%` <- ggplot2::`%+replace%`

  ggplot2::theme_bw() %+replace%    #replace elements we want to change

    ggplot2::theme(

      #grid elements
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = -3, r = 5, b = 0, l = 0),
        axis.ticks.length.y = ggplot2::unit(.65, "cm"),
        axis.ticks.y = ggplot2::element_blank(),

      #text elements
        axis.title.x = ggplot2::element_text(color="black"),
        axis.title.y.left = ggplot2::element_text(color="black"),
        axis.title.x.bottom = ggplot2::element_text(color="black"),
        axis.text.x = ggplot2::element_text(angle = -90, vjust = 0.14, color="black"),
        axis.text.x.top = ggplot2::element_text(angle = -90, hjust = 1.2, color="black", margin = ggplot2::margin(t = 0, r = 0, b = 5, l = 0)),
        axis.text.y.left = ggplot2::element_text(color="black", hjust=1),
        axis.text.x.bottom = ggplot2::element_text(color="black", hjust = 0, vjust=0.5, margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 0)),

      #legend elements
        legend.position = "bottom",
        legend.margin=ggplot2::margin(t = -0.05, b = 0.15, unit='cm'),
        legend.text = ggplot2::element_text(color="black"),
        legend.title = ggplot2::element_text(color="black", margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0))

    )
}

simplogo_ic_theme <- function(){
  font <- "Arial"   #assign font family up front

  `%+replace%` <- ggplot2::`%+replace%`

  ggplot2::theme_bw() %+replace%    #replace elements we want to change

    ggplot2::theme(

      #grid elements
      plot.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1.5),

      #text elements
      axis.title.x = ggplot2::element_text(color="black"),
      axis.title.y.left = ggplot2::element_text(color="black"),
      axis.text.x.top = ggplot2::element_text(angle = -90, hjust = 1.2, vjust=0.5, color="black", margin = ggplot2::margin(t = 0, r = 0, b = 5, l = 0)),
      axis.text.y.left = ggplot2::element_text(color="black", margin=ggplot2::margin(l=10)),

      #legend elements
      legend.position = "bottom",
      legend.margin=ggplot2::margin(t = 0, b = 0.05, unit='cm'),
      legend.text = ggplot2::element_text(color="black"),
      legend.title = ggplot2::element_text(color="black")

    )
}
