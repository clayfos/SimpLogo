SimpLogo
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![DOI](https://zenodo.org/badge/418649749.svg)](https://zenodo.org/badge/latestdoi/418649749)

<!-- badges: end -->

## Overview

The goal of SimpLogo is to develop a simplified, easy-to-read
visualization for a large group of multiple amino acid sequence
alignments that allows for convenient, simultaneous intergroup
comparisons with minimized graphical clutter.

## Installation

You can install SimpLogo from Github using the devtools package:

``` r
library(devtools)
install_github("clayfos/SimpLogo",build_vignettes = TRUE)
```

## Example

SimpLogo utilizes the ggplot2 system to convert sequence alignments (in
FASTA format) into pseudo-1D, color-coded rows. Because of this, all of
the customization options offered by ggplot2 can be used to modify the
resulting SimpLogo plot. Here is a basic example of analyzing a
directory of FASTA alignments together.

### Input data

``` r
library(SimpLogo)

## use available files (.fa) in the sample_alignments directory
## we're only using 5 alignments for this example, but the 23 different CheW-containing architectures referenced in the publication can be found in the chew_alignments directory
## all alignments must be of equal length (for the example, all sequences were aligned simultaneously to the CheW-like Pfam HMM, then split into separate sequence groups)

## check filenames
seq.files <- list.files(path="sample_alignments/", 
                        all.files=TRUE, 
                        full.names=TRUE, 
                        include.dirs = FALSE,
                        pattern = "*fa")
seq.files
#> [1] "sample_alignments/CheA.I.separated.hmm.aligned.chew.as.ref.fa"   
#> [2] "sample_alignments/CheV.I.separated.hmm.aligned.chew.as.ref.fa"   
#> [3] "sample_alignments/CheW.IB.separated.hmm.aligned.chew.as.ref.fa"  
#> [4] "sample_alignments/CheW.IC.separated.hmm.aligned.chew.as.ref.fa"  
#> [5] "sample_alignments/CheW.II.1.separated.hmm.aligned.chew.as.ref.fa"
#> [6] "sample_alignments/CheW.II.2.separated.hmm.aligned.chew.as.ref.fa"

## you can manually specify alignment group names (for plotting)
## default takes them from the filenames
## make sure you order them correctly (check the filenames, like above)
groups <- c("CheA.I", "CheV.I", "CheW.IB", "CheW.IC", "CheW.II.1", "CheW.II.2")

## manually specify lineages/seq groups (for grouping during plotting)
## must have 1 assignment for every seq alignment/group being analyzed
lineage.assignments <- c("CheA","CheV","CheW","CheW","CheW","CheW")

results <- SimpLogo("sample_alignments/",
                    group.names = groups,
                    lineage.names = lineage.assignments)
#> 
```

The **results** object is a list, with ***results.table*** being a
long-format dataframe containing the converted frequencies and blended
colors in hexadecimal and ***color.scheme*** showing the chosen residue
groups and ideal color codes. This object can be used as input for
SimpLogoPlot().

``` r
head(results$results.table, n=10)
#>      color position     gap.freq info.content top.type secondary.type top.color
#> 1  #4e330b        1 0.9153339605   3.17106897  Q/N/S/T            D/E   #ffe119
#> 2  #366043        2 0.3066792098   2.83531858      A/G        M/V/L/I   #2F4F4F
#> 3  #5cad57        3 0.0799623706   2.76753716  M/V/L/I          F/Y/W   #3cb44b
#> 4  #7cb94d        4 0.0159924741   1.59028786  M/V/L/I        Q/N/S/T   #3cb44b
#> 5  #53b451        5 0.0047036689   3.02199594  M/V/L/I          F/Y/W   #3cb44b
#> 6  #a2908d        6 0.0037629351   0.83436446    H/K/R        Q/N/S/T   #4363d8
#> 7  #5cb355        7 0.0037629351   2.22179693  M/V/L/I            A/G   #3cb44b
#> 8  #796f5c        8 0.0028222013   2.13808286      A/G        Q/N/S/T   #2F4F4F
#> 9  #bf7d69        9 0.0028222013   1.19827800      D/E            A/G   #e6194b
#> 10 #dda359       10 0.0028222013   1.46895956  Q/N/S/T            D/E   #ffe119
#>    secondary.color   arch lineage residue     dummy
#> 1          #e6194b CheA.I    CheA       1 IC (bits)
#> 2          #3cb44b CheA.I    CheA       2 IC (bits)
#> 3          #9a6324 CheA.I    CheA       3 IC (bits)
#> 4          #ffe119 CheA.I    CheA       4 IC (bits)
#> 5          #9a6324 CheA.I    CheA       5 IC (bits)
#> 6          #ffe119 CheA.I    CheA       6 IC (bits)
#> 7          #2F4F4F CheA.I    CheA       7 IC (bits)
#> 8          #ffe119 CheA.I    CheA       8 IC (bits)
#> 9          #2F4F4F CheA.I    CheA       9 IC (bits)
#> 10         #e6194b CheA.I    CheA      10 IC (bits)

results$color.scheme
#>     H/K/R       D/E       A/G   Q/N/S/T   M/V/L/I     F/Y/W         P         C 
#> "#4363d8" "#e6194b" "#2F4F4F" "#ffe119" "#3cb44b" "#9a6324"  "purple" "#46f0f0" 
#>         - 
#>        NA

## here, the group names defined above are stored in the "arch" column
```

You may specify a custom color scheme for the residue groups (though
this is not recommended). To do so, you must include the “res.color”
parameter by specifying a named vector (names are taken from the default
residue groups). The vector must have all 9 groups included.

### Plotting SimpLogo

SimpLogoPlot() uses ggplot2 to generate the plot objects (with
geom_tiles).

``` r
library(ggplot2)
library(Cairo)

plot <- SimpLogoPlot(results, plot.ic = TRUE)
```

It produces a list of 3 graphical results, the primary SimpLogo plot
(the thing with the colored boxes) as a ggplot2 object, shown here:

``` r
plot$primary.plot
```

<img src="man/figures/README-fig2-1.png" width="100%" style="display: block; margin: auto;" />

A line graph of position-wise information content (in bits) as a ggplot2
object, shown here:

``` r
plot$ic.plot
```

<img src="man/figures/README-fig3-1.png" width="100%" style="display: block; margin: auto;" />

And a final combined plot (for convenience) created using the
plot_grid() function of cowplot, shown here:

``` r
plot$final.plot
```

<img src="man/figures/README-fig4-1.png" width="100%" height="120%" style="display: block; margin: auto;" />

### Manually editing individual graphs

You can edit the individual plots as ggplot2 objects by modifying the.
Here we add a new overall title to our plot and change our axis text to
green using theme \[notice how we used the primary.plot object, the
final.plot can’t be edited due to the use of plot_grid()\].

``` r
## now plot
plot$primary.plot + 
  labs(subtitle="this is an example title.") + 
  theme(axis.title.x.bottom =element_text(color="forestgreen"), 
        axis.text.x.bottom = element_text(color="forestgreen"), 
        axis.text.y.left = element_text(color="forestgreen"),
        plot.margin = margin(t=2)) ## to accomodate new title
```

<img src="man/figures/README-fig5-1.png" width="100%" style="display: block; margin: auto;" />

### Recombining plots

If you decide to edit the ggplot2 objects, you can match the plot_grid()
call in SimpLogoPlot() using the following:

``` r
final.plot <- cowplot::plot_grid(plot$ic.plot, 
                                 NULL, ## add blank space, to modify space between plots
                                 plot$primary.plot,
                                 ncol = 1, align = 'v',
                                 axis = "lr", rel_heights = c(0.35, 0.0225, 1.4),
                                 labels = c("A","","B")) ## don't forget blanks for NULL elements
final.plot
```

<img src="man/figures/README-fig6-1.png" width="100%" style="display: block; margin: auto;" />

### Renumbering alignment positions

If you’d like to easily renumber your alignment, starting at position 18
for example, you can use the following option when calling
SimpLogoPlot() \[or just modify it yourself by manipulating the ggplot2
object\]

``` r
plot2 <- SimpLogoPlot(results, plot.ic = TRUE, position.start = 18)
plot2$final.plot
```

<img src="man/figures/README-fig7-1.png" width="100%" style="display: block; margin: auto;" />

### Reordering sequence groups/architectures

You might want to reorder your architecture/sequence alignment groups.
You can set the factor levels in the original results object to manually
specify the order to ggplot2.

``` r
results2 <- results
results2$arch <- factor(results2$arch, levels = c("CheV.I","CheA.I","CheW.IB","CheW.IC","CheW.II.1","CheW.II.2"), ordered = TRUE) ## put CheV first, then CheA, then CheW
plot3 <- SimpLogoPlot(results2, plot.ic = TRUE, position.start = 18)
plot3$final.plot
```

<img src="man/figures/README-fig9-1.png" width="100%" style="display: block; margin: auto;" />

### Weirder examples

Perhaps you want to visualize only a specific type of residue. Maybe you
want to highlight only negative charged side chains. One way to do this
is the use a custom color palette. Here, we can set all residue groups
to “black” except D/E, which we make “red”. This generates an
heatmap-like SimpLogo.

``` r
highlight.negative.res.results <- SimpLogo("sample_alignments/",
                                           group.names = groups,
                                           lineage.names = lineage.assignments,
                                           res.colors = c("H/K/R" = "black", 
                                                          "D/E" = "red", 
                                                          "A/G" = "black", 
                                                          "Q/N/S/T" = "black",
                                                          "M/V/L/I" = "black", 
                                                          "F/Y/W" = "black", 
                                                          "P" = "black", 
                                                          "C" = "black", 
                                                          "-" = NA))
plot4 <- SimpLogoPlot(highlight.negative.res.results, plot.ic = TRUE)
plot4$final.plot
```

<img src="man/figures/README-fig10-1.png" width="100%" style="display: block; margin: auto;" />

## Output

The output of SimpLogo() is a dataframe with 10 fields:

- **color**: blended final color for main tile (in hexadecimal)
- **position**: position in alignment (starting from 1)
- **gap.freq**: frequency of gaps in alignment at position (1 = 100%)
- **info.content**: estimated information content (in bits) at position
  from alignment
- **top.type**: the most abundant residue type at the given position
- **secondary.type**: the second most abundant residue type at the given
  position
- **top.color**: idealized (100%) color for the most abundant residue
  type (in hexadecimal)
- **secondary.color**: idealized (100%) color for the second most
  abundant residue type (in hexadecimal)
- **arch**: sequence group assignment (taken from filenames if not s)
- **lineage**: higher level group assignment (useful for plotting)
- **residue**: residue number (used for plotting function)

The output of SimpLogoPlot() is a list containing 3 plot objects:

- **final.plot**: merged SimpLogo and IC plot (not a true ggplot2
  object)
- **primary.plot**: main SimpLogo plot (editable ggplot2 object)
- **ic.plot**: line graph of estimated information content (editable
  ggplot2 object; NULL if plot.ic = FALSE)

## To do list

Custom residue groupings. Add more complex instances in example section,
including highlighting specific positions, etc.

## License information

This package is licensed under GPL3 - please refer to LICENSE.md for
more details
