metacell
========

The MetaCell package facilitate analysis of single cell RNA-seq UMI matrices through generation of cell-to-cell similarity graphs and computation of partitions of the graph into small (~20-200 typically) homogeneous groups of cells. The derived metacells are then used for deriving different representations of the data, and form the basis for analysis of cell types, subtypes, transcriptional gradients, cell-cycle variation, gene modules and their regulatory models and more. More details on the usage of the MetaCell pipeline is available in the package vignettes, and in papers using it.

#### References:

Method: Baran et al. 2018 (bioarxiv link comming up)

Examples of applications:

-   [Giladi et al., Nat Cell Biol 2018](http://www.nature.com/articles/s41556-018-0121-4) (Mouse hematopoiesis scRNA-seq)
-   [Sebe-pedros et al., NEE 2018](https://www.nature.com/articles/s41559-018-0575-6) (whole organisms scRNA-seq)
-   [Sebe-pedros et al., Cell 2018](https://www.cell.com/cell/abstract/S0092-8674(18)30596-8) (whole organism scRNA-seq)

#### Installation

``` r
install.packages('metacell', repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo'))
```
