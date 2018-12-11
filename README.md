---
output:
  md_document:
    variant: markdown_github
---



# metacell

The MetaCell R package facilitates analysis of single cell RNA-seq UMI matrices 
by computing partitions of a cell similarity graph into small (~20-200 typically) homogeneous groups of cells which are defined as metacells (MCs). The derived MCs are then used for building different representations of the data, allowing matrix or 2D graph visualization forming a basis for analysis of cell types, subtypes, transcriptional gradients, cell-cycle variation, gene modules and their regulatory models and more. More details on the usage of the MetaCell pipeline is available in the package vignettes, and in papers using it.

#### References:

Method: Baran et al. 2018 ([bioarxiv](https://www.biorxiv.org/content/early/2018/10/08/437665)).

Functions reference and usage vignettes are available in the package [homepage](https://tanaylab.bitbucket.io/metacell-r/index.html).

Examples of applications:

* [Giladi et al., Nat Cell Biol 2018](http://www.nature.com/articles/s41556-018-0121-4) (Mouse hematopoiesis scRNA-seq)
* [Sebe-pedros et al., NEE 2018](https://www.nature.com/articles/s41559-018-0575-6) (whole organisms scRNA-seq)
* [Sebe-pedros et al., Cell 2018](https://www.cell.com/cell/abstract/S0092-8674(18)30596-8) (whole organism scRNA-seq)
* [Bornstein et al., Nature 2018](https://www.nature.com/articles/s41586-018-0346-1) (thymic stroma scRNA-seq)
* [Cohen et al., Cell 2018](https://www.cell.com/cell/fulltext/S0092-8674(18)31181-4) (lung scRNA-seq)

R Package [homepage](https://tanaylab.bitbucket.io/metacell-r/index.html).

#### Installation


```r
install.packages('BiocManager') 
BiocManager::install('metacell',  site_repository = 'tanaylab.bitbucket.io/repo', update = FALSE)
```

**Note**: Metacell is implemented in R and C++. In particular it uses the Tanay group tgstat library that utilizes shared memory and distributed computing (as well as some specific optional CPU features). The package is tested on linux and macbooks, and is currently not compatible on Windows. A typical application will require at least 16G RAM. 
For heavier applications (100K cells) we recommend a dual CPU multi-core workstation with 128GM RAM or more.
