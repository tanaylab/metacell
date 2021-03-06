% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scmat.r
\name{scm_merge_mats}
\alias{scm_merge_mats}
\title{Merge multiple single cell matrix object.}
\usage{
scm_merge_mats(...)
}
\arguments{
\item{...}{tgScMat objects to merge.
Each parameter can be either a single tgScMat or a list of tgScMat (that will be merged)}
}
\description{
Merge multiple single cell matrix object.
Return the merged matrix, with merged meta data and issues an error if there are overlapping cell names between the two matrices.
}
\details{
In case genes sets differs between the matrices, the union is used, with zeros (not NAs!) filling up the missing genes in the respective matrix.
}
