########################################################################
#' R6 class representing a cell atlas
#'
#' @description
#'
#' @export
tgMAtlas <- R6::R6Class("tgMAtlas", cloneable = FALSE, public = list(
    #' @field mat UMI matrix (class tgScMat)
    #' @field mc Metacells (class tgMCCov)
    #' @field gset Gene set (class tgGeneSets)
    #' @field mc2d Cells/metacells 2D projection (class tgMC2D)
    #' @field ecdfs ECDFs of mathces of ranking variables
    mat = NULL,
    mc = NULL,
    gset = NULL,
    mc2d = NULL,
    ecdfs = list(),

    #' @description Create a new tgMAtlas object
    #'
    #' @details
    #' It is possible to load an atlas from the filesystem by passing as
    #' `mat` the path of either an atlas json file or a directory holding
    #' an atals. In the later case it is possible to pass a glob
    #' pattern in `mc` that will be used to match files in the directory.
    #'
    #' @param mat,mc,gset,mc2d
    #' Objects or filenames of stored objects (.Rda or .rds) that will
    #' be used for member initalization. XXX
    #'
    #' @param rankings
    #' A vector of metadata columns (from `mat` metdata) that can be
    #' used for ranking projected cells.
    #'
    #' @param precalc_ecdf
    #' Whether to immediately calculated matches ECDFs for the columns in
    #' `ranking`. If this value is FALSE, the ECDFs will be caclulated the
    #' first time they are required.
    #'
    initialize = function(mat, mc = NULL, gset = NULL, mc2d = NULL, rankings = NULL, precalc_ecdf = FALSE) {
        if (is.null(gset) && is.null(mc2d)) {
            path <- mat
            fname <- mc
            config <- matlas_load_path(path, fname)
            fnames <- matlas_norm_fnames(config$fnames, path)
            if (!is.null(fnames)) {
                mat <- fnames$mat
                mc <- fnames$mc
                gset <- fnames$gset
                mc2d <- fnames$mc2d
            }
            if (is.null(rankings)) {
                rankings <- unlist(config$rankings)
            }
        }

        mandatory <- list(mat = mat, mc = mc, gset = gset)
        if (any(sapply(mandatory, is.null))) {
            stop(sprintf("The following members must be defined: %s", paste0(names(mandatory), collapse = ", ")))
        }

        self$mat <- matlas_load_object(mat, "tgScMat")
        self$mc <- matlas_load_object(mc, "tgMCCov")
        self$gset <- matlas_load_object(gset, "tgGeneSets")
        self$mc2d <- matlas_load_object(mc2d, "tgMC2D")

        if (!is.null(rankings)) {
            if (precalc_ecdf) {
                for (ranking in rankings) {
                    self$ecdfs[[ranking]] <- self$.ranking_ecdf(ranking)
                }
            }
            else {
                for (ranking in rankings) {
                    self$ecdfs[ranking] <- list(NULL)
                }
            }
        }

        invisible(self)
    },


    #' @description Save the atlas into a directory
    #'
    #' @param dname
    #' Name of directory that will hold the atals. If needed, this
    #' directory will be created if needed.
    #'
    #' @param fname
    #' Filename pattern for the files that will be stored (e.g.
    #' if `fname` is `"xyz"` the UMI matrix will be saved as
    #' `mat.xyz.rds`)
    save = function(dname, fname = "atlas") {
        if (!dir.exists(dname)) {
            dir.create(dname, recursive = TRUE)
        }

        objs <- list(
            mat = self$mat,
            mc = self$mc,
            gset = self$gset,
            mc2d = self$mc2d
        )

        fnames <- sapply(names(objs), function(name) {
            if (is.null(objs[[name]])) {
                return(NULL)
            }
            path <- paste0(name, ".", fname, ".rds")
            saveRDS(objs[[name]], file.path(dname, path))
            return(file.path(".", path))
        })

        fnames <- fnames[!sapply(fnames, is.null)]
        config <- list(fnames = fnames)
        config$rankings <- I(names(self$ecdfs))

        path <- paste0(fname, ".json")
        jsonlite::write_json(config, file.path(dname, path), null = "null", auto_unbox = TRUE, pretty = TRUE)
    },


    #' @description For given query cells, find nearest atlas matches
    #'
    #' @param query A genes x cells matrix of UMIs
    #' 
    #' @return
    #' A vector with the names of the matching atlas cells. The length
    #' of this vector equals the number of columns in `query`.
    project = function(query) {
        genes_ref <- intersect(names(self$gset@gene_set), rownames(self$mc@e_gc))
        genes <- intersect(genes_ref, rownames(query))
        if (length(genes) / length(genes_ref) < 0.5) {
            warning("Less than half of the atlas feature genes can be founf in query.")
        }

        ref_feats <- self$mat@mat[genes, names(self$mc@mc)]
        feats <- query[genes, ]

        ref_only <- setdiff(colnames(ref_feats), colnames(feats))
        if (length(ref_only) < ncol(ref_feats)) {
            ref_feats <- ref_feats[, ref_only]
        }

        ref_cor <- tgs_cor(as.matrix(feats), as.matrix(ref_feats))
        ref_cor[is.na(ref_cor)] <- -1

        best_ref <- colnames(ref_feats)[max.col(ref_cor, "first")]

        return(best_ref)
    },


    #' @description
    #' Annotate an atlas projection
    #' 
    #' @param projection
    #' A set of names of atlas cells, as returned by `project()`.
    #' 
    #' @param column
    #' Metadata column that will be used as the annotation
    #' 
    #' @return
    #' A vector with the annotation values
    annotate = function(projection, column) {
        match <- self$mat@cell_metadata[projection, column]
        return(match)
    },


    #' @description
    #' Rank a projection of cells using a ranking variable
    #'
    #' @param projection
    #' A set of names of atlas cells, as returned by `project()`.
    #'
    #' @param column
    #' Metadata column that will be used as the ranking variable.
    #'
    #' @details
    #' The ECDF of the ranking varibales for `projection` is
    #' calculated. Then the correlations between this ECDF and all
    #' ECDFs achieved by projecting the atlas upon itself are returned
    #'
    #' @return
    #' A vector with the correlation values
    rank = function(projection, column) {

        ecdf <- self$ecdfs[[column]]
        if (is.null(ecdf)) {
            self$ecdfs[[column]] <- self$.ranking_ecdf(column)
            ecdf <- self$ecdfs[[column]]
        }

        match <- self$mat@cell_metadata[projection, column]
        match <- tabulate(match, nbins = as.integer(rownames(ecdf)[nrow(ecdf)]))
        match <- cumsum(match)
        names(match) <- as.character(seq_along(match))

        match <- match[rownames(ecdf)]
        dim(match) <- c(length(match), 1)
        match <- tgs_cor(ecdf, match)

        dim(match) <- NULL
        names(match) <- colnames(ecdf)

        return(match)
    },


    #' @description
    #' PRIVATE: Calcuate the ECDF of a ranking column

    #' @param column
    #' Metadata column that will be used as the ranking variable.
    .ranking_ecdf = function(column) {
        mat <- self$mat
        mc <- self$mc
        gset <- self$gset

        genes <- names(gset@gene_set)
        feats <- mat@mat[genes, names(mc@mc)]

        ranking <- mat@cell_metadata[colnames(feats), column]

        # delta equals 1 when the two cells have a different ranking
        delta <- 1 - diag(1, max(ranking))
        delta <- delta[ranking, ]
        delta <- delta[, ranking]

        cross_ref <- tgstat::tgs_cor(as.matrix(feats)) * delta
        rank_match <- ranking[max.col(cross_ref, "first")]
        rank_match <- table(ranking, rank_match)

        match_target <- t(matrixStats::rowCumsums(rank_match))
        rownames(match_target) <- colnames(rank_match)
        colnames(match_target) <- rownames(rank_match)

        return(match_target)
    }
))


########################################################################
matlas_load_path <- function(path, fname = NULL) {

    if (is.null(fname)) {
        fname <- "*"
    }

    if (!file.exists(path)) {
        return(NULL)
    }

    if (!dir.exists(path)) {
        if (!is.null(fname)) {
            warning("Atlas loaded from json file. Parameter 'fname' will be ignored.")
        }
        return(jsonlite::read_json(path))
    }

    config <- Sys.glob(file.path(path, paste0(fname, ".json")))
    if (length(config) > 1) {
        stop("More than one match for json file:\n", paste0(config, collapse = ","))
    }
    if (length(config) == 1) {
        return(jsonlite::read_json(config))
    }

    mat <- basename(Sys.glob(file.path(path, paste0("mat.", fname, "*"))))
    mc <- basename(Sys.glob(file.path(path, paste0("mc.", fname, "*"))))
    gset <- basename(Sys.glob(file.path(path, paste0("gset.", fname, "*"))))
    mc2d <- basename(Sys.glob(file.path(path, paste0("mc2d.", fname, "*"))))

    fnames <- list(mat = mat, mc = mc, gset = gset, mc2d = mc2d)
    mask <- sapply(fnames, length) > 0
    fnames <- fnames[mask]

    mandatory <- c("mat", "mc", "gset")
    mask <- !(mandatory %in% names(fnames))
    if (any(mask)) {
        stop(sprintf("Could not find matching file for members: %s", paste0(mandatory[mask], collapse = ",")))
    }

    mask <- sapply(fnames, length) > 1
    if (any(mask)) {
        dups <- sapply(fnames[mask], paste0, collapse = ", ")
        dups <- paste0("[", dups, "]", collapse = "\n")
        stop("More than one match for atlas file:\n", dups)
    }

    return(list(fnames = fnames))
}


########################################################################
matlas_norm_fnames <- function(fnames, path) {
    if (!dir.exists(path)) {
        path <- dirname(path)
    }
    fnames <- lapply(fnames, function(fn) {
        prefix <- ifelse(startsWith(fn, "/"), "", path)
        return(file.path(prefix, fn))
    })
    return(fnames)
}


########################################################################
matlas_load_object <- function(obj, expected_class = NULL) {
    if (is.null(obj)) {
        return(NULL)
    }

    if (!is.null(expected_class) && (expected_class %in% class(obj))) {
        return(obj)
    }

    if (!is.character(obj)) {
        if (is.null(expected_class)) {
            return(obj)
        }
        stop(sprintf("Wrong 'obj' class: '%s', Expected: '%s'", paste0(class(obj), collapse = ","), expected_class))
    }

    if (!file.exists(obj)) {
        stop(sprintf("Could not read 'obj'. Missing file: '%s'", obj))
    }

    if (endsWith(obj, ".Rda")) {
        pocket <- new.env()
        load(obj, pocket)
        if (length(pocket) == 0) {
            stop(sprintf("No object is found in '%s'", obj))
        }
        if (length(pocket) > 1) {
            stop(sprintf("More than one object is found in '%s'", obj))
        }
        obj <- pocket[[names(pocket)[1]]]
    }
    else {
        obj <- readRDS(obj)
    }

    if (!is.null(expected_class) && !(expected_class %in% class(obj))) {
        stop(sprintf("Wrong 'obj' class: '%s', Expected: '%s'", paste0(class(obj), collapse = ","), expected_class))
    }

    return(obj)
}


# ########################################################################
# matlas_ranking_ecdf <- function(atlas, column) {

#     mat <- atlas$mat
#     mc <- atlas$mc
#     gset <- atlas$gset

#     genes <- names(gset@gene_set)
#     feats <- mat@mat[genes, names(mc@mc)]

#     ranking <- mat@cell_metadata[colnames(feats), column]

#     # delta equals 1 when the two cells have a different ranking
#     delta <- 1 - diag(1, max(ranking))
#     delta <- delta[ranking, ]
#     delta <- delta[, ranking]

#     cross_ref <- tgstat::tgs_cor(as.matrix(feats)) * delta
#     rank_match <- ranking[max.col(cross_ref, "first")]
#     rank_match <- table(ranking, rank_match)

#     match_target <- t(matrixStats::rowCumsums(rank_match))
#     rownames(match_target) <- colnames(rank_match)
#     colnames(match_target) <- rownames(rank_match)

# 	return(match_target)
# }


# ########################################################################
# matlas_rank <- function(self, projection, column) {
#         # best_ref <- colnames(ref_feats)[max.col(ref_cor, "first")]
#     	# best_ref_type = self$mc@colors[self$mc@mc[best_ref]]

#         ecdf <- self$ecdfs[[column]]
#         if (is.null(ecdf)) {
#             self$ecdfs[[column]] <- matlas_ranking_ecdf(self, column)
#             ecdf <- self$ecdfs[[column]]
#         }

#         match <- self$mat@cell_metadata[projection, column]
#         match <- tabulate(match, nbins = as.integer(rownames(ecdf)[nrow(ecdf)]))
#         match <- cumsum(match)
#         names(match) <- as.character(seq_along(match))

#         match <- match[rownames(ecdf)]
#         dim(match) <- c(length(match), 1)
#         match <- tgs_cor(ecdf, match)

#         dim(match) <- NULL
#         names(match) <- colnames(ecdf)

#         return(match)
# }
