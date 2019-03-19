#!/usr/local/bin/Rscript

task <- dyncli::main()

# load libraries
library(dyncli, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)

library(destiny, warn.conflicts = FALSE)

#####################################
###           LOAD DATA           ###
#####################################
expression <- task$expression %>% as.matrix
parameters <- task$parameters
priors <- task$priors

start_cell <-
  if (!is.null(priors$start_id)) {
    sample(priors$start_id, 1)
  } else {
    NULL
  }

# TIMING: done with preproc
timings <- list(method_afterpreproc = Sys.time())

#####################################
###        INFER TRAJECTORY       ###
#####################################
# run diffusion maps
dm <- destiny::DiffusionMap(
  data = expression,
  sigma = parameters$sigma,
  distance = parameters$distance,
  n_eigs = parameters$ndim,
  density_norm = parameters$density_norm,
  n_local = parameters$n_local,
  vars = parameters$features_id
)

# run DPT
if (!is.null(start_cell)) {
  tips <- which(rownames(expression) %in% start_cell)
} else {
  tips <- destiny::random_root(dm)
}
dpt <- destiny::DPT(
  dm,
  w_width = parameters$w_width,
  tips = tips
)

# find DPT tips
tips <- destiny::tips(dpt)
tip_names <- rownames(expression)[tips]

# TIMING: done with trajectory inference
timings$method_aftermethod <- Sys.time()

#####################################
###     SAVE OUTPUT TRAJECTORY    ###
#####################################
cell_ids <- rownames(expression)

# construct grouping
grouping <- dpt@branch[,1] %>%
  ifelse(is.na(.), 0, .) %>%
  as.character() %>%
  paste0("Tip", .) %>%
  set_names(cell_ids)

group_ids <- sort(unique(grouping))

# collect distances from tips
tip_dists <- dpt[,tips] %>%
  magrittr::set_colnames(., paste0("Tip", seq_len(ncol(.)))) %>%
  magrittr::set_rownames(cell_ids)

# calculate progressions
outs <- map(
  group_ids,
  function(gid) {
    cat("Processing ", gid, "\n", sep = "")
    cixs <- which(grouping == gid)
    cids <- cell_ids[cixs]
    if (length(cids) > 0) {
      if (gid == "Tip0") {
        progr <- data_frame(
          cell_id = cids,
          from = gid,
          to = sample(setdiff(group_ids, gid), length(cids), replace = TRUE),
          percentage = 0
        )
        list(progr = progr, milnet = NULL)
      } else {
        # calculate min dist of gid to all other cells
        max_range <- min(tip_dists[-cixs, gid])

        # calculate percentage value of cells in gid
        percentage <- 1 - pmin(tip_dists[cixs, gid] / max_range, 1)
        progr <- data_frame(cell_id = cids, from = "Tip0", to = gid, percentage = percentage)
        milnet <- data_frame(from = "Tip0", to = gid, length = max_range, directed = FALSE)
        list(progr = progr, milnet = milnet)
      }
    } else {
      list()
    }
  }
)
progressions <- map_df(outs, ~ .$progr)
milestone_network <- map_df(outs, ~ .$milnet)

# collect dimred
dimred <- dm@eigenvectors %>%
  magrittr::set_colnames(., paste0("Comp", seq_len(ncol(.)))) %>%
  magrittr::set_rownames(cell_ids)

output <-
  wrap_data(
    cell_ids = cell_ids
  ) %>%
  add_grouping(
    group_ids = group_ids,
    grouping = grouping
  ) %>%
  add_trajectory(
    milestone_ids = group_ids,
    milestone_network = milestone_network,
    progressions = progressions
  ) %>% 
  add_dimred(
    dimred = dimred
  ) %>%
  add_timings(
    timings = timings
  )

dyncli::write_output(output, task$output)
