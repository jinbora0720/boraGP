# find geographically sensible neighbors
barrier_neighbor <- function(coords, coords.0 = NULL, ord = NULL,  
                             n.neighbors = 15, barrier, 
                             cores = 2, 
                             verbose = TRUE, 
                             debug = list(barrier_n.indx = NULL, 
                                          barrier_dist = NULL, 
                                          barrier_nn.indx.0_list = NULL, 
                                          barrier_dist0 = NULL, 
                                          ref_window = NULL, 
                                          nonref_window = NULL,
                                          ref_fill = TRUE, 
                                          nonref_fill = TRUE)) {
  
  # coords: reference coords
  # coords.0: other coords
  # ord: ordering for coords 
  # n.neighbors: number of neighbors 
  # barrier as an sf object
  
  #################
  # Prepare input #
  #################
  # set the number of cores
  doParallel::registerDoParallel(cores = cores)
  
  # make coords data frame
  coords <- data.frame(coords)
  coords.0 <- data.frame(coords.0)
  colnames(coords) = colnames(coords.0) <- 
    c("easting", "northing")
  
  # make coords as an sf object
  crs <- sf::st_crs(barrier)
  coords_sf <- st_as_sf(rbind(coords, coords.0), 
                        coords = c("easting", "northing"), 
                        crs = crs)
  
  # order coords
  n_tr <- nrow(coords)
  if (is.null(ord)) {
    ord <- order(coords[,1])
  }
  coords_tr_ord <- coords[ord,]
  
  if (!is.null(coords.0)) {
    n_tt <- nrow(coords.0)
    coords_sf_ord <- coords_sf$geometry[c(ord, n_tr+(1:n_tt))]
  } else {
    coords_sf_ord <- coords_sf$geometry[ord]
  }
  
  #####
  # R #
  #####
  # find first-order neighbors
  if (is.null(debug$barrier_n.indx)) {
    if (verbose) {
      cat("Finding first-order neighbors for reference locations... \n")
    }
    first_R <- 
      foreach (i = 2:n_tr) %dopar% {
        nn2res <- RANN::nn2(data = coords_tr_ord[1:(i-1),],                     # where
                            query = coords_tr_ord[i,],                          # whose
                            k = i-1)                                            # how many
        nindx <- as.numeric(nn2res$nn.idx)
        ndists <- as.numeric(nn2res$nn.dists)
        if (i < (n.neighbors+2)) {
          return(cbind(nindx, ndists))
        } else {
          barrier_nindx = barrier_d <- NULL
          j <- 1
          repeat {
            k <- nindx[j]
            test <- test_cross(barrier = barrier,
                               pt1 = coords_sf_ord[i],
                               pt2 = coords_sf_ord[k])
            if (is.na(test)) {
              barrier_nindx <- c(barrier_nindx, k)
              barrier_d <- c(barrier_d, ndists[j])
            } 
            
            if (length(barrier_nindx) == n.neighbors | j == i-1) break
            j <- j + 1
          }
          
          ## save
          return(cbind(barrier_nindx, barrier_d))
        }
      }
    
    barrier_n.indx = barrier_dist <- list()
    barrier_n.indx[[1]] = barrier_dist[[1]] <- NA
    barrier_n.indx[2:n_tr] <- lapply(first_R, function(x) as.numeric(x[,1]))
    barrier_dist[2:n_tr] <- lapply(first_R, function(x) as.numeric(x[,2]))
  } else {
    if (verbose) {
      cat("User specified BORA-GP neighbors for reference locations.\n")
    }
    barrier_n.indx <- debug$barrier_n.indx
    barrier_dist <- debug$barrier_dist
  }
  
  # check if any location needs fill-ins
  num_nb <- sapply(barrier_n.indx, length)
  im1 <- 1:n_tr-1
  im1[im1 > n.neighbors] <- n.neighbors
  
  if (sum(num_nb - im1 < 0) == 0) {
    if (verbose) {
      cat("BORA-GP neighbors for reference locations are complete.\n")
    }
  } else {
    if (debug$ref_fill) {
      needi_idx <- which(num_nb == 0)
      
      # find proxy of first-order neighbors for needi_idx 
      if (length(needi_idx) > 0) {
        if (verbose) {
          cat("Finding proxies of first-order neighbors for", 
              length(needi_idx), "reference location(s)...\n")
        }
        
        if (is.null(debug$ref_window)) {
          ref_window <- rep(n.neighbors, n_tr)
        } else {
          ref_window <- debug$ref_window
        }
        
        i_R <- foreach (i = needi_idx) %dopar% {
          ## make grid using n.neighbors nearest neighbors
          nn2res <- RANN::nn2(data = coords_tr_ord[1:(i-1),],                     # where
                              query = coords_tr_ord[i,],                          # whose
                              k = n.neighbors)                                    # how many
          nindx <- as.numeric(nn2res$nn.idx)
          ndists <- as.numeric(nn2res$nn.dists)
          
          ## grid size is the shortest overlapped length
          pt1 <- coords_sf_ord[i]
          gridsize <- min(sapply(nindx, function(j) cross_length(
            barrier = barrier, pt1 = pt1, pt2 = coords_sf_ord[j])))
          
          ## the half length of the grid box is the longest distance 
          ## from n.neighbors nearest neighbors
          radius <- max(ndists)
          gridbox <- st_bbox(pt1) + c(-radius, -radius, radius, radius)
          igrid <- expand.grid(easting = seq(gridbox$xmin, 
                                             gridbox$xmax, by = gridsize), 
                               northing = seq(gridbox$ymin, 
                                              gridbox$ymax, by = gridsize))
          igrid_sf <- st_as_sf(igrid, coords = c("easting", "northing"), 
                               crs = crs)
          n_igrid <- nrow(igrid)
          
          ## sort grid by shortest Euclidean distance 
          nn2_g <- RANN::nn2(data = igrid,                                        # where
                             query = coords_tr_ord[i,],                           # whose
                             k = n_igrid)                                         # how many
          nindx_g <- as.numeric(nn2_g$nn.idx)
          ndists_g <- as.numeric(nn2_g$nn.dists)
          
          barrier_nindx = barrier_d <- NULL
          j_g <- 1
          repeat {
            k_g <- nindx_g[j_g]
            test_g <- test_cross(barrier = barrier,
                                 pt1 = pt1,
                                 pt2 = igrid_sf$geometry[k_g])
            if (is.na(test_g)) {
              nn2res_gtoR <- RANN::nn2(data = coords_tr_ord[1:(i-1),],            # where
                                       query = igrid[k_g,],                       # whose
                                       k = ref_window[i])                         # how many
              nindx_gtoR <- as.numeric(nn2res_gtoR$nn.idx)
              ndists_gtoR <- as.numeric(nn2res_gtoR$nn.dists)
              
              pt1_gtoR <- igrid_sf$geometry[k_g]
              test_gtoR <- sapply(nindx_gtoR, function(j) 
                test_cross(barrier = barrier, pt1 = pt1_gtoR, 
                           pt2 = coords_sf_ord[j]))
              connect_gtoR <- which(is.na(test_gtoR))
              add <- min(length(connect_gtoR), n.neighbors)
              
              if (add > 0) {
                barrier_nindx <- c(barrier_nindx, nindx_gtoR[connect_gtoR[1:add]])
                barrier_d <- c(barrier_d, 
                               ndists_gtoR[connect_gtoR[1:add]] + ndists_g[j_g])
              }
            }
            if (length(barrier_nindx) > 0 | j_g == n_igrid) break
            j_g <- j_g+1
          }
          
          ## if still isolated, take the nearest neighbor as a proxy
          if (is.null(barrier_nindx)) {
            barrier_nindx <- nindx[1]
            barrier_d <- ndists[1]
            return(cbind(barrier_nindx, barrier_d, i))
          } else {
            return(cbind(barrier_nindx, barrier_d, NA))
          }
        }
        
        barrier_n.indx[needi_idx] <- lapply(i_R, function(x) as.numeric(x[,1]))
        barrier_dist[needi_idx] <- lapply(i_R, function(x) as.numeric(x[,2]))
        nn_idx <- na.omit(sapply(i_R, function(x) as.numeric(x[,3])[1]))
        
        # message if some reference locations were not expandable
        if (length(nn_idx) > 0) {
          message("!Caution! ", length(nn_idx), 
                  " reference location(s), ", paste0(nn_idx, ", "), 
                  "in the ordered coords failed to find proxies for ",
                  "first-order neighbors using imaginary grid. ", 
                  "The nearest neighbor is taken as the proxy.")
          # Recommend to increase the number of neighbors to consider at each grid point.
        }
      }
      
      # check if any location needs second-order neighbors
      num_nb <- sapply(barrier_n.indx, length)
      need2_idx <- which(num_nb - im1 < 0)
      
      if (length(need2_idx) > 0) {
        if (verbose) {
          cat("Finding second-order neighbors for", length(need2_idx), 
              "reference location(s)... \n")
        }
        
        res_hn_R <- find_hn_R(need_idx = need2_idx, 
                              barrier_n.indx = barrier_n.indx, 
                              barrier_dist = barrier_dist, 
                              coords_tr_ord = coords_tr_ord, 
                              n.neighbors = n.neighbors)
        barrier_n.indx <- res_hn_R$barrier_n.indx   
        barrier_dist <- res_hn_R$barrier_dist
        
        # check if any location needs third-order neighbors
        num_nb <- sapply(barrier_n.indx, length)
        need3_idx <- which(num_nb - im1 < 0)
        
        if (length(need3_idx) > 0) {
          warning(length(need3_idx), " reference locations, ",
                  paste0(need3_idx, ", "), "in the ordered coords needs", 
                  " third-order neighbors.\n")
        } else {
          if (verbose) 
            cat("BORA-GP neighbors for reference locations are complete.\n")
        }
      } else {
        if (verbose) 
          cat("BORA-GP neighbors for reference locations are complete.\n")
      }
    }
  }
  
  #####
  # U #
  #####
  if (is.null(coords.0)) {
    out <- list(barrier_n.indx = barrier_n.indx,
                barrier_dist = barrier_dist)
  } else {
    # find first-order neighbors
    if (is.null(debug$barrier_nn.indx.0_list)) {
      if (verbose) {
        cat("Finding first-order neighbors for non-reference locations... \n")
      }
      first_U <- 
        foreach (i = 1:n_tt) %dopar% {
          nn2res <- RANN::nn2(data = coords,                                    # where
                              query = coords.0[i,],                             # whose
                              k = n_tr)                                         # how many
          nindx <- as.numeric(nn2res$nn.idx)
          ndists <- as.numeric(nn2res$nn.dists)
          barrier_nn0 = barrier_d0 <- NULL
          j <- 1
          repeat {
            k <- nindx[j]
            test <- test_cross(barrier = barrier,
                               pt1 = coords_sf$geometry[n_tr+i],
                               pt2 = coords_sf$geometry[k])
            if (is.na(test)) {
              barrier_nn0 <- c(barrier_nn0,k)
              barrier_d0 <- c(barrier_d0,ndists[j])
            }
            
            if (length(barrier_nn0) == n.neighbors | j == n_tr) break
            j <- j + 1
          }
          
          ## save
          return(cbind(barrier_nn0, barrier_d0))
        }
      barrier_nn.indx.0_list <- lapply(first_U, function(x) as.numeric(x[,1]))
      barrier_dist0 <- lapply(first_U, function(x) as.numeric(x[,2]))
    } else {
      if (verbose) {
        cat("User specified BORA-GP neighbors for non-reference locations.\n")
      }
      barrier_nn.indx.0_list <- debug$barrier_nn.indx.0_list
      barrier_dist0 <- debug$barrier_dist0
    }
    
    # check if any location needs fill-ins
    num_nb.0 <- sapply(barrier_nn.indx.0_list, length)
    
    if (sum(num_nb.0 != n.neighbors) == 0) {
      if (verbose) {
        cat("BORA-GP neighbors for non-reference locations are complete.\n")
      }
    } else {
      if (debug$nonref_fill) {
        # check if any location needs imaginary neighbors
        needi_idx.0 <- which(num_nb.0 == 0)
        
        if (length(needi_idx.0) > 0) {
          if (verbose) {
            cat("Finding proxies of first-order neighbors for", 
                length(needi_idx.0), "non-reference location(s)...\n")
          }
          
          if (is.null(debug$nonref_window)) {
            nonref_window <- rep(n.neighbors, n_tt)
          } else {
            nonref_window <- debug$nonref_window
          }
          
          i_U <- foreach (i = needi_idx.0) %dopar% {
            ## make grid using n.neighbors nearest neighbors
            nn2res <- RANN::nn2(data = coords,                                    # where
                                query = coords.0[i,],                             # whose
                                k = n.neighbors)                                  # how many    
            nindx <- as.numeric(nn2res$nn.idx)
            ndists <- as.numeric(nn2res$nn.dists)
            
            ## grid size is the shortest overlapped length
            pt1 <- coords_sf$geometry[n_tr+i]
            gridsize <- min(sapply(nindx, function(j) cross_length(
              barrier = barrier, pt1 = pt1, pt2 = coords_sf$geometry[j])))
            
            ## the half length of the grid box is the longest distance 
            ## from n.neighbors nearest neighbors
            radius <- max(ndists)
            gridbox <- st_bbox(pt1) + c(-radius, -radius, radius, radius)
            igrid <- expand.grid(easting = seq(gridbox$xmin, 
                                               gridbox$xmax, by = gridsize), 
                                 northing = seq(gridbox$ymin, 
                                                gridbox$ymax, by = gridsize))
            igrid_sf <- st_as_sf(igrid, coords = c("easting", "northing"), 
                                 crs = crs)
            n_igrid <- nrow(igrid)
            
            ## sort grid by shortest Euclidean distance 
            nn2_g <- RANN::nn2(data = igrid,                                      # where
                               query = coords.0[i,],                              # whose
                               k = n_igrid)                                       # how many
            nindx_g <- as.numeric(nn2_g$nn.idx)
            ndists_g <- as.numeric(nn2_g$nn.dists)
            
            barrier_nindx = barrier_d <- NULL
            j_g <- 1
            repeat {
              k_g <- nindx_g[j_g]
              test_g <- test_cross(barrier = barrier,
                                   pt1 = pt1,
                                   pt2 = igrid_sf$geometry[k_g])
              if (is.na(test_g)) {
                nn2res_gtoR <- RANN::nn2(data = coords,                           # where
                                         query = igrid[k_g,],                     # whose
                                         k = nonref_window[i])                    # how many
                nindx_gtoR <- as.numeric(nn2res_gtoR$nn.idx)
                ndists_gtoR <- as.numeric(nn2res_gtoR$nn.dists)
                
                pt1_gtoR <- igrid_sf$geometry[k_g]
                test_gtoR <- sapply(nindx_gtoR, function(j) 
                  test_cross(barrier = barrier, pt1 = pt1_gtoR, 
                             pt2 = coords_sf$geometry[j]))
                connect_gtoR <- which(is.na(test_gtoR))
                add <- min(length(connect_gtoR), n.neighbors)
                
                if (add > 0) {
                  barrier_nindx <- c(barrier_nindx, nindx_gtoR[connect_gtoR[1:add]])
                  barrier_d <- c(barrier_d, 
                                 ndists_gtoR[connect_gtoR[1:add]] + ndists_g[j_g])
                }
              }
              if (length(barrier_nindx) > 0 | j_g == n_igrid) break
              j_g <- j_g+1
            }
            
            ## if still isolated, take the nearest neighbor as a proxy
            if (is.null(barrier_nindx)) {
              barrier_nindx <- nindx[1]
              barrier_d <- ndists[1]
              return(cbind(barrier_nindx, barrier_d, i))
            } else {
              return(cbind(barrier_nindx, barrier_d, NA))
            }
          }
          
          barrier_nn.indx.0_list[needi_idx.0] <- 
            lapply(i_U, function(x) as.numeric(x[,1]))
          barrier_dist0[needi_idx.0] <- lapply(i_U, function(x) as.numeric(x[,2]))
          nn_idx.0 <- na.omit(sapply(i_U, function(x) as.numeric(x[,3])[1]))
          
          # message if some reference locations were not expandable
          if (length(nn_idx.0) > 0) {
            message("!Caution! ", length(nn_idx.0), 
                    " non-reference location(s), ", paste0(nn_idx.0, ", "), 
                    "in coords.0 failed to find proxies for ",
                    "first-order neighbors using imaginary grid. ", 
                    "The nearest neighbor is taken as the proxy.")
          }
        } 
        
        # check if any location needs second-order neighbors
        num_nb.0 <- sapply(barrier_nn.indx.0_list, length)
        need2_idx.0 <- which(num_nb.0 != n.neighbors)
        
        if (length(need2_idx.0) > 0) {
          if (verbose) {
            cat("Finding second-order neighbors for", length(need2_idx.0), 
                "non-reference location(s)... \n")
          }
          
          res_hn_U <- find_hn_U(need_idx = need2_idx.0, 
                                barrier_nn.indx.0_list = barrier_nn.indx.0_list, 
                                barrier_dist0 = barrier_dist0,
                                barrier_n.indx = barrier_n.indx, 
                                barrier_dist = barrier_dist, 
                                ord = ord, 
                                coords.0 = coords.0, 
                                barrier = barrier, 
                                coords_sf_ord = coords_sf_ord,
                                coords_tr_ord = coords_tr_ord, 
                                n.neighbors = n.neighbors)
          
          barrier_nn.indx.0_list <- res_hn_U$barrier_nn.indx.0_list   
          barrier_dist0 <- res_hn_U$barrier_dist0  
          
          # check if any location needs third-order neighbors
          num_nb.0 <- sapply(barrier_nn.indx.0_list, length)
          need3_idx.0 <- which(num_nb.0 != n.neighbors)
          
          if (length(need3_idx.0) > 0) {
            warning(length(need3_idx.0), " non-reference locations, ",
                    paste0(need3_idx.0, ", "), "in coords.0 needs", 
                    " third-order neighbors.\n")
          } else {
            if (verbose) 
              cat("BORA-GP neighbors for non-reference locations are complete.\n")
          }
        } else {
          if (verbose) 
            cat("BORA-GP neighbors for non-reference locations are complete.\n")
        }
      }
    }
    
    out <- list(barrier_n.indx = barrier_n.indx,
                barrier_dist = barrier_dist, 
                barrier_nn.indx.0_list = barrier_nn.indx.0_list, 
                barrier_dist0 = barrier_dist0)
  }
  return(out)
}