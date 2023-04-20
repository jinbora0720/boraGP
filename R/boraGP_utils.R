# util functions 
# find higher-order neighbor for reference locations
find_hn_R <- function(need_idx, barrier_n.indx, barrier_dist, 
                      coords_tr_ord, n.neighbors) {
  for (i in need_idx) {                                                         # cannot parallelize
    inneed <- n.neighbors - length(barrier_n.indx[[i]])
    
    idx1 <- barrier_n.indx[[i]]                                                 # first-order neighbors
    idx2_all <- unique(unlist(barrier_n.indx[idx1]))                          
    idx2 <- idx2_all[!idx2_all %in% idx1]                                       # second-order neighbors
    
    if (length(idx2) == 0) stop(paste0("The neighbor set of the ", i, "th 
    reference location in the ordered coords cannot expand. 
    Recommend to examine the location on a map. Different ordering may help."))
    
    geodist <- data.frame(idx2 = unlist(barrier_n.indx[idx1]), 
                          dist = rep(barrier_dist[[i]], 
                                     times = sapply(barrier_n.indx[idx1], 
                                                    length)) +
                            unlist(barrier_dist[idx1])) %>% 
      group_by(idx2) %>% 
      summarize(geodist = min(dist)) %>% 
      na.omit()                                                                 # there may be NA if 1 in idx1
    
    idx2 <- na.omit(idx2)                                                       # delete NA in idx2 after computing geodist
    nn2res <- RANN::nn2(data = coords_tr_ord[idx2,],                            # where
                        query = coords_tr_ord[i,],                              # whose
                        k = length(idx2))                                       # how many
    Eucdist <- data.frame(idx2 = idx2[as.numeric(nn2res$nn.idx)], 
                          Eucdist = as.numeric(nn2res$nn.dists))
    
    pool <- right_join(geodist, Eucdist, by = "idx2") %>% 
      mutate(abs_diff = abs(geodist - Eucdist)) %>% 
      arrange(abs_diff, geodist) 
    
    idx2_pool <- pool %>% pull(idx2)
    dist_pool <- pool %>% pull(geodist)
    add <- min(length(idx2_pool), inneed)
    
    barrier_n.indx[[i]] <- c(idx1, idx2_pool[1:add])
    barrier_dist[[i]] <- c(barrier_dist[[i]], dist_pool[1:add])
  }
  
  return(list(barrier_n.indx = barrier_n.indx, 
              barrier_dist = barrier_dist))
}

# find higher-order neighbor for non-reference locations
find_hn_U <- function(need_idx, barrier_nn.indx.0_list, barrier_dist0,
                      barrier_n.indx, barrier_dist, 
                      ord, coords.0, barrier, 
                      coords_sf_ord, coords_tr_ord, 
                      n.neighbors) {
  
  n_tr <- nrow(coords_tr_ord)
  
  second0 <- foreach (i = need_idx) %dopar% {                                   
    idx1 <- barrier_nn.indx.0_list[[i]]                                         # first-order neighbors
    idx1_ord <- order(ord)[idx1]                                                # convert in ord
    
    m1 <- length(idx1)
    inneed <- n.neighbors - m1
    
    idx2_ord_list <- barrier_n.indx[idx1_ord]
    idx2_dist <- barrier_dist[idx1_ord]
    for (ii in 1:m1) {
      iii <- idx1_ord[ii]
      if (iii < n_tr) {
        idx2_which <- sapply(barrier_n.indx[(iii+1):n_tr], function(x)
          ifelse(sum(x == iii) == 0, 0, which(x == iii)))
        idx2_after <- which(idx2_which > 0) + iii
        if (length(idx2_after) > 0) {
          idx2_ord_list[[ii]] <- c(idx2_ord_list[[ii]], idx2_after)
          idx2_dist[[ii]] <- c(idx2_dist[[ii]], 
                               unlist(sapply((iii+1):n_tr, function(j) 
                                 barrier_dist[[j]][idx2_which[j-iii]])))
        }
      }
    }
    geodist <- data.frame(idx2_ord = unlist(idx2_ord_list), 
                          dist = rep(barrier_dist0[[i]], 
                                     times = sapply(idx2_ord_list, length)) +
                            unlist(idx2_dist)) %>% 
      group_by(idx2_ord) %>% 
      summarize(geodist = min(dist)) %>% 
      na.omit()                                                                 # there may be NA if 1 in idx1
    
    idx2_ord_all <- unique(unlist(idx2_ord_list))                          
    idx2_ord <- idx2_ord_all[!idx2_ord_all %in% idx1_ord]                       # second-order neighbors
    idx2_ord <- na.omit(idx2_ord)                                               # delete NA in idx2 after computing geodist
    
    nn2res <- RANN::nn2(data = coords_tr_ord[idx2_ord,],                        # where
                        query = coords.0[i,],                                   # whose
                        k = length(idx2_ord))                                   # how many
    Eucdist <- data.frame(idx2_ord = idx2_ord[as.numeric(nn2res$nn.idx)], 
                          Eucdist = as.numeric(nn2res$nn.dists))
    
    pool <- right_join(geodist, Eucdist, by = "idx2_ord") %>% 
      mutate(abs_diff = abs(geodist - Eucdist)) %>% 
      arrange(abs_diff, geodist) 
    
    idx2_ord_pool <- pool %>% pull(idx2_ord)
    dist_pool <- pool %>% pull(geodist)
    add <- min(length(idx2_ord_pool), inneed)
    
    return(cbind(c(idx1, ord[idx2_ord_pool[1:add]]),
                 c(barrier_dist0[[i]], dist_pool[1:add])))
  }
  
  barrier_nn.indx.0_list[need_idx] <- 
    lapply(second0, function(x) as.numeric(x[,1]))
  barrier_dist0[need_idx] <- 
    lapply(second0, function(x) as.numeric(x[,2]))
  
  return(list(barrier_nn.indx.0_list = barrier_nn.indx.0_list, 
              barrier_dist0 = barrier_dist0))
}

cross_length <- function(barrier, pt1, pt2) {
  line <- st_cast(st_combine(c(pt1, pt2)), "LINESTRING")
  cross_length <- st_length(st_intersection(line, st_buffer(barrier, 0)))
  return(ifelse(length(cross_length) == 0, 0, sum(cross_length)))
}

test_cross <- function(barrier, pt1, pt2) {
  line <- st_cast(st_combine(c(pt1, pt2)), "LINESTRING")
  return(as.numeric(st_intersects(line, barrier)))
}

create_Ctilde <- function(coords, neighbor.info, sig_sq, phi, nu = NULL,
                          base_cov = "exponential", tol = 0, cores = 2) {
  # coords: coords for reference locations
  # neighbor.info: list containing ord and n.indx
  # sig_sq: spatial variance
  # phi: decay parameter
  # nu: smoothness parameter, applicable only for base_cov = "matern"
  # base_cov = "exponential" or "matern"
  # tol: tolerance for solve
  # cores: number of cores for parallelization

  # set the number of cores
  registerDoParallel(cores = cores)

  ord <- neighbor.info$ord
  n.indx <- neighbor.info$n.indx

  # number of reference locations
  k <- length(n.indx)
  idx0 <- which(is.na(n.indx))

  # sort coords by ord
  coords <- coords[ord,]

  # create B and F
  if (base_cov == "exponential") {
    BFlist <- foreach (i = 2:k) %dopar% {
      if (i %in% idx0) {
        B <- 0
        Fmat <- sig_sq
      } else {
        coords_tmp <- rbind(coords[i,], coords[n.indx[[i]],])
        Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
        if (tol > 0) {
          CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1] + diag(tol, nrow(Cor)-1))
        } else {
          CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1])
        }
        B <- CinvC
        Fmat <- sig_sq*(Cor[1,1]-CinvC%*%Cor[-1,1])
      }
      return(list(B = B, Fmat = Fmat))
    }
  } else {
    BFlist <- foreach (i = 2:k) %dopar% {
      if (i %in% idx0) {
        B <- 0
        Fmat <- sig_sq
      } else {
        coords_tmp <- rbind(coords[i,], coords[n.indx[[i]],])
        Cor <- Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                          phi = phi, nu = nu)
        if (tol > 0) {
          CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1] + diag(tol, nrow(Cor)-1))
        } else {
          CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1])
        }
        B <- CinvC
        Fmat <- sig_sq*(Cor[1,1]-CinvC%*%Cor[-1,1])
      }
      return(list(B = B, Fmat = Fmat))
    }
  }

  Felem <- c(sig_sq, sapply(BFlist, function(x) as.numeric(x$Fmat)))
  Fmat <- diag(Felem)
  B <- matrix(0, k, k)
  for (i in 2:k) {
    if (!i %in% idx0) {
      B[i, n.indx[[i]]] <- BFlist[[i-1]]$B
    }
  }

  invImB <- solve(diag(1, k) - B)
  Ctilde <- invImB %*% Fmat %*% t(invImB)
  return(Ctilde[order(ord), order(ord)]) # put back to original ordering
}

# Neighbor GP covariance single (location vs. location)
# old name : NNGPcov (changed 20220214)
NGPcov_s <- function(v1, v2,
                     coords, neighbor.info, sig_sq, phi, nu = NULL,
                     base_cov = "exponential", Ctilde,
                     coords.0, nn.indx.0_ord) {
  # v1, v2: locations to compute covariance
  # coords: coords for reference locations
  # neighbor.info: list containing ord and n.indx
  # sig_sq: spatial variance
  # phi: decay parameter
  # nu: smoothness parameter, applicable only for base_cov = "matern"
  # base_cov = "exponential" or "matern"
  # Ctilde: result of create_Ctilde
  # coords.0: coords for non-reference locations
  # nn.indx.0_ord: nearest neighbor index (ordered by neighbor.info$ord) for coords.0, either a list or a matrix

  list <- is.list(nn.indx.0_ord)
  ord <- neighbor.info$ord
  n.indx <- neighbor.info$n.indx

  # number of reference locations
  k <- length(n.indx)

  # sort coords by ord
  coords <- coords[ord,]
  Ctilde <- Ctilde[ord, ord]

  coords_name <- apply(coords, 1, function(x) paste(x[1], x[2], sep = ","))
  coords0_name <- apply(coords.0, 1, function(x) paste(x[1], x[2], sep = ","))
  v1_name <- paste(v1[1], v1[2], sep = ",")
  v2_name <- paste(v2[1], v2[2], sep = ",")
  si <- which(coords_name == v1_name)
  ui <- which(coords0_name == v1_name)
  sj <- which(coords_name == v2_name)
  uj <- which(coords0_name == v2_name)

  if (length(si) != 0) {
    if (length(sj) != 0) {
      return(as.numeric(Ctilde[si, sj]))
    } else {
      if (list) {
        N_v2 <- nn.indx.0_ord[[uj]]
      } else {
        N_v2 <- nn.indx.0_ord[uj,]
      }

      if (length(N_v2) == 0) {
        return(0)
      } else {
        coords_tmp <- rbind(v2, coords[N_v2,])

        if (base_cov == "exponential") {
          Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
        } else {
          Cor <- Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                            phi = phi, nu = nu)
        }

        B <- Cor[1,-1]%*%solve(Cor[-1,-1])
        return(as.numeric(Ctilde[si, N_v2] %*% t(B)))
      }
    }
  } else {
    if (length(sj) != 0) {
      if (list) {
        N_v1 <- nn.indx.0_ord[[ui]]
      } else {
        N_v1 <- nn.indx.0_ord[ui,]
      }

      if (length(N_v1) == 0) {
        return(0)
      } else {
        coords_tmp <- rbind(v1, coords[N_v1,])

        if (base_cov == "exponential") {
          Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
        } else {
          Cor <- Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                            phi = phi, nu = nu)
        }

        B <- Cor[1,-1]%*%solve(Cor[-1,-1])
        return(as.numeric(B %*% Ctilde[N_v1, sj]))
      }
    } else {
      if (list) {
        N_v1 <- nn.indx.0_ord[[ui]]
        N_v2 <- nn.indx.0_ord[[uj]]
      } else {
        N_v1 <- nn.indx.0_ord[ui,]
        N_v2 <- nn.indx.0_ord[uj,]
      }

      if (length(N_v1) == 0 | length(N_v2) == 0) {
        return(0)
      } else {
        coords_tmp1 <- rbind(v1, coords[N_v1,])

        if (base_cov == "exponential") {
          Cor1 <- exp(-phi*as.matrix(dist(coords_tmp1)))
        } else {
          Cor1 <- Cov_matern(as.matrix(dist(coords_tmp1)), sigmasq = 1,
                             phi = phi, nu = nu)
        }

        CinvC <- Cor1[1,-1]%*%solve(Cor1[-1,-1])
        B1 <- CinvC
        F1 <- sig_sq*(Cor1[1,1]-CinvC%*%Cor1[-1,1])
        coords_tmp2 <- rbind(v2, coords[N_v2,])

        if (base_cov == "exponential") {
          Cor2 <- exp(-phi*as.matrix(dist(coords_tmp2)))
        } else {
          Cor2 <- Cov_matern(as.matrix(dist(coords_tmp2)), sigmasq = 1,
                             phi = phi, nu = nu)
        }

        B2 <- Cor2[1,-1]%*%solve(Cor2[-1,-1])
        return(as.numeric(B1 %*% Ctilde[N_v1, N_v2] %*% t(B2) +
                            (v1_name == v2_name)*F1))
      }
    }
  }
}

# Neighbor GP covariance multiple (location vs. locations)
# old name : NNGPcov2 (changed 20220214)
NGPcov_m <- function(v1, v2_mat,
                     coords, neighbor.info, sig_sq, phi, nu = NULL,
                     base_cov = "exponential", Ctilde,
                     coords.0, nn.indx.0_ord) {
  # v1: location to compute covariance
  # v2_mat: matrix of locations to compute pairwise covariance with v1
  # coords: coords for reference locations
  # neighbor.info: list containing ord and n.indx
  # sig_sq: spatial variance
  # phi: decay parameter
  # nu: smoothness parameter, applicable only for base_cov = "matern"
  # base_cov = "exponential" or "matern"
  # Ctilde: result of create_Ctilde
  # coords.0: coords for non-reference locations
  # nn.indx.0_ord: nearest neighbor index (ordered by neighbor.info$ord) for coords.0, either a list or a matrix

  list <- is.list(nn.indx.0_ord)
  ord <- neighbor.info$ord
  n.indx <- neighbor.info$n.indx

  # number of reference locations
  k <- length(n.indx)

  # sort coords by ord
  coords <- coords[ord,]
  Ctilde <- Ctilde[ord, ord]

  coords_name <- apply(coords, 1, function(x) paste(x[1], x[2], sep = ","))
  coords0_name <- apply(coords.0, 1, function(x) paste(x[1], x[2], sep = ","))
  v1_name <- paste(v1[1], v1[2], sep = ",")
  v2_name <- apply(v2_mat, 1, function(x) paste(x[1], x[2], sep = ","))

  si <- which(coords_name == v1_name)
  ui <- which(coords0_name == v1_name)

  res <- rep(0, nrow(v2_mat))
  for (k in 1:nrow(v2_mat)) {
    sj <- which(coords_name == v2_name[k])
    uj <- which(coords0_name == v2_name[k])

    if (length(si) != 0) {
      if (length(sj) != 0) {
        res[k] <- as.numeric(Ctilde[si, sj])
      } else {
        if (list) {
          N_v2 <- nn.indx.0_ord[[uj]]
        } else {
          N_v2 <- nn.indx.0_ord[uj,]
        }

        if (length(N_v2) == 0) {
          res[k] <- 0
        } else {
          coords_tmp <- rbind(v2_mat[k,], coords[N_v2,])

          if (base_cov == "exponential") {
            Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
          } else {
            Cor <- Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                              phi = phi, nu = nu)
          }

          B <- Cor[1,-1]%*%solve(Cor[-1,-1])
          res[k] <- as.numeric(Ctilde[si, N_v2] %*% t(B))
        }
      }
    } else {
      if (length(sj) != 0) {
        if (list) {
          N_v1 <- nn.indx.0_ord[[ui]]
        } else {
          N_v1 <- nn.indx.0_ord[ui,]
        }

        if (length(N_v1) == 0) {
          res[k] <- 0
        } else {
          coords_tmp <- rbind(v1, coords[N_v1,])

          if (base_cov == "exponential") {
            Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
          } else {
            Cor <- Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                              phi = phi, nu = nu)
          }

          B <- Cor[1,-1]%*%solve(Cor[-1,-1])
          res[k] <- as.numeric(B %*% Ctilde[N_v1, sj])
        }
      } else {
        if (list) {
          N_v1 <- nn.indx.0_ord[[ui]]
          N_v2 <- nn.indx.0_ord[[uj]]
        } else {
          N_v1 <- nn.indx.0_ord[ui,]
          N_v2 <- nn.indx.0_ord[uj,]
        }

        if (length(N_v1) == 0 | length(N_v2) == 0) {
          res[k] <- 0
        } else {
          coords_tmp1 <- rbind(v1, coords[N_v1,])

          if (base_cov == "exponential") {
            Cor1 <- exp(-phi*as.matrix(dist(coords_tmp1)))
          } else {
            Cor1 <- Cov_matern(as.matrix(dist(coords_tmp1)), sigmasq = 1,
                               phi = phi, nu = nu)
          }

          CinvC <- Cor1[1,-1]%*%solve(Cor1[-1,-1])
          B1 <- CinvC
          F1 <- sig_sq*(Cor1[1,1]-CinvC%*%Cor1[-1,1])
          coords_tmp2 <- rbind(v2_mat[k,], coords[N_v2,])

          if (base_cov == "exponential") {
            Cor2 <- exp(-phi*as.matrix(dist(coords_tmp2)))
          } else {
            Cor2 <- Cov_matern(as.matrix(dist(coords_tmp2)), sigmasq = 1,
                               phi = phi, nu = nu)
          }

          B2 <- Cor2[1,-1]%*%solve(Cor2[-1,-1])
          res[k] <- as.numeric(B1 %*% Ctilde[N_v1, N_v2] %*% t(B2) +
                                 (v1_name == v2_name[k])*F1)
        }
      }
    }
  }
  return(res)
}


