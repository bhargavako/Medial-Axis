#' Medial Axis using Delaunay Triangulation
#'
#' This function computes the medial axis for the input polygon
#' @param input
#' @keywords polygon
#' @export
#' @examples
#' medial()
medial <- function(Shape) {
  require(deldir)
  require(plyr)
  require(rgdal)
  require(iterators)
  require(doParallel)
  require(parallel)
  require(foreach)
  require(sp)
  #Detecting the number of cores available
  num_cores <- detectCores() - 1
  clust <- makeCluster(num_cores, type = "PSOCK")
  registerDoParallel(clust)

  idpoly <- lapply(Shape@polygons , slot , "ID")
  holpoly <- sapply(Shape@polygons[[1]]@Polygons , slot , "hole")

  results <- foreach(
    j = 1:iter(idpoly)$length,
    .combine = 'list',
    .multicombine = TRUE,
    .packages = c('deldir', 'plyr', 'sp', 'iterators')
  ) %dopar% {
    orderpoly <- Shape@polygons[[j]]@plotOrder
    if (length(orderpoly) == 1) {
      coord <- Shape@polygons[[j]]@Polygons[[1]]@coords
      coord <- coord[-1,]

      voro.tri <- deldir(coord[, 1], coord[, 2])
      voro.tile <- tile.list(voro.tri)

      voro.vert <-
        plyr::ldply(voro.tile, function(coord)
          cbind(coord$x, coord$y))

      ## Remove dublicated vertices
      voro.vert <- unique(voro.vert, MARGIN = 1)

      ## Select vertices inside the polygon
      idx <-
        point.in.polygon(voro.vert[, 1], voro.vert[, 2], coord[, 1], coord[, 2])
      voro.vert <- voro.vert[idx == 1,]

      ## "Minimization" part of optimization.
      ## Distance to closest shape-point
      min.dist <-
        apply(voro.vert, 1, function(x)
          min(spDistsN1(coord, x)))
      polytile <-
        list("tiles" = voro.tile,
             "ctroid" = voro.vert[which.max(min.dist), ],
             "coord" = coord)

    } else {
      #Used 'for' loop to maintain sequence
      for (i in 1:iter(orderpoly)$length) {
        holechk <- Shape@polygons[[1]]@Polygons[[i]]@hole
        coord <- Shape@polygons[[1]]@Polygons[[i]]@coords
        coord <- coord[-1,]
        if (!holechk) {
          ocoord <- coord
        } else {
          hcoord <- coord
          lpoly <- list("hole" = coord)
        }
        if (!exists("mcoord")) {
          mcoord <- NULL
        }
        mcoord <- rbind(mcoord, coord)
      }

      voro.tri <- deldir(mcoord[, 1], mcoord[, 2])
      voro.tile <- tile.list(voro.tri)

      voro.vert <-
        plyr::ldply(voro.tile, function(ocoord)
          cbind(ocoord$x, ocoord$y))

      ## Remove dublicated vertices
      voro.vert <- unique(voro.vert, MARGIN = 1)

      ## Select vertices inside the polygon
      pip <- function(voro.vert, xim) {
        point.in.polygon(voro.vert[, 1], voro.vert[, 2], xim[, 1], xim[, 2])
      }

      idx <- pip(voro.vert, ocoord)
      idx2 <- pip(voro.vert, lpoly$hole)

      idx3 <- idx - idx2

      voro.vert <- voro.vert[idx3 == 1,]

      ## "Minimization" part of optimization.
      ## Distance to closest shape-point
      min.dist <-
        apply(voro.vert, 1, function(x)
          min(spDistsN1(ocoord, x)))
      polytile <-
        list("tiles" = voro.tile,
             "ctroid" = voro.vert[which.max(min.dist), ],
             "coord" = mcoord)
    }
  }

  stopCluster(clust)

  plot(Shape, axes = T, asp = 1)
  if (iter(results$tiles)$length != 0) {
    for (res in 1:iter(results)$length) {
      plot(results$tiles, close = FALSE, add = T)
      if (TRUE %in% holpoly) {
        plot(SpatialPolygons(lapply(Shape@polygons,
                                    function(x) {
                                      x@Polygons = x@Polygons[2]
                                      x
                                    })),
             col = "white",
             add = T)
      }
      points(results$coord,

             col = "steelblue4",
             cex = 1)
      points(results$ctroid,
             col = "orangered2",
             pch = 20,
             cex = 2)
      box()
    }

  } else {
    for (res in 1:iter(results)$length) {
      plot(results[[res]]$tiles, close = FALSE, add = T)
      points(results[[res]]$coord,
             col = "steelblue4",
             cex = 1)
      points(results[[res]]$ctroid,
             col = "orangered2",
             pch = 20,
             cex = 2)
      box()
    }
  }
  return(list("ctroid" = results$ctroid,
                 "coord" = results$coord))
}
