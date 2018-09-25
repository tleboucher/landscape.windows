#' Landscape windows: random method
#'
#' This function creates landscape windows into an area defined by a shapefile (shape format) using the 'random method'. n windows are placed randomly in the defined area. With default values, to be validated, a window must contain at least 25 sites distributed over at least 80 percents of the window’s area (at least 1 site in 8 cells of a 3 x 3 grid).
#' @param data Informations about sites (spatial coordinates). 'data' must have 3 columns: site ID, x and y.
#' @param shape.path Absolute path to the shapefile (study area).
#' @param proj Projection system. 'proj' must be expressed in EPSG code or in CRS.
#' @param windows Informations about the desired number of windows. 'windows' must contain 2 values: the desired number of windows per scale and the number of attempts (should be a multiplicative factor to multiply the desired number of windows).
#' @param scales Spatial scales to test. 'scales' must contain 3 values: the minimum windows' size, the maximum windows' size and the gap/break between two scales (in kilometers).
#' @param min.no.sites Minimum number of sites in each window (default = 25).
#' @param area Minimum percentage for the sites distributed over the window’s area to ensure equal sampling effort across categories (default = 80).
#' @param nrowcol Area is calculated with a grid of n rows/cols. 'nrowcol' determines n value (default = 3).
#' @param plot TRUE or FALSE. Creation of a map of the study area with the centers of the validated windows (default = FALSE).
#' @return Output data is a list of 3 objects: a table containing all the validated windows, for each scale, with a list of sites present in each window; a SpatialPointsDataFrame containing the data (sites); and a SpatialPolygon of the study area.
#' @export

landscape.windows.random <- function(data, shape.path, proj, windows, scales, min.no.sites = 25, area = 80, nrowcol = 3, plot = FALSE){

  if(base::length(scales)!=3)stop("scales must contain 3 values: the minimum windows' size, the maximum windows' size and the gap/break between two scales (in kilometers)")
  if(base::length(windows)!=2)stop("windows must contain 2 values: the desired number of windows per scale and the number of attempts (should be a multiplicative factor to multiply the desired number of windows)")
  nrowcol <- base::ceiling(x = nrowcol)
  area <- base::ceiling(x = area)
  method = "random"

  ##############################
  presence.absence.raster <- function(mask.raster, points.data, raster.label = "") {
    # http://www.r-bloggers.com/creating-a-presence-absence-raster-from-point-data/
    # Author: Amy Whitehead
    # May 27, 2013
    mask.raster[!is.na(mask.raster)] <- 0
    pointsRaster <- raster::rasterize(x = points.data, y = mask.raster, field = 1)
    pointsRaster <- raster::merge(x = pointsRaster, y = mask.raster)
    base::names(x = pointsRaster) <- raster.label
    return(pointsRaster)
  }
  ##############################

  oldw <- base::getOption(x = "warn")
  base::options(warn = -1)
  if(base::nchar(x = proj)>10){
    shape <- maptools::readShapeSpatial(fn = base::paste0(shape.path), proj4string = sp::CRS(projargs = proj))
  }else{
    if(rgdal::checkCRSArgs(uprojargs = base::paste0("+init=epsg:",base::as.numeric(x = base::as.character(x = proj))))[[1]]){}else{stop("Undefined EPSG argument.")}
    shape <- maptools::readShapeSpatial(fn = paste0(shape.path), proj4string = sp::CRS(projargs = paste0("+init=epsg:",base::as.numeric(x = base::as.character(x = proj)))))
  }

  shape.save <- shape

  base::options(warn = oldw)

  robinson<-"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs"
  proj.save <- shape@proj4string@projargs

  if(base::dim(x = data)[2]!=3)stop("Data must have 3 columns: site ID, x and y.")
  base::colnames(x = data) <- c("site","x","y")
  base::rownames(x = data) <- data$site
  data <- data[,-1]
  sp::coordinates(obj = data) <- c("x","y")
  sp::proj4string(obj = data) <- sp::CRS(projargs = proj.save)
  data.save <- data
  data <- sp::spTransform(x = data, CRSobj = sp::CRS(projargs = robinson))
  shape <- sp::spTransform(x = shape, CRSobj = sp::CRS(projargs = robinson))

  if(shape@bbox[1,1]>data@bbox[1,1])stop("Data seems to be out of the shape's area or data and shape are not in the same projection system.")
  if(shape@bbox[1,2]<data@bbox[1,2])stop("Data seems to be out of the shape's area or data and shape are not in the same projection system.")
  if(shape@bbox[2,1]>data@bbox[2,1])stop("Data seems to be out of the shape's area or data and shape are not in the same projection system.")
  if(shape@bbox[2,2]<data@bbox[2,2])stop("Data seems to be out of the shape's area or data and shape are not in the same projection system.")

  no.windows <- base::as.numeric(base::as.character(windows[1]))
  no.test <- base::as.numeric(base::as.character(windows[2]))
  aa <- base::as.numeric(base::as.character(scales[1]))
  bb <- base::as.numeric(base::as.character(scales[2]))
  cc <- base::as.numeric(base::as.character(scales[3]))
  size <- base::seq(from = aa, to = bb, by = cc)

  result <- NULL
  for(l in size){

    temp2 <- 0
    temp3 <- 0
    while(temp3<(no.windows*no.test)){

      if((temp3>=base::ceiling(x = no.windows*(0.5*no.test)) & temp2==0)){
        base::cat(base::paste0("\nStop scale ",l," kilometers: there are no windows that match the defined conditions after ",base::ceiling(x = no.windows*(0.5*no.test))," attempts.\n"))
        break()}

      x_win <- stats::runif(n = 1, min = shape@bbox[1,1], max = shape@bbox[1,2])
      y_win <- stats::runif(n = 1, min = shape@bbox[2,1], max = shape@bbox[2,2])

      point <- base::as.data.frame(x = cbind(x_win,y_win))
      sp::coordinates(obj = point) <- c("x_win","y_win")
      sp::proj4string(obj = point) <- sp::CRS(projargs = robinson)
      test <- sp::over(x = point, y = shape)
      if(!base::is.na(test[,1])){

        temp3 <- temp3 + 1

        xmn <- xmx <- ymn <- ymx <- xy <- rast <- NULL
        xmn <- x_win-(l*500)
        xmx <- x_win+(l*500)
        ymn <- y_win-(l*500)
        ymx <- y_win+(l*500)
        xy <- base::matrix(data = 0, nrow = nrowcol, ncol = nrowcol)
        rast <- raster::raster(x = xy)
        raster::extent(x = rast) <- c(xmn,xmx,ymn,ymx)
        raster::projection(x = rast) <- sp::CRS(projargs = robinson)

        pa.raster <- area.test <- NULL
        pa.raster <- presence.absence.raster(mask.raster = rast, points.data = base::as.data.frame(x = data@coords), raster.label = base::paste0("ID_",temp3))
        area.test <- base::round(x = (base::length(x = pa.raster@data@values[pa.raster@data@values==1])/(nrowcol*nrowcol))*100,2)

        if(area.test>=area){

          list.sites <- base::rownames(x = raster::intersect(x = data, y = rast)@coords)

          if(base::length(x = list.sites)>min.no.sites){

            point <- sp::spTransform(x = point, CRSobj = sp::CRS(projargs = proj.save))
            x_win <- point@coords[1,1]
            y_win <- point@coords[1,2]
            base::rm(point)

            temp <- base::cbind(base::paste0(l,temp3),x_win,y_win,l,list.sites)
            temp2 <- temp2+1
            result <- base::rbind(result,temp)

            if(temp2>(no.windows-1)){
              break()}

          }
        }
      }
    }
  }

  if(base::is.null(base::dim(x = result))){stop("There is no window that match the defined conditions.")}else{

    result <- base::as.data.frame(x = result)
    base::colnames(x = result) <- c("no.window","x.window","y.window","scale","site")
    result$x.window <- base::as.numeric(x = base::as.character(x = result$x.window))
    result$y.window <- base::as.numeric(x = base::as.character(x = result$y.window))
    result$scale <- base::as.numeric(x = base::as.character(x = result$scale))
    result$no.window <- base::as.numeric(x = base::as.character(x = result$no.window))
    result$site <- base::as.character(x = result$site)
    final <- base::list(result,shape.save,data.save)
    base::names(x = final) <- c("result","shape","sites")

    base::cat("\n")
    base::cat(base::paste0("==========================================================================================================================================\n"))
    base::cat(base::paste0("|                                                       Landscape windows function                                                       |\n"))
    base::cat(base::paste0("==========================================================================================================================================\n"))
    base::cat(base::paste0("\n\n"))
    base::cat(base::paste0("==>  Minimum number of sites per window = ", min.no.sites,"\n"))
    base::cat(base::paste0("==>  Sites must be distributed over at least ",area,"% of the window’s area.\n"))
    base::cat(base::paste0("==>  The ",method," method has been selected.\n"))

    base::cat(base::paste0("==>  The centers of each windows are randomly placed and the desired number of windows per scale is ",no.windows,".\n"))
    base::cat(base::paste0("==>  Creation of windows is tested in these conditions for windows'size ranging between ",aa," and ",bb," kilometers (every ",cc," kilometers).\n"))
    base::cat(base::paste0("==>  There are windows that match the defined conditions between ", base::min(result$scale)," and ",base::max(result$scale)," kilometers.\n"))
    base::cat(base::paste0("==>  After ",(no.windows*no.test)," attemps at most, there are:\n\n"))
    for(i in base::unique(x = result$scale)[base::order(base::unique(x = result$scale))]){
      base::cat(base::paste0("            ",base::length(x = base::unique(x = result$no.window[result$scale==i]))," windows that match the defined conditions for the windows'size ",i," kilometers.\n"))
    }

    base::cat(base::paste0("\n\n"))
    base::cat(base::paste0("==========================================================================================================================================\n"))

    if(plot==TRUE){
      plot(x = shape.save)
      points(data.save, pch = 21, bg = "grey80", cex = 0.6)
      pal<-base::rep(x = c("#A50026","#D73027","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#D9EF8B","#A6D96A","#66BD63","#1A9850","#006837"), times = 10)
      temp.col <- 0
      for(i in base::unique(x = result$scale)[base::order(base::unique(x = result$scale))]){
        temp.col <- temp.col+1
        points(result$x.window[result$scale==i],result$y.window[result$scale==i], pch = 21, bg = pal[temp.col])
      }
    }
    return(final)
  }
}
