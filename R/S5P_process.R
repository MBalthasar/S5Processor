#' Sentinel-5P L2 Processing Tool
#'
#' This function automatically converts the Sentinel-5P L2 NC-data files into TIFF files.
#' The user can define a additional aoi, the preferred resolution as well as the resampling method.
#'
#' @param input A vector or single character string containing the path(s) to the NetCDF file(s).
#'
#' @param product (optional) Define a number for the product that should be transformed into
#'                a raster file. If not defined, the user will be prompted to choose from a list when
#'                executing the function.
#'
#' @param my_res (optional) Define geometric resolution in m for the final raster file.
#'               The actual resolution will likely differ slightly from the users preset, since
#'               the extent of the NetCDF data or provided aoi/reference data will be used and turned
#'               into a raster file. The dimensions of the extent will be devided by the resolution and
#'               the resulting decimal number will be rounded to an integer to provide the number of
#'               rows and columns for the resulting raster file -> Default is 20000m.
#'
#' @param my_aoi (optional) Shapefile which will be used to crop the data.
#'
#' @param extent_only (optional) If TRUE it will only use the extent of the aoi or reference raster to
#'                    crop the data, if FALSE (default) the NetCDF data will be masked to the
#'                    aoi/ref_raster.
#'
#' @param ref_raster (optional) Raster file to which the S5P data will be cropped, projected,
#'                   including resolution.
#'
#' @param interpol_method (optional) Either 'ngb' (nearest neighbor) or 'bilinear'
#'                        (bilinear interpolation; the default value).
#'
#' @param apply_scale_factor (optional) If TRUE it will convert the pixel values to molecules per cm2
#'                           if a multiplication factor is available for the defined product. If there
#'                           is no multiplication factor available or set to FALSE (default) the original
#'                           unit of the product will be used.
#'
#' @return A Tiff file.
#'
#' @references This program was written at Remote Sensing Solutions RSS GmbH,
#'             Dingolfinger Str. 9, 81673 Munich, Germany \url{https://rssgmbh.de/}.
#'
#' @examples
#' # Load required packages
#' library(ncdf4)
#' library(ggplot2)
#' library(dismo)
#' library(maptools)
#' library(raster)
#' library(geosphere)
#' library(rgdal)
#' library(rgeos)
#' library(sp)
#' library(S5Processor)
#'
#'
#' # Load sample NetCDF files
#' x1 <- system.file(package = "S5Processor", "extdata", "S5P_NRTI_L2__NO2_1.nc")
#' x2 <- system.file(package = "S5Processor", "extdata", "S5P_NRTI_L2__NO2_2.nc")
#'
#' # Load sample shapefile including the borders of vietnam
#' vnm_shp <- raster::shapefile(system.file(package = "S5Processor",
#'                              "extdata", "vietnam_borders.shp"))
#'
#' # Create vector from both NetCDF files.
#' my_files <- c(x1,x2)
#'
#' # Most basic case: provide path to single NetCDF file.
#' # The user will be promted to choose from a list of product names.
#' S5P_1 <- S5P_process(input = my_files[1])
#'
#' # Execute function on several files at once; additionally provide product number.
#' S5P_2 <- S5P_process(input = my_files, product = 39)
#'
#' # As above but this time also mask data to aoi and define resolution.
#' S5P_3 <- S5P_process(input = my_files, my_res = 10000,
#'                      product = 39, my_aoi = vnm_shp,
#'                      extent_only = FALSE)
#'
#' # This time also apply scale factor to convert units to molecules per cm2.
#' S5P_4 <- S5P_process(input = my_files, my_res = 10000,
#'                      product = 39, my_aoi = vnm_shp,
#'                      extent_only = FALSE,
#'                      apply_scale_factor = T)
#'
#' @export
S5P_process <- function(input, product, my_res, my_aoi, extent_only,
                        ref_raster, interpol_method, apply_scale_factor){
  # Check if variables are correctly defined
  if (missing(my_aoi)){
    print(paste0("Variables checked"))
  } else {
    if (missing(ref_raster)){
      print(paste0("Variables checked"))
    } else {
      stop("You can only define EITHER an aoi OR a reference raster")
    }
  }

  #############################
  ### get product from user ###
  #############################
  print(paste0("Receiving product information"))
  # open first file
  first_nc <- ncdf4::nc_open(input[1])
  # get product number
  if (missing(product)){
    # display all the available variables
    print(attributes(first_nc$var)$names)
    # Prompt user to specify which product should be selected
    user_input <- readline(prompt="Enter an integer for product selection: ")
    user_input <- as.integer(user_input)
  } else {
    user_input <- product
  }

  ########################################################################
  ### Start for loop creating a data frame with coordinates and values ###
  ########################################################################
  # Measure time
  start.time <- Sys.time()
  for (i in 1:length(input)){
    print(paste0("Getting values from file ", i, " of ", length(input)))
    # open NetCDF file
    nc <- ncdf4::nc_open(input[i])
    # save the multiplication factor, fill value, unit and name of product
    mfactor = ncdf4::ncatt_get(nc, attributes(nc$var)$names[user_input],
                        "multiplication_factor_to_convert_to_molecules_percm2")
    fillvalue = ncdf4::ncatt_get(nc, attributes(nc$var)$names[user_input],
                          "_FillValue")
    my_unit = ncdf4::ncatt_get(nc, attributes(nc$var)$names[user_input],
                        "units")
    my_product_name = ncdf4::ncatt_get(nc, attributes(nc$var)$names[user_input],
                                "long_name")
    # Check if multiplication factor and fill value are present
    if (mfactor$hasatt == TRUE){
      mfactor <- mfactor$value
    } else {
      mfactor <- 1
    }
    if (fillvalue$hasatt == TRUE){
      fillvalue <- fillvalue$value
    } else {
      fillvalue <- 9.96921e+36
    }
    # read the data values, along with the latitude and longitude information
    vals <- ncdf4::ncvar_get(nc, attributes(nc$var)$names[user_input])
    lat <- ncdf4::ncvar_get(nc, "PRODUCT/latitude")
    lon <- ncdf4::ncvar_get(nc, "PRODUCT/longitude")
    # set fill values to NA
    vals[vals == fillvalue] <- NA
    # Check if scale factor should be applied
    if (missing(apply_scale_factor)){
      vals <- vals
    } else if (apply_scale_factor == FALSE){
      vals <- vals
    } else if (apply_scale_factor == TRUE){
      vals <- vals*mfactor
    }
    # convert information into data frame
    vals_df = NULL
    # apply multiplication factor for unit conversion
    # vals <- vals*mfactor
    # combine lon lat and values
    vals_df <- rbind(vals_df, data.frame(lat=as.vector(lat),
                                         lon=as.vector(lon),
                                         vals=as.vector(vals)))
    if (i == 1){
      final_df <- vals_df
    } else {
      final_df <- rbind(final_df, vals_df)
    }
  }

  #############################
  ### Covert data to points ###
  #############################
  print(paste0("Convert data frame to spatial points"))
  # Convert df into spatialpointsdataframe #
  pts <- final_df
  sp::coordinates(pts) <- ~ lon + lat
  # Give projection to points
  my_projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  sp::proj4string(pts) <- sp::CRS(my_projection)

  ##############################################
  ### Crop points to aoi or reference raster ###
  ##############################################
  # check if aoi is missing
  if (missing(my_aoi)){
    # check if reference raster is missing
    if (missing(ref_raster)){
      pts <- pts
    } else {
      # if reference raster is present use it to crop the data but keep lon lat projection for now
      ref_raster_layer <- ref_raster[[1]]
      print(paste0("Crop points to reference raster"))
      p <- methods::as(raster::extent(ref_raster_layer), 'SpatialPolygons')
      sp::proj4string(p) <- sp::CRS(as.character(raster::crs(ref_raster_layer)))
      p <- sp::spTransform(p, CRS=my_projection)
      pts <- raster::crop(pts, p)
    }
    # if aoi is defined, use it instead to crop the data
  } else {
    # check projection of aoi
    crs_test <- raster::compareCRS(pts, my_aoi)
    if (crs_test == FALSE){
      my_aoi <- sp::spTransform(my_aoi, CRS=as.character(raster::crs(pts)))
    } else {
      my_aoi <- my_aoi
    }
    print(paste0("Crop points to aoi"))
    # crop pts to extent of aoi
    p <- methods::as(raster::extent(my_aoi), 'SpatialPolygons')
    sp::proj4string(p) <- sp::CRS(my_projection)
    pts <- raster::crop(pts, p)
  }

  #####################################
  ### Convert points to raster data ###
  #####################################
  print(paste0("Calculate number of rows and columns for raster file"))
  # calculate ymin/ymax and xmin/xmax distances of extent
  # for lat its easy since it doesn't matter from which longitutde we measure the vertical distance
  extent_distance_vertical <- geosphere::distm(c(raster::extent(pts)[1], raster::extent(pts)[3]),
                                    c(raster::extent(pts)[1], raster::extent(pts)[4]),
                                    fun = geosphere::distHaversine)
  # BUT: longitude distance depends on latitude
  # create new lat coordinate in the vertical middle for the western and easter boundary of the extent
  vertical_mid_distance <- (raster::extent(pts)[4] - raster::extent(pts)[3])/2
  # add vertical middle to minimum latitude <- get mid latitude of extent
  lat_mid <- raster::extent(pts)[3] + vertical_mid_distance
  # get total longitudinal distance
  horizontal_distance <- raster::extent(pts)[2] - raster::extent(pts)[1]
  # Check if horizontal distance is > 180; if yes we need to manually calculate the distance since it always takes the
  # shortest route to two points
  if (horizontal_distance > 180){
    # get distance for one horizontal degree in given latitude
    one_degree_horizontal_distance <- geosphere::distm(c(1, lat_mid),
                                            c(2, lat_mid),
                                            fun = geosphere::distHaversine)
    # multiply with total longitudal degrees
    extent_distance_horizontal <- one_degree_horizontal_distance * horizontal_distance
  } else {
    # calcualte longitudinal distance
    extent_distance_horizontal <- geosphere::distm(c(raster::extent(pts)[1], lat_mid),
                                        c(raster::extent(pts)[2], lat_mid),
                                        fun = geosphere::distHaversine)
  }
  # Define resolution
  if (missing(my_res)){
    # use default resolution of 20km
    my_res <- 20000
  } else {
    my_res <- my_res
  }
  # calculate number of rows and columns based on resolution
  ncol_rast <- as.integer(extent_distance_horizontal/my_res)
  nrow_rast <- as.integer(extent_distance_vertical/my_res)
  print(paste0("Create raster file from points"))
  # create raster
  rast <- raster::raster(nrows=nrow_rast, ncols=ncol_rast, crs=as.character(raster::crs(pts)),
                 ext= raster::extent(pts), vals=NULL)
  # Rasterize the points
  final <- raster::rasterize(pts, rast, pts$vals, fun=mean)

  ############################################################################
  ### Change features of final raster to match reference raster if present ###
  ############################################################################
  if (missing(ref_raster)){
    final <- final
  } else {
    print(paste0("Project the values to match reference raster"))
    # Check which interpolation method should be used; default is bilinear
    if (missing(interpol_method)){
      final <- raster::projectRaster(from = final, to = ref_raster_layer,
                             crs=as.character(raster::crs(ref_raster_layer)), method = 'bilinear')
    } else {
      final <- raster::projectRaster(from = final, to = ref_raster_layer,
                             crs=as.character(raster::crs(ref_raster_layer)), method = interpol_method)
    }
  }
  # Check if raster file should be masked
  if (missing(my_aoi)){
    if (missing(ref_raster)){
      final <- final
    } else {
      # Check if masking it only to the extent
      if (missing(extent_only)){
        print(paste0("Mask final raster to reference raster"))
        final <- raster::mask(final, ref_raster_layer)
      } else if (extent_only == FALSE){
        print(paste0("Mask final raster to reference raster"))
        final <- raster::mask(final, ref_raster_layer)
      } else if (extent_only == TRUE){
        final <- final
      }
    }
  } else {
    print(paste0("Mask final raster to aoi"))
    if (missing(extent_only)){
      my_aoi <- sp::spTransform(my_aoi, CRS=raster::crs(final))
      final <- raster::mask(final, my_aoi)
    } else if (extent_only == FALSE){
      my_aoi <- sp::spTransform(my_aoi, CRS=raster::crs(final))
      final <- raster::mask(final, my_aoi)
    } else if (extent_only == TRUE){
      final <- final
    }
  }

  ##############################################################
  ### Stop time measurement and print additional information ###
  ##############################################################
  stop.time <- Sys.time()
  time.taken <- stop.time - start.time
  print(paste0("====================== Done! ======================"))
  print(paste0("Total processing time: ", time.taken, " seconds"))
  if (missing(apply_scale_factor)){
    if (mfactor != 1){
      print(paste0(my_product_name$value, " in the unit ", my_unit$value))
      print(paste0("Multiplication factor to convert to molecules per cm2: ", mfactor))
    } else {
      print(paste0(my_product_name$value, " in the unit ", my_unit$value))
    }
  } else if (apply_scale_factor == FALSE){
    if (mfactor != 1){
      print(paste0(my_product_name$value, " in the unit ", my_unit$value))
      print(paste0("Multiplication factor to convert to molecules per cm2: ", mfactor))
    } else {
      print(paste0(my_product_name$value, " in the unit ", my_unit$value))
    }
  } else if (apply_scale_factor == TRUE){
    if (mfactor != 1){
      print(paste0(my_product_name$value))
      print(paste0("Data has been converted from ", my_unit$value,
                   " to molecules per cm2 by multiplying with the factor ", mfactor))
    } else {
      print(paste0(my_product_name$value))
      print(paste0("No multiplicatin factor has been applied to convert to molecules per cm2, since there is non."))
    }
  }
  raster::plot(final)
  return(final)
}

