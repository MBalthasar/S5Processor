# S5Processor
A R package for processing Sentinel-5P L2 data.

## Installation
To install the current version, use `devtools`.

```R
devtools::install_github("MBalthasar/S5Processor")
library(S5Processor)
```

## Available Functions
The following function is currently available and tested on Windows 10.

* `S5P_process()` This function automatically converts the Sentinel-5P L2 NC-data files into TIFF files. The user can define an additional aoi, the preferred resolution as well as the resampling method.

## Example

```R
# Load required packages
library(ncdf4)
library(ggplot2)
library(dismo)
library(maptools)
library(raster)
library(geosphere)
library(rgdal)
library(rgeos)
library(sp)


# Load sample NetCDF files
x1 <- system.file(package = "S5Processor", "extdata", "S5P_NRTI_L2__NO2_1.nc")
x2 <- system.file(package = "S5Processor", "extdata", "S5P_NRTI_L2__NO2_2.nc")

# Load sample shapefile including the borders of Vietnam
vnm_shp <- raster::shapefile(system.file(package = "S5Processor",
                             "extdata", "vietnam_borders.shp"))

# Create vector from both NetCDF files.
my_files <- c(x1,x2)

# Most basic case: provide path to single NetCDF file.
# The user will be promted to choose from a list of product names.
S5P_1 <- S5P_process(input = my_files[1])

# Execute function on several files at once; additionally provide product number.
S5P_2 <- S5P_process(input = my_files, product = 39)

# As above but this time also mask data to aoi and define resolution.
S5P_3 <- S5P_process(input = my_files, my_res = 10000,
                     product = 39, my_aoi = vnm_shp,
                     extent_only = FALSE)

# This time also apply scale factor to convert units to molecules per cm2.
S5P_4 <- S5P_process(input = my_files, my_res = 10000,
                     product = 39, my_aoi = vnm_shp,
                     extent_only = FALSE,
                     apply_scale_factor = T)
```

## References
This program was developed at Remote Sensing Solutions RSS GmbH, Dingolfinger Str. 9, 81673 Munich, Germany, https://rssgmbh.de/
