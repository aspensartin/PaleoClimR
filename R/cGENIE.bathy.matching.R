#' cGENIE Bathymetry Matching
#'
#' This function matches a cGENIE 3D tracer field to bathymetry from a (palaeo)
#' digital elevation model (DEM). It interpolates values for the chosen tracer
#' at the depth of each DEM point from the nearest cGENIE grid cell using a
#' locally-smoothed depth profile of the tracer in that cell. The function
#' returns a dataframe of depths and interpolated values.
#'
#' @param var A string specifying the variable name to extract from the NetCDF
#'              file. Examples include "ocn_temp", "ocn_O2", etc.
#' @param experiment A string specifying the path to the folder or file prefix
#'                   where the experiment's data is stored.
#' @param demfname A string specifying the path to the DEM NetCDF file.
#' @param demzname A string specifying the name of the depth variable in the
#'                 DEM file.
#'
#' @return A dataframe (`df.sum`) containing:
#'   - `lon`: Longitude of DEM points.
#'   - `lat`: Latitude of DEM points.
#'   - `depth`: Depth at DEM points.
#'   - `var`: Interpolated tracer values (e.g., [O2])
#'
#' @import RNetCDF, dplyr, reshape2
#' @export
#'
cGENIE.bathy.matching <- function(var, experiment, demfname, demzname) {

    library(RNetCDF)
    library(dplyr)
    library(reshape2)

    # Get dataframe of cGENIE grid and tracer values
    df.cGENIE <- cGENIE.data.3D(var, experiment, "default", "biogem")
    
    # Extract DEM lat/lon/depth
    nc <- open.nc(demfname)
    lat <- var.get.nc(nc, "lat")
    lon <- var.get.nc(nc, "lon")
    depth <- var.get.nc(nc, demzname)

    # Adjust DEM grid for longitude projection (to 0-360 degrees if needed)
    if (mean(between(dem.lon, -180, 180)) < 1) {
        dem.lon[dem.lon <= -180] <- dem.lon[dem.lon <= -180] + 360                            
    }

    # Create a dataframe for the DEM data
    df.DEM <- as.data.frame(cbind(
        rep(lon, times = length(lat), each = 1),
        rep(lat, times = 1, each = length(lon)),
        as.data.frame(melt(depth))$value))
    names(df.DEM) <- c("lon", "lat", "depth")




                                }