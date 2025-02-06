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

    experiment <- "../f19b_biogemtest.SPIN"
    var <- "ocn_O2"
    demfname <- "../Map01_PALEOMAP_1deg_Holocene_0Ma.nc"
    demzname <- "z"

    library(RNetCDF)
    library(dplyr)
    library(reshape2)
    library(PaleoClimR)

    # Get dataframe of cGENIE grid and tracer values
    df_cGENIE <- cGENIE.data.3D(var, experiment,
                                "default", "biogem")
    
    # Extract DEM lat/lon/depth
    nc <- open.nc(demfname)
    lat <- var.get.nc(nc, "lat")
    lon <- var.get.nc(nc, "lon")
    depth <- var.get.nc(nc, "z")
    # Set land to NA
    depth[which(depth>=0)] <- NA
    # Convert topography to depths
    depth <- -depth
    # Clip depths to range of cGENIE mid-layer depths
    depth[depth < min(df_cGENIE$depth)] <- min(df_cGENIE$depth)
    depth[depth > max(df_cGENIE$depth)] <- max(df_cGENIE$depth)

    # Create a dataframe for the DEM data
    df_DEM <- as.data.frame(cbind(
        rep(lon, times = length(lat), each = 1),
        rep(lat, times = 1, each = length(lon)),
        as.data.frame(melt(depth))$value))
    names(df_DEM) <- c("lon", "lat", "depth")

    # Adjust DEM grid for longitude projection (to 0-360 degrees if needed)
    if (mean(between(df_DEM$lon, -180, 180)) < 1) {
        df_DEM$lon[df_DEM$lon <= -180] <- df_DEM$lon[df_DEM$lon <= -180] + 360                            
    }

    # Generate LOESS models for each unique lat/lon
    loess_models <- df_cGENIE %>%
    group_by(lat.mid, lon.mid) %>%
    # Exclude cells where var is NA over all depths (i.e. land cells)
    filter(!all(is.nan(var))) %>%
    summarize(
        model = list(loess(var ~ depth, 
            data = cur_data(), 
            na.action = na.exclude,
            span = 0.5,
            degree = 2,
            )
        )
    )

    # Estimate var at each DEM point using closest LOESS model in lat/lon space
    df_DEM <- df_DEM %>%
    rowwise() %>%
    mutate(
        var = {
            predict(loess_models$model[[
                which.min(
                    sqrt(
                        (loess_models$lat.mid - lat)^2 
                        + (loess_models$lon.mid - lon)^2
                    )
                )
            ]], 
            depth)
        }
    )
      
