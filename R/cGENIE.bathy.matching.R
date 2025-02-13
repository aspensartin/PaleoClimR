#' cGENIE Bathymetry Matching
#'
#' [function description]
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
    library(interp)
    library(PaleoClimR)
    library(doParallel)

    registerDoParallel(cores = 4)

    # Get dataframe of cGENIE grid and tracer values
    df_cGENIE <- cGENIE.data.3D(var, experiment,
                                "default", "biogem")

    df_cGENIE <- df_cGENIE %>%
        filter(is.finite(var))                            
    
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

    # Create DEM dataframe
    df_DEM <- as.data.frame(cbind(
        rep(lon, times = length(lat), each = 1),
        rep(lat, times = 1, each = length(lon)),
        as.data.frame(melt(depth))$value))
    names(df_DEM) <- c("lon", "lat", "depth")

    # Adjust DEM grid for longitude projection (to 0-360 degrees if needed)
    if (mean(between(df_DEM$lon, -180, 180)) < 1) {
        df_DEM$lon[df_DEM$lon <= -180] <- df_DEM$lon[df_DEM$lon <= -180] + 360                            
    }

    df_cGENIE$lon.mid <- df_cGENIE$lon.mid + rnorm(nrow(df_cGENIE), 0, 0.1)
    df_cGENIE$lat.mid <- df_cGENIE$lat.mid + rnorm(nrow(df_cGENIE), 0, 0.1)

    # Interpolate cGENIE var field to DEM grid for each depth level
    df_interp_all <- foreach(d = unique(df_cGENIE$depth), .combine = rbind, .packages = c('dplyr', 'interp')) %dopar% {
        df_filt <- filter(df_cGENIE, depth == d)

        interp_result <- interp(
            x = df_filt$lon.mid,
            y = df_filt$lat.mid,
            z = df_filt$var,
            xo = df_DEM$lon,
            yo = df_DEM$lat,
            output = "points",
            method = "akima",
            extrap = TRUE
        )

        df_interp <- data.frame(
            lon = as.vector(interp_result$x),
            lat = as.vector(interp_result$y),
            var = as.vector(interp_result$z),
            depth = d
        )

        write.table(df_interp, 
            file = paste0("interp_", var, "_", as.integer(d), "m.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE) 

        df_interp
}
