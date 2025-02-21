#' cGENIE DEM Interpolation
#'
#' This function interpolates a cGENIE 3D output field to the points of a
#' digital elevation model (DEM). If a point is outside the range of cGENIE
#' mid-layer depths, it is assigned the value from the closest valid top/bottom
#' layer cell. If a point is within that range, a value is interpolated using
#' the two closest valid cells immediately above and below it. This function
#' returns a dataframe 'df_DEM' and saves a plot.
#' 
#' @param var A string specifying the variable name to extract from the NetCDF
#'              file. Examples include "ocn_temp", "ocn_O2", etc.
#' @param experiment A string specifying the path to the folder or file prefix
#'                   where the experiment's data is stored.
#' @param dem A string specifying the path to the DEM NetCDF file.
#' @param zname A string specifying the name of the depth variable in the
#'                 DEM file.
#'
#' @return A dataframe (`df.sum`) containing:
#'   - `lon`: Longitude of DEM points.
#'   - `lat`: Latitude of DEM points.
#'   - `depth`: Depth at DEM points.
#'   - `var`: Interpolated tracer values (e.g., [O2])
#'
#' @import RNetCDF, dplyr, reshape2, ggplot2, RNetCDF, sf
#' @export

cGENIE.DEM.interp <- function(experiment, var, dem, zname){

    library(dplyr)
    library(reshape2)
    library(PaleoClimR)
    library(ggplot2)
    library(RNetCDF)
    library(sf)

    # Extract cGENIE data
    df_cGENIE <- cGENIE.data.3D(var, experiment, "default", "biogem")
    depths_cGENIE <- as.data.frame(sort(unique(df_cGENIE$depth)))

    # Extract DEM data and compile to dataframe
    nc <- open.nc(dem)
    lat <- var.get.nc(nc, "lat")
    lon <- var.get.nc(nc, "lon")
    depth <- var.get.nc(nc, "z")
    close.nc(nc)
    depth <- depth * -1
    df_DEM <- as.data.frame(
        cbind(
            rep(lon, times = length(lat), each = 1),
            rep(lat, times = 1, each = length(lon)),
            as.data.frame(melt(depth))$value
        )
    )
    names(df_DEM) <- c("lon", "lat", "depth")

    # Filter out land from DEM dataframe
    df_DEM <- df_DEM %>%
        filter(depth >= 0)

    # Adjust DEM grid for longitude projection (to 0-360 degrees if needed)
    if (mean(between(df_DEM$lon, -180, 180)) < 1) {
        df_DEM$lon[df_DEM$lon <= -180] <- df_DEM$lon[df_DEM$lon <= -180] + 360                            
    }

    # Process cGENIE dataframe into three objects:

    #   srf_cGENIE: cells in top layer with valid value of var
    srf_cGENIE <- df_cGENIE %>%
        filter(depth == min(df_cGENIE$depth) & !is.na(var))

    #   ben_cGENIE: cells in bottom layer with valid value of var
    ben_cGENIE <- df_cGENIE %>%
        filter(depth == max(df_cGENIE$depth) & !is.na(var))

    #   layers_cGENIE: list of dataframes, each containing two consecutive
    #   depth layers, with cells where var is valid across both layers
    layers_cGENIE <- lapply(1:(length(unique(df_cGENIE$depth))-1), function(i) {
        depths <- sort(unique(df_cGENIE$depth))
        d1 <- depths[i]
        d2 <- depths[i+1]
        
        df_cGENIE %>%
            # Filter to just these two depths
            filter(depth %in% c(d1, d2)) %>%
            # Group by lat/lon pairs
            group_by(lon.mid, lat.mid) %>%
            # Keep only pairs where both depths have non-NA values
            filter(n() == 2, !any(is.na(var))) %>%
            ungroup()
    })

    names(layers_cGENIE) <- sapply(1:(length(unique(df_cGENIE$depth))-1), function(i) {
        depths <- sort(unique(df_cGENIE$depth))
        sprintf("%g", depths[i])
    })

    # Assign a value of var for each DEM point
    df_DEM <- df_DEM %>%
        rowwise() %>%
        mutate(
            var_interp = case_when(
                # Shallow case:
                depth <= min(df_cGENIE$depth) ~ {
                    layer <- srf_cGENIE %>%
                        # Assign distance from point to cells
                        mutate(dist = sqrt((lon.mid - lon)^2 + (lat.mid - lat)^2)) %>%
                        # Slice closest point
                        slice_min(dist, n = 1) %>%
                        # Pull var from closest point
                        pull(var)
                },
                
                # Deep case:
                depth >= max(df_cGENIE$depth) ~ {
                    layer <- ben_cGENIE %>%
                        # Assign distance from point to cells
                        mutate(dist = sqrt((lon.mid - lon)^2 + (lat.mid - lat)^2)) %>%
                        # Slice closest point
                        slice_min(dist, n = 1) %>%
                        # Pull var from closest point
                        pull(var)
                },

                # Intermediate case:
                TRUE ~ {
                    # Find name of correct layer in layers_cGENIE
                    target_depth <- sprintf("%g", max(filter(depths_cGENIE, depths_cGENIE < depth)))
                    layer <- layers_cGENIE[[target_depth]]
                    if (is.null(layer)){
                        NA
                    }
                    else {
                    # Assign distance from point to cells
                    mutate(layer, dist = sqrt((lon.mid - lon)^2 + (lat.mid - lat)^2)) %>%
                        # Slice closest point
                        slice_min(dist, n = 1)
                    # Interpolate using the closest point
                    approx(
                                x = layer$depth,
                                y = layer$var,
                                xout = depth
                    )$y
                }
                }
            )
        ) %>%
        ungroup()

    # Make plot
    ggplot(df_DEM, aes(x = lon, y = lat, fill = var_interp)) +
        geom_raster() +
        geom_rect( 
                aes(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
                fill = alpha("grey", 0), linewidth = 0.1, colour = "black"
                ) +
        scale_fill_viridis_c(limits = c(0, 0.0003), oob = scales::squish, name=var) +
        coord_quickmap() +
        theme_bw()
    ggsave(paste(experiment, "interp", var, ".png", sep="_" ))

    return(df_DEM)
}
