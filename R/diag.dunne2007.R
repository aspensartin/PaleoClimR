# organic carbon export flux to burial flux scheme
# after Dunne, J.P., Sarmiento, J.L., and Gnanadesikan, A. (2007) 
# Global Biogeochem. Cy., 21, GB4006, doi:10.1029/2006GB002907

diag.dunne2007 <- function(export){
    # convert export flux from μmol m^-2 yr^-1 to mmol m^-2 d^-1
    export <- export / 365.25e3
    # apply eq. 3 Dunne et al. (2007)
    burial <- export * (0.013 + 0.53 * export^2 / (7.0 + export)^2)
    return(burial)
}

