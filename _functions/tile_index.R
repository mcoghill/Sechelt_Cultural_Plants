#' Tile index.
#'
#' Creates a polygon tile index of a raster.  Attributes of this polygon layer can be used to define \emph{pseudo-tiles} or perhaps better named \emph{virtual-tiles}
#'
#' @param rFile  Character, location of the raster file to index
#' @param pxCount   Specify how many pixels the tile should be in both x and y directions (long. and lat.)
#' @keywords raster, tiles
#' @export
#' @examples
# Create an raster example
#' r1 <- raster::raster(ncol=200, nrow=200, xmn=547000, xmx=552000, ymn=5985000, ymx=5990000)
#' raster::crs(r1) <- "+init=EPSG:3157"
#' r1 <- raster::setValues(r1, seq(1:length(r1)))
#' raster::writeRaster(r1, "e:/tmp/tmp.tif", overwrite = TRUE)
#' # call function
#' tiles <- tile_index("e:/tmp/tmp.tif", pxCount = 50)
#' mapview::mapview(tiles, alpha = 0.5)

# GSIF requires a GDAL object
tile_index <- function(rFile, pxCount) {
  r <- rgdal::GDALinfo(rFile)
  dist <- pxCount * r["res.x"] ## GSIF uses a distance for the tile ... not pixel count
  
  # Generate spatial object and create index
  p_tiles <- GSIF::getSpatialTiles(r, block.x = dist, return.SpatialPolygons = TRUE)
  p_tiles <- sf::st_as_sf(p_tiles) %>% 
    cbind(do.call(rbind, lapply(sf::st_geometry(.) , sf::st_bbox))) %>% 
    dplyr::rename(xl = xmin, yl = ymin, xu = xmax, yu = ymax)
  
  # Suppress further messages
  suppressMessages({
    # Generate associated columns and merge with sf dataframe
    t_tiles <- GSIF::getSpatialTiles(r, block.x = dist, return.SpatialPolygons = FALSE)
    p_tiles <- merge(p_tiles, t_tiles) %>% 
      tibble::rowid_to_column("id")
  })
  
  return(p_tiles)
}
