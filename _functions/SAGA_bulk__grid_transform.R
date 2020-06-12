# This document isn't really reproducable, but the idea is to use SAGA to transform multiple input rasters
# all at once. SAGA is used to transform the grids, those grids are saved, and then the CRS is set using the
# raster package in R.
saga_path <- file.path("C:/SAGA/saga_cmd.exe")
if (Sys.info()['sysname'] == "Windows") {
  saga_cmd = saga_path
} else {
  saga_cmd = "saga_cmd"
}
system(paste(saga_cmd, "-v"))



#batch1 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1961_1990S", full.names = TRUE)[1:20]
#batch7 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1961_1990S", full.names = TRUE)[21:40]
#batch7 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1961_1990S", full.names = TRUE)[41:56]
#batch7 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1961_1990Y", full.names = TRUE)
#batch7 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1981_2010S", full.names = TRUE)[1:20]
#batch7 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1981_2010S", full.names = TRUE)[21:40]
#batch7 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1981_2010S", full.names = TRUE)[41:56]
#batch7 <- list.files("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/Normal_1981_2010Y", full.names = TRUE)
dem <- file.path("./Sechelt_AOI/1_map_inputs/covariates/4m/dem.tif")


batch_transform_xml <- paste(deparse(substitute(paste0(
  "<?xml version='1.0' encoding='UTF-8'?>
<toolchain saga-version='7.4.0'>
<group>toolchains</group>
<identifier>transform</identifier>
<name>grid_transform</name>
<description>Satellite Image Processing for PEM</description>

<parameters>
  <input varname='INPUT' type='grid_list'>
    <name>Input Grids</name>
  </input>
  <option varname='REFERENCE_GRID_SYSTEM' type='grid_system'>
    <name>Reference Grid System</name>
  </option>
  <input varname='DEM_REFERENCE' type='grid' parent='REFERENCE_GRID_SYSTEM'>
    <name>Input DEM</name>
  </input>
  </parameters>

<tools>
  <tool library='pj_proj4' tool='0' name='Set Coordinate Reference System'>
    <option id='CRS_PROJ4'>", st_crs(4326)$proj4string, "</option>
    <input id='GRIDS'>INPUT</input>
</tool>

    <tool library='pj_proj4' tool='3' name='Coordinate Transformation (Grid List)'>
      <option id='CRS_METHOD'>0</option>
      <input id='CRS_GRID.PICK'>DEM_REFERENCE</input>
      <input id='SOURCE'>INPUT</input>
      <input id='TARGET_TEMPLATE'>DEM_REFERENCE</input>
      <option id='TARGET_DEFINITION'>1</option>
      <option id='RESAMPLING'>3</option>
      <option id='KEEP_TYPE'>false</option>
      <output id='GRIDS'>TRANSFORM_LIST</output>
    </tool>",
  mask_xml,
  "</tools>
</toolchain>"))), collapse = '')

mask_xml <- NULL
for(j in 1:length(batch7)){
  k <- basename(file_path_sans_ext(batch7[j]))
  mask_xml <- paste(mask_xml, paste0("
  <tool library='grid_tools' tool='32' name='Select Grid from List'>
    <input id='GRIDS'>TRANSFORM_LIST</input>
    <output id='GRID' parent='REFERENCE_GRID_SYSTEM'>", k, "</output>
    <option id='INDEX'>", j-1, "</option>
  </tool>
  <tool library='io_gdal' tool='2' name='Export GeoTIFF'>
        <input id='GRIDS'>", k, "</input>
        <option id='FILE'>", normalizePath(file.path("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/test", paste0(k, ".tif"))), "</option>
      </tool>"), collapse = ' ')
}

write_xml(read_xml(eval(parse(text = batch_transform_xml))), file.path(dirname(saga_path), "tools", "toolchains", "Transform.xml"))
sys_cmd <- paste(saga_cmd, "toolchains transform", 
                 "-DEM_REFERENCE", dem,
                 "-INPUT", paste(batch7, collapse = ";")
)
system(sys_cmd)

file_list <- list.files(file.path("F:/Sechelt/ClimateBC/DEM_Inputs/DEM_4m_Cropped/test"), full.names = TRUE)
for(i in file_list) {
  rast <- raster(i)
  crs(rast) <- crs(raster(dem))
  writeRaster(rast, file.path("./Sechelt_AOI/1_map_inputs/covariates/4m", paste0("normal_1981_2010y_", basename(i))))
}
rm(rast)
gc()
unlink(file_list)
