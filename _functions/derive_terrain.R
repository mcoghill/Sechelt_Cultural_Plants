# Below is the function which generates all SAGA covariate terrain layers used in 
# this project. The first part of this chunk is the XML code for creating the tool 
# chain. Within the function, this only needs to be ran once since it gets changed 
# based on resolution with a much simpler function below it.

# Certain parameters are defined within the toolchain which are important for some 
# of the SAGA layers in order to keep things relatively efficient as well as 
# accurate. These parameters include:

# min_sunrise/max_sunset: Calculations for the minimum sunrise and maximum sunset 
# times over the course of one year. This helps to keep the time for potential 
# incoming solar radiation layer at reasonable levels so that it doesn't take 
# forever creating the layer.

# scale_param: This defines the TRI (ruggedness) search area and channel 
# network length. Essentially, the function scales up the resulting value at 
# low resolutions and scales it down at higher resolutions. Since ruggedness
# looks at an area around the focal cell, it's important to scale this based 
# on resolution since differences will be miniscule at low resolutions, but 
# more dramatic at higher resolutions. This is also used in the "channel 
# network" grid tool to define how long a channel needs to be before it is no 
# longer considered a channel, thus at high resolutions channels need to be 
# continuous for a longer segment before they are considered a channel.

# mrvbf_param: This defines the "initial threshold for slope" parameter in the 
# MRVBF algorithm. The math is based on the paper that this originally came 
# from as well as a very helpful infographic retrieved from here: 
# https://www.nrcs.usda.gov/wps/PA_NRCSConsumption/download?cid=stelprdb1258050&ext=pdf 

# tpi_param: This defines a search radius around a cell. I figured since it has 
# to do with topographic position, search radius should be the same value in 
# cells at each scale. This will create different results at each resolution 
# as well which might help indicate why certain resolutions are more important 
# than others.

# openness_param: This parameter is for using the multi-scaled approach to 
# defining openness. I found through exhaustive testing that this parameter 
# is best defined by half the number of columns, which produced the best 
# looking (least blocky) result out of the plethora of other trials ran.

# The workflow for this function is as follows:
#   I. Perform sanity checks on the function inputs
#  II. Define the tool specific parameters using R language
# III. Build the XML code used for the SAGA toolchain which includes the following tools:

# 1. Fill Sinks in a DEM using the Sink Drainage Route Detection, and subsequent 
# Sink removal tools: For filling sinks in the DEM and only returning a DEM

# 2. Slope, Aspect, Curvature: returns Slope, Aspect, General Curvature, and 
# Total Curvature rasters

# 3. Flow Accumulation (Recursive): Creates the "total catchment area" (TCA) 
# intermediate layer. This layer is only used in the generation of the TWI, 
# and not in the final raster stack.

# 4. Topographic wetness index (TWI): Creates the TWI for the study area. I 
# found these to be the best methods based on exhaustive testing of many other 
# methods. Rather than using the specific catchment area, TCA is used and 
# converted all at once within this tool thus it is more efficient. My methods 
# are based on those found here: 
# https://gracilis.carleton.ca/CUOSGwiki/index.php/Enhanced_Wetness_Modelling_in_SAGA_GIS 

# 5. Channel Network: This creates a channel network grid using the filled DEM 
# and the TCA layer as an initiation grid. I found that using the initiaion 
# value of 1,000,000 worked to create channels consistently at all scales, 
# thus I propose using that here as well. scale_param is used here to define 
# how long a stream needs to be in pixels before it is considered a stream. 
# Source: https://sourceforge.net/projects/saga-gis/files/SAGA%20-%20Documentation/SAGA%20Documents/SagaManual.pdf/download 
# Note: This is an intermediate layer and is not used in the final raster stack

# 6. Overland Flow Distance to Channel Network: This tool draws the overland 
# flow distance to the channel network created in step 5. In order to create 
# a raster that is used to the borders of the grid, the "boundary" option had
# to be set to "true" This may not give a realistic representation of this 
# variable at a given study area, but without this it constricts the final grid.

# 7. Multiresolution Index of Valley Bottom Flatness (MRVBF): Looks at valley 
# bottom flatness and ridge top flatness and uses the mrvbf_param to define 
# the "initial threshold for slope" parameter. The input DEM's used from here 
# all use the original, unfilled DEM.

# 8. Terrain Ruggedness Index (TRI): Looks at ruggedness at a specified distance 
# away from a focal cell.

# 9. Convergence index: This tool is unchanged from the parameters defined by 
# Lucas and Nicholas (uses 3x3 gradient to determind convergence)

# 10. Topographic Openness: Calculates openness across a landscape at multiple
# scales. The advantage to the multi-scale approach is that it calculates the
# values at the edges of the map, whereas the line-tracing method does not. 
# Multi-scale keeps things more consistent I would wager. This outputs both 
# negative and positive openness.

# 11. Diurnal Anisotropic Heat: The settings are unchanged from the default.

# 12. Topographic Position Index (TPI): This calculates the position of a cell 
# relative to neighboring cells in a defined vicinity. This value could 
# easily change depending on the scale of interest. When running this tool, 
# it pastes the default search parameters; however, it actually runs the ones 
# defined by the tpi_param code. This is a bug from SAGA, but has no effect 
# on the outcome of the TPI, it outputs it properly.

# 13. Potential Incoming Solar Radiation: This calculates direct and diffuse 
# insolation from a range of days spanning a whole year. A customized function
# uses the data from https://www.nrc-cnrc.gc.ca/eng/services/sunrise/advanced.html
# to calculate the minimum sunrise and maximum sunset times at a given point 
# throughout the span of one year (rounded to the nearest 30 minutes).

# 14. Change Grid Value: Changes the grid values of 0 to -99999 (i.e.: no data) 
# for the insolation grids

# 15. Export Raster: The final output(s)

# This function requires that you specify an input file DEM either as a raster 
# file or a file path to a raster file. The files will be written in a path either 
# specified or relative to the input DEM.

dem_derived_layers <- function(
  dem_input, # Input DEM either as file path or raster/rast layer
  out_dir,  # output file path, will default to dem_input directory
  saga_cmd # path to saga_cmd
) {
  # Note: All SAGA tools and associated documentation can be found here:
  # http://www.saga-gis.org/saga_tool_doc/7.3.0/index.html
  
  # Data input checks/error handling
  if(missing(dem_input)) stop("dem_input is missing with no default")
  if(missing(saga_cmd)) stop("A path to saga_cmd is missing with no default")
  
  if(class(dem_input) %in% "SpatRaster") {
    reference <- dem_input
    dem_input <- terra::sources(reference)[, "source"]
  } else if(base::class(dem_input) %in% "RasterLayer") {
    reference <- terra::rast(dem_input)
    dem_input <- terra::sources(reference)[, "source"]
  } else if(base::is.character(dem_input) && base::length(dem_input) == 1) {
    reference <- terra::rast(dem_input)
  } else {
    stop(paste(
      "\rError: 'dem_input' must either be a character string to a valid raster 
      \rfile, or a single raster layer from either the 'raster' or 'terra' package"
    ))
  }
  
  if(missing(out_dir)) out_dir <- base::dirname(dem_input)[1]
  if(!base::is.character(out_dir)) stop("out_dir is not a valid directory")
  if(base::length(out_dir) != 1) warning("Multiple strings passed to out_dir, only the first one will be used")
  out_dir <- out_dir[1]
  base::dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  if(!base::file.exists(saga_cmd)) stop("Please provide a valid path to 'saga_cmd'")
  
  # Define inputs used in some of the SAGA tools below
  reference_res <- base::max(terra::res(reference))
  base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Sunrise/sunset time calculation for solar radiation tool:
  base::suppressWarnings(base::suppressMessages({
    centroid <- data.frame(x = mean(c(terra::xmax(reference), terra::xmin(reference))),
                           y = mean(c(terra::ymax(reference), terra::ymin(reference)))) %>%
      sf::st_as_sf(coords = c("x", "y"), crs = crs(reference)) %>%
      sf::st_transform(4326) %>% # needs lat/long for calculating sunlight times
      sf::st_coordinates() %>%
      as.data.frame() %>% 
      dplyr::rename_at(vars(names(.)), ~c("Lon", "Lat"))
    
    # Using rvest to gather sunset/sunrise information
    # Use previous year to calculate full year of data
    yr <- base::as.numeric(format(Sys.Date(), "%Y")) - 1
    session <- rvest::html_session(
      "https://www.nrc-cnrc.gc.ca/eng/services/sunrise/advanced.html")
    form <- rvest::html_form(session)
    form <- form[lapply(form, function(x) x[["name"]]) == "sunrise-sunset"][[1]]
    filled_form <- rvest::set_values(
      form, calcSubject = "ss-year-txt", calcMethod = "byLatLong",
      calcDate = yr, olat = centroid$Lat, elong = centroid$Lon, timezone = "PST|8",
      latitude = "North", longitude = "West")
    filled_form$url <- session$url
    
    result <- rvest::submit_form(session, filled_form)
    
    # Remove everything up to and including the first line break (encoded as 0a)
    # This makes the table making much simpler
    lb <- base::which(result$response$content == "0a")[1] + 1
    result$response$content <- result$response$content[lb : length(result$response$content)]
    
    response <- httr::content(result$response, type = "text/csv") %>% 
      na.omit() %>%
      dplyr::mutate(Month = base::match(trimws(gsub("[[:digit:]]+|/", "", Date)), month.abb),
                    Day = readr::parse_number(Date)) %>%
      dplyr::mutate(rise = lubridate::as_datetime(paste0(
        yr, "-", Month, "-", Day, " ", lubridate::hms(`Sun rise`)), tz = "Canada/Pacific"),
        set = lubridate::as_datetime(paste0(
          yr, "-", Month, "-", Day, " ", lubridate::hms(`Sun set`)), tz = "Canada/Pacific")) %>%
      dplyr::mutate(rise = lubridate::floor_date(rise, "30 minutes"),
                    set = lubridate::ceiling_date(set, "30 minutes")) %>%
      dplyr::mutate(rise = base::strftime(rise, format = "%H:%M:%S"),
                    set = base::strftime(set, format = "%H:%M:%S"))
  }))
  
  min_sunrise <- min(response$rise)
  max_sunset <- max(response$set)
  min_sunrise_numeric <- as.numeric(lubridate::hms(min_sunrise)) / 3600
  max_sunset_numeric <- as.numeric(lubridate::hms(max_sunset)) / 3600
  
  time_res <- (max_sunset_numeric - min_sunrise_numeric) / 4
  
  ## Scaling parameter for a few of the tools:
  scale_param <- ifelse(
    round((reference_res*(25 / reference_res)) / reference_res) == 0, 1, 
    round((reference_res*(25 / reference_res)) / reference_res)
  )
  
  ## Tool specific parameters
  mrvbf_param <- 116.57 * (reference_res ^ -0.62)
  tpi_param <- reference_res * 5
  openness_param <- as.numeric(ncol(reference) / 2)
  
  #############################################################################
  # Build the xml toolchain piece by piece (easier than previously)
  
  # Does this by first creating a list of tools and their associated parameters
  # If a tool is not going to be ran, it can simply be commented out/deleted, and
  # that's it, no mucking around anywhere else! 
  
  # Note: There are multiple ways to fill sinks for a given DEM. The way 
  # currently coded here may not be the most efficient but it produces a 
  # sinkroute layer which is used in multiple tools later on and appears to be 
  # the most accurate.
  tool_list <- list(
    fill_sinks = paste0(
      "<tool library='ta_preprocessor' tool='1' name='Sink Drainage Route Detection'>
        <input  id='ELEVATION'>dem</input>
        <output id='SINKROUTE'>sinkroute</output>
    </tool>
    <tool library='ta_preprocessor' tool='2' name='Sink Removal'>
        <input  id='DEM'>dem</input>
        <input  id='SINKROUTE'>sinkroute</input>
        <output id='DEM_PREPROC'>dem_preproc</output>
        <option id='METHOD'>1</option>
        <option id='THRESHOLD'>0</option>
    </tool>"
    ),
    
    # Can script more output curvature types, see:
    # http://www.saga-gis.org/saga_tool_doc/7.3.0/ta_morphometry_0.html
    slope = paste0(
      "<tool library='ta_morphometry' tool='0' name='Slope, Aspect, Curvature'>
        <input id='ELEVATION'>dem_preproc</input>
        <output id='SLOPE'>slope</output>
        <output id='ASPECT'>aspect</output>
        <output id='C_GENE'>gencurve</output>
        <output id='C_TOTA'>totcurve</output>
        <option id='METHOD'>6</option>
        <option id='UNIT_SLOPE'>0</option>
        <option id='UNIT_ASPECT'>0</option>
    </tool>"
    ),
    
    #   Following this method for calculating topographic wetness index:
    #    https://gracilis.carleton.ca/CUOSGwiki/index.php/Enhanced_Wetness_Modelling_in_SAGA_GIS
    #    See this paper as well for discussion on different ways to calculate TWI:
    #    https://link.springer.com/article/10.1186/s40965-019-0066-y
    flow_accum = paste0(
      "<tool library='ta_hydrology' tool='1' name='Flow Accumulation (Recursive)'>
         <input id='ELEVATION'>dem_preproc</input>
         <output id='FLOW'>tca</output>
     </tool>"
    ),
    
    
    twi = paste0(
      "<tool library='ta_hydrology' tool='20' name='Topographic Wetness Index (TWI)'>
      <input id='SLOPE'>slope</input>
      <input id='AREA'>tca</input>
      <output id='TWI'>twi</output>
      <option id='CONV'>1</option>
      <option id='METHOD'>1</option>
  </tool>"
    ),
    
    # The threshold value of 1,000,000 on the TCA grid appears to give a good
    # representation of the channel network. See
    # https://sourceforge.net/projects/saga-gis/files/SAGA%20-%20Documentation/SAGA%20Documents/SagaManual.pdf/download
    chan_net = paste0(
      "<tool library='ta_channels' tool='0' name='Channel Network'>
       <input id='ELEVATION'>dem</input>
       <input id='SINKROUTE'>sinkroute</input>
       <output id='CHNLNTWRK'>cnetwork</output>
       <input id='INIT_GRID'>tca</input>
       <option id='INIT_METHOD'>2</option>
       <option id='INIT_VALUE'>1000000</option>
       <option id='DIV_CELLS'>5</option>
       <option id='MINLEN'>", scale_param, "</option>
   </tool>"
    ),
    
    # Since multiple flow direction (MFD) algorithms were used throughout, I
    # thought it would be consistent to keep that trend in this tool here
    # (METHOD = 1). Keeping the "BOUNDARY" attribute at TRUE will model the flow
    # distance to the boundaries of the study area (by making the study area
    # boundary a "channel" in itself).
    
    flow_dist = paste0(
      "<tool library='ta_channels' tool='4' name='Overland Flow Distance to Channel Network'>
            <input id='ELEVATION'>dem_preproc</input>
            <input id='CHANNELS'>cnetwork</input>
            <output id='DISTANCE'>hdist</output>
            <output id='DISTVERT'>vdist</output>
            <option id='METHOD'>1</option>
            <option id='BOUNDARY'>true</option>
        </tool>"
    ),
    
    # See page 2-3 of:
    # https://www.nrcs.usda.gov/wps/PA_NRCSConsumption/download?cid=stelprdb1258050&ext=pdf
    # Which makes reference to the following article:
    # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2002WR001426
    # So initial slope threshold needs to be changed depending on the input
    # resolution, which is defined by the mrvbf_param variable. All other values
    # are SAGA default
    
    mrvbf = paste0(
      "<tool library='ta_morphometry' tool='8' name='Multiresolution Index of Valley Bottom Flatness (MRVBF)'>
           <input id='DEM'>dem</input>
           <output id='MRVBF'>mrvbf</output>
           <output id='MRRTF'>mrrtf</output>
           <option id='T_SLOPE'>", mrvbf_param, "</option>
           <option id='T_PCTL_V'>0.400000</option>
           <option id='T_PCTL_R'>0.350000</option>
           <option id='P_SLOPE'>4.000000</option>
           <option id='P_PCTL'>3.000000</option>
           <option id='UPDATE'>false</option>
           <option id='CLASSIFY'>false</option>
           <option id='MAX_RES'>100.000000</option>
       </tool>"
    ),
    
    # Terrain Ruggedness Index (TRI): Looks at ruggedness at a specified distance
    # away from a focal cell. The scale_param will adjust the search radius
    # distance in order to get a consistent idea across scales of the terrain
    # ruggedness. Ex: at 2.5m, it should look at a distance of 25m (i.e.: radius
    # of 10 cells), while at 5m if we want to look at a 25m distance away, it
    # should have a radius of 5 cells, etc.
    
    tri = paste0(
      "<tool library='ta_morphometry' tool='16' name='Terrain Ruggedness Index (TRI)'>
           <input id='DEM'>dem</input>
           <output id='TRI'>tri</output>
           <option id='MODE'>0</option>
           <option id='RADIUS'>", scale_param, "</option>
           <option id='DW_WEIGHTING'>0</option>
        </tool>"
    ),
    
    convergence = paste0(
      "<tool library='ta_morphometry' tool='1' name='Convergence Index'>
           <input id='ELEVATION'>dem</input>
           <output id='RESULT'>convergence</output>
           <option id='METHOD'>1</option>
           <option id='NEIGHBOURS'>1</option>
       </tool>"
    ),
    
    # The default openness parameters create a fairly blocky grid. With the Eagle
    # Hills data, changing the radius parameter to be half the amount of rows in
    # a grid helped to alleviate that. I don't have documentation for that, that
    # was just me trying things out and seeing what worked. The original paper for
    # calculating openness is:
    # https://www.asprs.org/wp-content/uploads/pers/2002journal/march/2002_mar_257-265.pdf
    
    openness = paste0(
      "<tool library='ta_lighting' tool='5' name='Topographic Openness'>
           <input id='DEM'>dem</input>
           <output id='POS'>open_pos</output>
           <output id='NEG'>open_neg</output>
           <option id='RADIUS'>", openness_param, "</option>
           <option id='METHOD'>0</option>
           <option id='DLEVEL'>3.000000</option>
           <option id='NDIRS'>8</option>
      </tool>"
    ),
    
    # These are default parameters
    dah = paste0(
      "<tool library='ta_morphometry' tool='12' name='Diurnal Anisotropic Heat'>
          <input id='DEM'>dem</input>
          <output id='DAH'>dah</output>
          <option id='ALPHA_MAX'>202.500000</option>
      </tool>"
    ),
    
    # The tpi_param here will calculate topographic position based on its position
    # relative to the cells in its neighborhood. For a 25m grid, 5 cells away is
    # 125m, while for a 2.5m grid 5 cells away is 12.5m This was scaled this way
    # to promote faster processing at higher resolutions. I also think this is a
    # better representation of topographic position at different scales anyway.
    
    tpi = paste0(
      "<tool library='ta_morphometry' tool='18' name='Topographic Position Index (TPI)'>
          <input id='DEM'>dem</input>
          <output id='TPI'>tpi</output>
          <option id='STANDARD'>false</option>
          <option id='RADIUS'>0.000000;", tpi_param, "</option>
          <option id='DW_WEIGHTING'>0</option>
      </tool>"
    ),
    
    # Using the rvest package, the centroid of the raster is found and sent off to
    # https://www.nrc-cnrc.gc.ca/eng/services/sunrise/advanced.html, where we get
    # official accurate representations of daily sunrise and sunset times
    # throughout a year at that centroid. Minimum sunrise and maximum sunset times
    # are calculated for a given study area in order to limit the amount of
    # calculations that occur while there would be no solar radiation at any given
    # time throughout a year. It also helps to decrease processing time.
    
    solar = paste0(
      "<tool library='ta_lighting' tool='2' name='Potential Incoming Solar Radiation'>
           <input id='GRD_DEM'>dem</input>
           <output id='GRD_DIRECT'>direinso</output>
           <output id='GRD_DIFFUS'>diffinso</output>
           <option id='SOLARCONST'>1367.000000</option>
           <option id='LOCALSVF'>1</option>
           <option id='UNITS'>0</option>
           <option id='SHADOW'>0</option>
           <option id='LOCATION'>1</option>
           <option id='PERIOD'>2</option>
           <option id='DAY'>2019-01-15</option>
           <option id='DAY_STOP'>2019-12-15</option>
           <option id='DAYS_STEP'>30</option>
           <option id='HOUR_RANGE'>", min_sunrise_numeric, "; ", max_sunset_numeric, "</option>
           <option id='HOUR_STEP'>", time_res, "</option>
           <option id='METHOD'>2</option>
           <option id='LUMPED'>70.000000</option>
           <option id='UPDATE'>0</option>
       </tool>
  
       <tool library='grid_tools' tool='12' name='Change Grid Values'>
         <input id='INPUT'>direinso</input>
         <option id='METHOD'>0</option>
         <option id='IDENTITY'>
           <OPTION type='static_table' id='IDENTITY' name='Lookup Table'>
             <FIELDS>
               <FIELD type='DOUBLE'>New Value</FIELD>
               <FIELD type='DOUBLE'>Value</FIELD>
             </FIELDS>
             <RECORDS>
               <RECORD>
                 <FIELD>-99999.000000</FIELD>
                 <FIELD>0.000000</FIELD>
               </RECORD>
             </RECORDS>
           </OPTION>
         </option>
       </tool>
  
     <tool library='grid_tools' tool='12' name='Change Grid Values'>
         <input id='INPUT'>diffinso</input>
         <option id='METHOD'>0</option>
         <option id='IDENTITY'>
           <OPTION type='static_table' id='IDENTITY' name='Lookup Table'>
             <FIELDS>
               <FIELD type='DOUBLE'>New Value</FIELD>
               <FIELD type='DOUBLE'>Value</FIELD>
             </FIELDS>
             <RECORDS>
               <RECORD>
                 <FIELD>-99999.000000</FIELD>
                 <FIELD>0.000000</FIELD>
               </RECORD>
             </RECORDS>
           </OPTION>
         </option>
       </tool>"
    )
  )
  
  # End list of tools, don't adjust code after here
  #############################################################################
  
  tools <- lapply(tool_list, function(x) {
    inputs = gsub(".*>", "", regmatches(
      x, gregexpr("(?<=<input).*?(?=</input>)", x, perl = TRUE))[[1]])
    outputs = unique(gsub(".*>", "", regmatches(
      x, gregexpr("(?<=<output).*?(?=</output>)", x, perl = TRUE))[[1]]))
    
    inputs = unique(inputs[!inputs %in% outputs])
    tc = sum(length(inputs), length(outputs))
    
    tibble(tool = x, inputs = list(inputs), outputs = list(outputs), tc = tc)
  })
  
  # Dynamically split processing based on the amount of available RAM a PC has.
  # The following code will determine the best way to process the rasters
  # based on the amount of RAM on a users machine. It does this by summing the
  # number of inputs and outputs and comparing them to the maximum number of 
  # files your system can hold based on the size of the input DEM.
  raster_size <- file.size(dem_input) / 1024 ^ 2
  memory <- memory.limit() - raster_size - 1024 # Leave ~1GB of RAM available?
  n_files <- floor(memory / raster_size)
  
  # Adjust tools list
  for(i in 1:length(tools)) {
    try({
      if(tools[[i]]$tc < n_files) {
        while(sum(tools[[i]]$tc, tools[[i + 1]]$tc) < n_files) {
          y <- rbind(tools[[i]], tools[[i + 1]])
          y[1, "tc"] <- sum(y$tc)
          y[1, "tool"] <- paste(y$tool, collapse = "\n")
          y[1, "inputs"][[1]] <- list(unlist(y$inputs))
          y[1, "outputs"][[1]] <- list(unlist(y$outputs))
          y <- y[1, ]
          tools[[i]] <- y
          tools[[i + 1]] <- NULL
        }
      }
    }, silent = TRUE)
  }
  
  # Fix inputs and outputs columns, add function to call the whole thing in cmd
  xml_layout <- lapply(tools, function(x) {
    x$inputs <- list(
      unique(unlist(x$inputs)[!unlist(x$inputs) %in% unlist(x$outputs)])
    )
    
    x$out_files <- foreach(k = x$outputs) %do% {
      tempfile(pattern = paste0("spat_", k, "_"), fileext = ".tif")
    }
    
    x$input_xml <- foreach(i = x$inputs, .combine = paste) %do% {
      paste0("<input varname='", i, "' type='grid' parent='GRID_SYSTEM'>
          <name>", i, "</name>
      </input>\n", collapse = " ")}
    
    x$output_xml <- foreach(k = 1:length(unlist(x$outputs)), .combine = paste) %do% {
      paste0(
        "<tool library='io_gdal' tool='2' name='Export GeoTIFF'>
        <input id='GRIDS'>", unlist(x$outputs)[k], "</input>
        <option id='FILE'>", unlist(x$out_files)[k], "</option>
        
      </tool>\n", collapse = " ")
    }
    x$header <- paste0(
      "<?xml version='1.0' encoding='UTF-8'?>
      <toolchain saga-version='7.3.0'>
      <group>toolchains</group>
      <identifier>Derived</identifier>
      <name>Derived Layers (one step)</name>
      <description>
        Common DEM derivatives in SAGA GIS
      </description>
    
      <parameters>
        <option varname='GRID_SYSTEM' type='grid_system'>
          <name>Grid System</name>
        </option>
        ", x$input_xml,
      "</parameters>
      <tools>", sep = "\n")
    
    x$footer <- paste0(
      "</tools>
  </toolchain>"
    )
    
    x$call <- paste(x$header, x$tool, x$output_xml, x$footer, sep = "\n")
    return(x)
  })
  
  # Define text for cmd input
  cmd_text <- lapply(xml_layout, function(x) {
    foreach(k = unlist(x$inputs), .combine = paste) %do% {
      if(k != "dem") {
        k_out_id <- unlist(sapply(xml_layout, function(x) which(x$outputs[[1]] == k)))
        k_dir <- xml_layout[[names(k_out_id)]]$out_files[[1]][k_out_id]
        paste0("-", k, " ", k_dir)
      } else {
        paste0("-dem ", dem_input)
      }
    }
  })
  
  # Determine the toolchain directory based on your system
  dem_derived_xml <- ifelse(
    Sys.info()[["sysname"]] == "Windows", 
    file.path(dirname(saga_cmd), "tools", "toolchains", "dem_derived_layers.xml"), 
    file.path(dirname(dirname(saga_cmd)), "share", "saga", 
              "toolchains", "dem_derived_layers.xml")
  )
  
  # Process SAGA toolchains in the list
  for(p in 1:length(xml_layout)) {
    write_xml(read_xml(xml_layout[[p]]$call), dem_derived_xml)
    sys_cmd <- paste("toolchains Derived", cmd_text[[p]])
    system2(saga_cmd, sys_cmd)
  }
  
  # Remove intermediate files. Files aren't deleted now, because everything 
  # is saved as a temp file. Rather, files are selected for final processing
  # that don't match the id's you specify
  rem <- lapply(xml_layout, function(x) {
    rem_id <- which(unlist(x$outputs) %in% c("cnetwork", "tca", "dem_preproc", "sinkroute"))
    out_files <- unlist(x$out_files)[!unlist(x$out_files) %in% unlist(x$out_files)[rem_id]]
    out_names <- unlist(x$outputs)[!unlist(x$outputs) %in% unlist(x$outputs)[rem_id]]
    list(out_files = out_files, out_names = out_names)
  })
  out_list <- foreach(i = 1:length(rem), .combine = rbind) %do% {
    data.frame(out_files = rem[[i]]$out_files, 
               out_names = rem[[i]]$out_names)
  }
  
  out <- terra::rast(out_list$out_files) %>% magrittr::set_names(out_list$out_names)
  crs(out) <- terra::crs(reference)
  out <- terra::writeRaster(out, overwrite = TRUE, 
                            filename = file.path(out_dir, paste0(names(out), ".tif")))
  
  return(out)
}
