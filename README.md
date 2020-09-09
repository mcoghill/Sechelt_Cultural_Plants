# Sechelt_Cultural_Plants

This repository contains coding used to create predictive maps of plant presence and absence for culturally important plants to the Sechelt Nation. The study area of interest is of the Sechelt Spipiyus swiya, Northwest of Sechelt, BC. This project is headed by Sari Saunders, Heather Klassen, Michael Ryan, and members of the Sechelt Nation. My role is the modeler, where I build the predictive maps using the code in this repository. Below, I will attempt to explain what each piece of the coding puzzle is and does.

## 01: Spatial Data Preparation Stage
### [01a_DEM_SpatialLayer_Prep](https://github.com/mcoghill/Sechelt_Cultural_Plants/blob/master/01a_DEM_SpatialLayer_Prep.Rmd)
In this script, the initial DEM and CHM files are properly transformed to the BC Albers projection and resampled to the desired resolutions. These files were originally provided by the Sechelt team and due to the size of files as well as agreements we have, they will not be uploaded to GitHub. Then, the TRIM elevation is downloaded as well for comparison's sake. Using [```SAGA GIS```](https://sourceforge.net/projects/saga-gis/), multiple covariate rasters are derived from elevation. Depending on your system and the desired resolutions you produce, this can take a fairly long time. For this project, a 4m resolution was chosen to do all of the modelling.
It's also important to note here that I exported DEM's as a .asc file in the Lat/Long projection (EPSG 4326). This was done in order to create climate layers from [ClimateBC](http://climatebc.ca/) software. There is no CLI for this program, so using the point-and-click software was necessary. Once the layers are created, they are saved as .tif files in the BC Albers projection.

### [01b_BCData_Layer_Prep](https://github.com/mcoghill/Sechelt_Cultural_Plants/blob/master/01b_BCDataLayer_Prep.Rmd)
This file uses the newly developed [```bcdata```](https://github.com/bcgov/bcdata) package to download spatial layers that are used for modelling later on. These layers include BEC zones, TEM, VRI, cutblocks, roads, waterbodies, and fire polygons. When a layer is initially downloaded, it comes in as the bounding box of our area of interest (AOI). Then, it uses computing power to trim the box to the AOI. In earlier iterations, I tried to download the area of interest and skip a step; however, there were too many vertices on the AOI polygon and it resulted in an error. Ultimately, this script doesn't take too long to run anyways so it doesn't matter which way it completes the process.

### [01c_Satellite_SpatialLayer_Prep](https://github.com/mcoghill/Sechelt_Cultural_Plants/blob/master/01c_Satellite_SpatialLayer_Prep.Rmd)
At this point, the raster datasets include terrain derived and climate based layers. One missing piece would be satellite imagery, which is freely available at a 10m resolution if you know how to access it. The way that I do this uses an interface to [Googel Earth Engine Code](https://code.earthengine.google.com/). You need a Google account in order to use it and it's free to sign up and use. This open source tool is able to access satellite imagery from a variety of satellite sources. In our case, we are only interested in Sentinel-2 imagery since it has the highest naturally available resolution of all available satellites. 

Google Earth Engine code uses python scripting in order to process any request. The [```reticulate```](https://github.com/rstudio/reticulate) package is able to download and use python and associated packages for our purposes here. Once that is completed, we then ask Google Earth Engine to composite images from a range of dates (e.g.: get all images from 2019), remove clouds from all pixels where there are clouds detected, and then take the median pixel value from the remaining composite. We want each of the bands from the processed median image, and those get downloaded to our Google Drive folder (there is no way to directly download those files). Once they are in our Google Drive folder, we can download them one at a time to continue processing.

Once the bands are downloaded, they get converted to satellite indices using the [```RStoolbox```](https://github.com/bleutner/RStoolbox) package. These satellite indices (NDVI, MNDWI, etc.) are saved into the same folder as the rest of the covariates and are what gets used for modelling (we do NOT use the raw bands, that is an incorrect process).

**Note: This script hasn't been run in a while and errors may occur because of that. 

## 02: Point Data Preparation and Modelling
### [02a_PointData_Prep](https://github.com/mcoghill/Sechelt_Cultural_Plants/blob/master/02a_PointData_Prep.Rmd)
At this stage, we want to gather all of the available points that we can use for modelling purposes. In 2019, random points within the Sechelt Spipiyus swiya were generated and data was collected at 84 total points. The types of data collected included site specific ecological data (site series, structural stage, canopy cover, dominant plants, etc.), as well as presence/absence and cover data for each of the plant species of interest. This wasn't exactly a lot of data in order to predict across the entire study area, but luckily a [TEM project](http://a100.gov.bc.ca/pub/acat/public/viewReport.do?reportId=35895) was carried out in 2008 and there are publicly accessible points with additional observations of plant presence/absence and cover, as well as site series information and structural stage recordings as well. Unfortunately, the TEM polygons do not carry any structural stage attributes on them. Instead, that data is gathered by inferring from VRI datasets

The outputs of this script are as follows: 
1. Unattributed points for site series calls ("observations", used in the next script)
2. TEM polygons with only compositional data attached to it ("polygons", used in the next script)
3. Compositional site series data for each TEM polygon ID ("composition", used in the next script)
4. Unattributed points for structural stage calls ("observations", used in the next script)
5. VRI polygons with only compositional data attached to it ("polygons", used in the next script)
6. Compositional structural stage data for each VRI polygon ID ("composition", used in the next script)
7. Attributed points for plant presence/absence as well as cover for each plant species of interest (used two scripts ahead)

Numbers 1-3 and 4-6 above are discrete sets of data used in the DSMART modelling algorithm to produce rasters of site series and structural stage, which can then be used in presence/absence and cover modelling. Number 7 creates points that will be used in the presence/absence and cover modelling only.

### [02b_DSMART](https://github.com/mcoghill/Sechelt_Cultural_Plants/blob/master/02b_DSMART.Rmd)
Disaggregating and harmonising soil map units through resampled classification trees (DSMART): this is an algorithm that takes spatial polygons with classified compositional data and rasterizes them through a modelling process. The original method uses the ```raster``` package, which has been superceded by the ```terra``` package; however, this is not incorporated in the [```rdsmart```](https://bitbucket.org/brendo1001/dsmart/src/master/) package. Instead, I have gone through and rewritten the code to conform to the ```terra``` package standards ([custom DSMART code](https://github.com/mcoghill/Sechelt_Cultural_Plants/tree/master/_functions/dsmart_custom)). The output includes a classified raster as well as probability rasters for each classification, both of which can be used in presence/absence modelling.

### [02c_Models](https://github.com/mcoghill/Sechelt_Cultural_Plants/blob/master/02c_Models.Rmd)
Since there were many large rasters being used in the modelling process, I needed a process that would be efficient. I turned to the [```mlr3```](https://github.com/mlr-org/mlr3) package in order to do this, since it can also perform a spatial cross validation modelling process which is more accurate than standard cross validation techniques.

We had a set of objectives with the modelling process as well: Does the TEM add any value to the predictive mapping? Should TEM be incorporated as a classified raster or as probabilities? Do climate variables matter? In order to answer all of these questions, various models needed to be created.

### [03_Model_Agreement](https://github.com/mcoghill/Sechelt_Cultural_Plants/blob/master/03_Model_Agreement.Rmd)
After generating the ~28 different models for each plant species, we found that the models produced fairly similar looking results. The next step was to then combine them to see where models agreed that plants were present. Since the model outputs a raster with values between 0 and 1 representing the probability of occurrence at a given pixel, a cutoff value was applied to say that a pixel was either present or absent of a given plant. Two cutoffs were used: 0.5 and 0.75 (i.e.: if a pixel value was 0.8, the plant was present; if the value was 0.2, the plant was absent). This was done for each predicted map, and then the maps were summed in order to produce an overall "agreement" map, whose values represented "number of maps that agreed a plant was present at a 0.5/0.75 cutoff value". It's an odd metric but worked quite well.

## 04: Point Data Enhancement from 2020 field collections
The 04 series of scripts follows the same idea as the 02 series of scripts, just adding in data from the 2020 field season in order to attempt to improve on the predictive maps as well as add in some extra species to try and create models for. The sampling process was different as well and allowed for more data to be collected in a shorter period of time.

This repository is currently active and will be getting regular updates as more coding is developed. If you have any questions about this please send me an email or submit an issue within this repository. Thanks!
