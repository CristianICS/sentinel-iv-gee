/*
===============================================================================
========               GET METRICS OVER CULTIVATED CROPS               ========
========                         SENTINEL 2                            ========
===============================================================================
The code extracts mean values of different VI over cultivated crops.

VI (Vegetation Indices);

- NDVI
- NDRE
- IRECI

----------------------------------- INPUTS ------------------------------------
1)
SHP with plots from one agricultural holding. It have to contain columns with
the production units of each plot, and the column names need to be the final
year of the season in which production values are obtained.

Example:

| ID | 2014 | 2015 | 2016 |
|----|------|------|------|
| 1  | 400  | 0    | 250  |
| 2  | 0    | 800  | 0    |
| 3  | 300  | 0    | 100  |

The above agricultural holding have three fields which are cultivated
alternatively over the years. The even years two ones are used (the fallow ones
has 0 production), and odd years only one is used.
    
2) Sentinel 2 image collections:

- "COPERNICUS/S2_HARMONIZED" (from 2015-06-27)
- "COPERNICUS/S2_SR_HARMONIZED" (from 2017-03-28)

Note: In order to increase the number of scenes, images between 2015-06-27 and
2017-03-28 (collection S2_HARMONIZED) have been atmospherically corrected with
6S and included inside S2_SR_HARMONIZED collection. The S2_HARMONIZED collection
only have 1C products (TOA reflectance).

--------------------------------- FUNCTIONS -----------------------------------
1) Define parameters and functions to perform a shadow mask, plus functions to
compute the VIs.

2) Initialize S2 collections.

3) Define main function (L359), which creates a FeatureCollection with VI mean
values per image over cultivated plots in the current season.

Note: Inside the example agricultural holding plots (barley), the season starts
in early October and ends in late July.

--------------------------------- REFERENCES -----------------------------------
Schmitt, M., Hughes, L.H., Qiu, C., Zhu, X.X., 2019. Aggregating cloud-free
Sentinel-2 images with Google Earth Engine, in: ISPRS Annals of Photogrammetry,
Remote Sensing and Spatial Information Sciences. Presented at the ISPRS ICWG II/III
PIA19+MRSS19 - Photogrammetric Image Analysis & Munich Remote Sensing 
Symposium: Joint ISPRS conference (Volume IV-2/W7) - Copernicus GmbH, pp. 145–152.
https://doi.org/10.5194/isprs-annals-IV-2-W7-145-2019
https://www.isprs-ann-photogramm-remote-sens-spatial-inf-sci.net/IV-2-W7/145/2019/

Wilson, R. T. (2013). Py6S: A Python interface to the 6S radiative transfer model.
Comput. Geosci. - UK, 51(2), 166-171.
----------------------------------------------------------------------------------
*/

var AOI = ee.FeatureCollection('users/iranzocristian/explotacion');
// Link to obtain the above layer
// https://code.earthengine.google.com/?asset=users/iranzocristian/explotacion

// Thresholds for Cloud and Shadow masks:
// 1. Cloud pixels are >= CLOUD_PROB
var CLOUD_PROB = 45; 
// 2. Shadow pixels are >= SHADOW_PROB
var SHADOW_PROB = 20; 
// 3. Water pixels are <= NDVI_THRESHOLD
var NDVI_THRESHOLD = -0.1;
// 4. Dark pixels are <= IR_THRESHOLD
var IR_THRESHOLD = 0.3; 

// Parameters for extending shadow mask over near pixels (Schmitt et al., 2019)
var ERODE_PIXELS = 1.5;
var DILATION_PIXELS = 3;
// Cloud heights to proyect shadow masks
var CLOUD_HEIGHTS = ee.List.sequence(200,10000,250);

/**
 * Compute cloud shadow mask
 * ===========================================================================
 * From (Schmitt et al., 2019)
 * https://www.isprs-ann-photogramm-remote-sens-spatial-inf-sci.net/IV-2-W7/145/2019/
 * @param {EEImage} img It must contain a CLOUD_PROBABILITY mask 
 * @returns 
 */
function projectShadows(img){
  
  // Get solar angles (perform illumination geometry calculations)
  var mean_azi = img.get('MEAN_SOLAR_AZIMUTH_ANGLE');
  var mean_zen = img.get('MEAN_SOLAR_ZENITH_ANGLE');
  
  // Get cloud probability mask
  var clouds = ee.Image(img.get('cloud_mask')).select('probability');
  var cloud_mask = clouds.gt(CLOUD_PROB);
    
  // Find dark pixels (infrared region)
  var dark_pixels = img.select(['B8','B11','B12']).reduce(ee.Reducer.sum());
  var dark_mask = dark_pixels.lt(IR_THRESHOLD);

  // Find dark pixels linked with water (not cloud shadows)
  var ndvi = img.normalizedDifference(['B8','B4']);
  var water_mask = ndvi.lt(NDVI_THRESHOLD);
  
  // Get dark pixels (probably) related with cloud shadows
  var shadow_mask = dark_mask.and(water_mask.not()); // Exclude water
  shadow_mask = shadow_mask.and(cloud_mask.not()); // Exclude clouds
     
  // Detect shadows by illumination geometry 
  // Angles in radians: multiply the number of degrees by pi/180
  var mean_azi_radians = ee.Number(mean_azi).add(180)
    .multiply(Math.PI).divide(180.0);
  var mean_zen_radians = ee.Number(mean_zen)
    .multiply(Math.PI).divide(180.0);
  // Note: Math.PI is a JavaScript Object
      
  /**
   * Find shadows
   * =========================================================================
   * Move cloud pixels to its likely shadow position.
   * 
   * @param {Integer} cloud_height 
   * @returns A list with EEImages (one per height value)
   */
  var findShadows = function(cloud_height){
    
    var height = ee.Number(cloud_height);
    
    // Compute projected shadow distance related with height
    var shadow_casted_dist = mean_zen_radians.tan().multiply(height);
    // Coordinate X shadow distance
    var x = mean_azi_radians.sin().multiply(shadow_casted_dist)
      .multiply(-1);
    // Coordinate Y shadow distance
    var y = mean_azi_radians.cos().multiply(shadow_casted_dist)
      .multiply(-1);
    
    // Get and image containing displacement values.
    var shadow_img = ee.Image.constant(x).addBands(ee.Image.constant(y));
    // Shift clouds values over their shadows
    return clouds.displace(shadow_img);
  };
  
  // Obtain expected shadow pixels from cloud values by height
  var shadows = CLOUD_HEIGHTS.map(findShadows);
    
  // Group the above shadow images inside a collection
  var shadow_coll = ee.ImageCollection.fromImages(shadows);
  // Get mean shadow probability over all cloud heights
  var shadow_height_mask = shadow_coll.mean();
    
  // Update shadow mask (it has been initialized with dark pixels)
  shadow_mask = shadow_height_mask.multiply(shadow_mask);
  
  /**
   * Apply circular kernel
   * =========================================================================
   * Smooth shadow mask.
   * 
   * @param {EEImage} shadow_prob EEImage with shadow probability
   * @returns 
   */
  function dilatedErossion(shadow_prob) {

    // Apply kernels
    var dilated = shadow_prob
    // .reproject('EPSG:4326', null, 20)
    .focal_min({
      radius: ERODE_PIXELS, kernelType: 'circle', iterations:3
    })
    .focal_max({
      radius: DILATION_PIXELS, kernelType: 'circle', iterations:3
    })
    // .reproject('EPSG:4326', null, 20);
                  
    return(dilated);
  }
  
  var shadow_mask_dilated = dilatedErossion(shadow_mask);
   
  // Perform a reduce neighbourhood operation
  var shadow_score = shadow_mask_dilated.reduceNeighborhood(
      {
          reducer: ee.Reducer.max(),
          kernel: ee.Kernel.square(1)
      });
  
  // Add shadow propability mask to the original image
  return img.addBands(shadow_score.rename(['MSK_SDWPRB']));
}

/**
 * Mask pixels with clouds
 * ===========================================================================
 * 
 * @param {EEImage} img It must contain 'cloud_mask' property. 
 * @returns 
 */
var applyCloudMask = function(img){
  // Get probability band from COPERNICUS/S2_CLOUD_PROBABILITY collection
  var cloud_prob = ee.Image(img.get('cloud_mask')).select('probability');
  // Create band with 1 values for pixels below cloud probability threshold  
  var cloud_mask = cloud_prob.lt(CLOUD_PROB);
  
  return img.updateMask(cloud_mask);
};

/**
 * Mask pixels with cloud shadows
 * ===========================================================================
 * 
 * @param {EEImage} img
 * @returns 
*/
var applyShadowMask = function(img){

  // Compute a cloud shadow probability mask inside the image
  var simg = projectShadows(img);  
  // Create band with 1 values for pixels below shadow probability threshold  
  var shadow_mask = simg.select('MSK_SDWPRB').lt(SHADOW_PROB);

  return img.updateMask(shadow_mask);
};

/**
 * Apply scale factor
 * ===========================================================================
 * Add a constant over stored S2 values to return reflectance plus REMAIN the
 * original image properties.
 * 
 * Note: This code could be executed with original INT values, because the VI
 * compute a ratio between several bands (whatever the value is). The S2 imgs
 * are corrected by the constant because the 6S correction yields reflectance
 * values.
 * 
 * The 'system:index' is written inside fileID property inside the 6S
 * corrected images. 'system:index' property in S2_SR collection is set inside
 * 'fileID' property in order to join both collections with
 * COPERNICUS/S2_CLOUD_PROBABILITY collection.
 * 
 * @param {EEImage} img 
 * @returns 
 */
var reflectance = function(img){
  var r = img.divide(10000);
  // Set new property
  r = r.set('fileID', img.id());
  // Add properties in the original image
  var props = img.propertyNames();
  
  return r.copyProperties(img, props);
};

/**
 * Compute Vegetation Indices
 * ==========================================================================
 * Add a new band per computed VI. Each VI band must contain a prefix 'VI_'
 * for selecting then inside getSeries() function.
 * 
 * Note: NDVI index is masked to remove saturated pixels and no vegetation.
 * 
 * @param {EEImage} img 
 * @returns 
 */
var computeVI = function(img){
  
  // Compute NDVI   
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('VI_NDVI');
  // Mask values with saturated pixels (gt 0.9) and water/soil (lt 0.1)
  var quality_mask = ndvi.lt(0.9).and(ndvi.gt(0.1));
  var ndvi_masked = ndvi.updateMask(quality_mask);
  
  // Compute NDRE
  var ndre = img.normalizedDifference(['B7', 'B4']).rename('VI_NDRE');
  // Compute IRECI
  var ireci = img.expression('(re3-r)/(re1/re2)',
    {
      re3: img.select('B7'),
      r: img.select('B4'),
      re1: img.select('B5'),
      re2: img.select('B6')
  }).rename('VI_IRECI');

  return img.addBands(ndvi_masked).addBands(ndre).addBands(ireci)
    .copyProperties(img, img.propertyNames());
};

// Load Sentinel 2 collection (2A)
var sentinel2 = ee.ImageCollection("COPERNICUS/S2_SR")
  .filterBounds(AOI)
  .map(reflectance);

// Load S2 images corrected with 6S
var assetList = ee.data.listAssets("users/iranzocristian/s2_6s")['assets']
  .map(function(asset) { return asset.name });
var sentinel2_6s = ee.ImageCollection(assetList);

// Merge the two above collections
sentinel2 = sentinel2.merge(sentinel2_6s);

// Load collection with Cloud Probability (same ID as S2 original images)
var s2clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(AOI);

// Join S2 with cloud probability dataset to add cloud mask.
// Note: cloud_mask is added to the result as an additional property
var sentinel2_clouds = ee.Join.saveFirst('cloud_mask').apply({
  primary: sentinel2,
  secondary: s2clouds,
  condition: ee.Filter.equals({
    leftField: 'fileID', rightField: 'system:index'})
});

// Mask clouds and shadows and compute VI
var sentinel2_mask = ee.ImageCollection(sentinel2_clouds)
  .map(applyCloudMask).map(applyShadowMask).map(computeVI);

// View CLOUD_PROB and SHADOW_PROB values in action
// Map.addLayer(ee.ImageCollection(sentinel2_clouds).first(), {bands:['B4','B3','B2'],min:0,max:0.3}, 'Original');
// Map.addLayer(sentinel2_mask.first(), {bands:['B4','B3','B2'],min:0,max:0.3}, 'Masked');

/**
 * Get VI mean values over predefined crop seasons
 * ===========================================================================
 * The function iterate over season years and compute VI mean values over the
 * images inside an image collection. Each VI band is reduced with cultivated
 * crops in the image date crop season.
 * 
 * E.g.: A crop season from 2015-10-01 to 2016-07-31. One image captured in
 * 2015-04-13 will be reduced using plots with more than 0 production units in
 * 2016 campaign (2016 column inside the SHP with the plots).
 * 
 * Chart errors:
 * 
 * The function execution could be too low because of reduceRegion function.
 * If this happens, it is possible to select only one VI band by 'vifilter'
 * parameter. Write the name of the band which contains the VI, e.g.,
 * VI_NDVI. These are calculated inside computeVI function.
 * 
 * If the chart wont be printed, added only one year inside 'season_years'
 * parameter like this: 
 * 
 * getSeries(sentinel2_mask, ee.List([2016]), start, end);
 * 
 * The export function should take 4/5 minutes to download data for the
 * three VIs and 4 crop seasons.
 *
 * @param {EEImageCollection} col It must contain images with VI bands. 
 * @param {EEList} season_years Season start years
 * @param {EEList} start MONTH and DAY numbers of the season start
 * @param {EEList} end MONTH and DAY nummbers of the season end
 * @param {String} vifilter String with vi band to select
 * @returns FeatureCollection with CR mean value and image date.
 */
var getSeries = function(col, season_years, start, end, vifilter){
  
  // Default parameter
  vifilter = vifilter || 'VI_[A-Za-z]+';
  
  // Iterate over season start years, returning a list.
  var mean_vi = season_years.map(function(year){
    
    // Get season periods
    var end_year = ee.Number(year).add(1);
    var season_start = ee.Date.fromYMD(year, start.get(0), start.get(1));
    var season_end = ee.Date.fromYMD(end_year, end.get(0), end.get(1));
    
    // Select images inside current season
    var filtered_col = col.filterDate(season_start, season_end);
    
    // Filter cultivated fields inside the season last year
    var str_year = season_end.get('year').format('%.0f'); // String mandatory
    var cultivated_plots = AOI.filter(ee.Filter.gt(str_year, 0));
    
    /**
     * Reduce VI over season cultivated plots
     * =========================================
     * Get bands with VI and obtain the mean value over the pixels
     * masked inside the cultivated plots in a current season.
     * 
     * Note: It's mandatory having a FeatureCollection called
     * "cultivated_plots".
     * 
     * @param {EEImage} img Image with VI bands.  
     * @returns Feature with an image date and CR mean value.
     */
    var reduceVI = function(img){
      // Select the VI bands (by regexp)
      // var vi = img.select('VI_[A-Za-z]+');
      var vi = img.select(vifilter);

      // Reduce VI values over cultivated plots (it returns a dict)
      var vi_mean = vi.reduceRegion(ee.Reducer.mean(), cultivated_plots);
      // Note: The scale is computed automatically to avoid time consuming
      
      // Add image date
      vi_mean = vi_mean.set('date', img.date());
      
      // Construct a Feature with the image date and the CR mean value
      return ee.Feature(null, vi_mean);
      
    };

    // Apply above function over all images inside filtered collection
    var fc = filtered_col.map(reduceVI);
    // The above line returns a FeatureCollection with equal number of
    // features than filtered images

    // Important: A map() method applied over an ee.List must return a list.
    // In order to return a list, each feature inside prior FeatureCollection
    // is passed inside a list: [<Feature>,<Feature>,...]
    return fc.toList(100); // The max entities number is mandatory.
    
  });
  
  // The object mean_cr is a list of lists with Features
  // [[<Feature>,<Feature>,...],[],...]

  // Obtain a FeatureCollection with all Features
  // 1) Flatten the list (get only one list with all features)
  // 2) Transform its to a FeatureCollection
  return ee.FeatureCollection(mean_vi.flatten());
};

// Apply getSeries with custom season parameters
// Note: Inside the example agricultural holding plots (barley), 
// the season starts in early October and ends in late July.

// Season start years
var years = ee.List.sequence(2015, 2019);
// Season start: month and day
var start = ee.List([10,1]);
// Season end: month and day
var end = ee.List([8,1]);
// End is exclusive in ee.ImageCollection.filterDate(start,end)

// Get FeatureCollection with all images reduced by cultivated fields
var vi_mean = getSeries(sentinel2_mask, years, start, end).sort('date');

print('Number of computed images:', vi_mean.size());

// Plot values in a chart
var chart = ui.Chart.feature.byFeature(vi_mean, 'date', 'VI_NDVI');
print(chart);

// Export FeatureCollection to CSV
Export.table.toDrive({
  collection: vi_mean,
  description: 'Export_vi_mean_over_cultivated_plots',
  folder: 'data',
  fileNamePrefix: 'vi_mean'
});
// Note: system:index property is the image ID from the VI mean is computed.
