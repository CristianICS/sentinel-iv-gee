/*
===============================================================================
========               GET METRICS OVER CULTIVATED CROPS               ========
========                         SENTINEL 1                            ========
===============================================================================

The code extracts mean values of VH/VV (CR) index over cultivated crops.

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
    
2) Sentinel 1 image collection ("COPERNICUS/S1_GRD").
    
--------------------------------- FUNCTIONS -----------------------------------
1) Perform CR index. First, dB values are transformed to sigma units. Then,
CR ratio is computed and, finally, sigma CR values are converted to dB units.

2) Define main function, which creates a FeatureCollection with CR mean values
per image over cultivated plots in the current season.

Note: Inside the example agricultural holding plots (barley), the season starts
in early October and ends in late July.

-------------------------------------------------------------------------------
*/

// Load SHP with the agricultural holding plots
var AOI = ee.FeatureCollection('users/iranzocristian/explotacion');

// Link to obtain the above layer
// https://code.earthengine.google.com/?asset=users/iranzocristian/explotacion

// Load Sentinel1 collection
var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  // Filter images with VV or VH polarisation
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  // 'Interferometric Wide Swath'
  .filter(ee.Filter.eq('instrumentMode', 'IW'))             
  // High resolution images
  .filter(ee.Filter.eq('resolution', 'H'))
  .filterBounds(AOI);


// Uncomment the next code chuck to compare the two orbits indicent angles
// ---------------------------------- Init chunk
// var s1_des = sentinel1
//   .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
//   .select('angle');
// var s1_asc = sentinel1
//   .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
//   .select('angle');
// var des_chart = ui.Chart.image.series(s1_des, AOI, ee.Reducer.mean(), 2000, 'system:time_start');
// print(des_chart.setOptions({'title':'Mean Angle, Descending orbit'}));
// var asc_chart = ui.Chart.image.series(s1_asc, AOI, ee.Reducer.mean(), 2000, 'system:time_start');
// print(asc_chart.setOptions({'title':'Mean Angle, Ascending orbit'}));
// ---------------------------------- Finish chunk

// Select descending orbits, less incident angle oscilations (see above chunk)
var sentinel1_des = sentinel1.filter(
    ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

// Define functions to apply CR ratio

/**
 * Apply inverse log transformation
 * ============================================================================
 * 
 * From dB to lineal sigma units.
 * 
 * @param {EEImage} img
 * @returns 
*/
var invLog = function(img){
    
    var expression = 'pow(10,(b/10))';
    var vh_sg = img.expression(expression, {b: img.select('VH')}).toDouble();
    var vv_sg = img.expression(expression, {b: img.select('VV')}).toDouble();
    
    var transformed_img = vh_sg.select([0], ['VH_sigma'])
    .addBands(vv_sg.select([0], ['VV_sigma']))

    return img.addBands(transformed_img);
};

/**
 * Apply log transformation
 * ============================================================================
 * Convert band values from lineal sigma units to dB.
 * 
 * @param {EEImage} img 
 * @param {String} bandname CR 
 * @returns 
 */
var toDB = function(img, bandname){

    var transformed = img.expression('10 * log10(b)', 
    {b: img.select(bandname)}).toDouble();
    
    // Update band with transformed one
    return transformed;
};

/**
 * Compute VH/VV coeficient (Cross Ratio)
 * ============================================================================
 * 
 * Compute CR index. 
 * 
 * The VH and VV bands are transformed to sigma values prior performing the
 * calculation. Then, the CR are converted again into db units.
 * 
 * @param {EEImage} img 
 * @returns 
*/
var applyCR = function(img){
  
  // First perform the inverse log transformation (obtain sigma)
  var sigma = invLog(img);
  
  // Then, apply CR index to above bands
  var cr = sigma.expression('vh / vv',
  {
      vh: sigma.select('VH_sigma'),    
      vv: sigma.select('VV_sigma'),    
  }).toDouble()
  // Rename the band with the index
  .select([0], ['CR']);
  
  // Transform CR values (in sigma) to dB
  var cr_db = toDB(cr, 'CR');
  
  // Add index to the original image
  return img.addBands(cr_db.select([0], ['CR']));
};

// Define the function to perform the CR series.

/**
 * Get CR mean values over predefined crop seasons
 * ===========================================================================
 * The function iterate over season years and compute CR mean values over the
 * images inside an image collection. Each CR band is reducing with cultivated
 * crops in the image date crop season.
 * 
 * E.g.: A crop season from 2015-10-01 to 2016-07-31. One image captured in
 * 2015-04-13 will be reduced using plots with more than 0 production units in
 * 2016 campaign (2016 column inside the SHP with the plots).
 * 
 * @param {EEImageCollection} col It must contain images with CR band. 
 * @param {EEList} season_years Season start years
 * @param {EEList} start MONTH and DAY numbers of the season start
 * @param {EEList} end MONTH and DAY nummbers of the season end
 * @returns FeatureCollection with CR mean value and image date.
 */
var getSeries = function(col, season_years, start, end){
  
  // Iterate over season start years, returning a list.
  var mean_cr = season_years.map(function(year){
    
    // Get season periods
    var end_year = ee.Number(year).add(1);
    var season_start = ee.Date.fromYMD(year, start.get(0), start.get(1));
    var season_end = ee.Date.fromYMD(end_year, end.get(0), end.get(1));
    
    // Select images inside current season
    var filtered_col = col.filterDate(season_start, season_end);
    
    // Compute CR from above filtered collection
    var cr_db = filtered_col.map(applyCR);
    
    // Filter cultivated fields inside the season last year
    var str_year = season_end.get('year').format('%.0f'); // String mandatory
    var cultivated_plots = AOI.filter(ee.Filter.gt(str_year, 0));
    
    /**
     * Reduce CR over season cultivated plots
     * =========================================
     * Get an image with CR band and obtain the mean value over the pixels
     * masked inside the cultivated plots in a current season.
     * 
     * Note: It's mandatory having a FeatureCollection called
     * "cultivated_plots".
     * 
     * @param {EEImage} img Image with CR band.  
     * @returns Feature with an image date and CR mean value.
     */
    var reduceCR = function(img){
      // Select the CR band
      var cr = img.select('CR');

      // Reduce CR values over cultivated plots
      var cr_mean = cr.reduceRegion(ee.Reducer.mean(), cultivated_plots);
      // Note: The scale is computed automatically to avoid time consuming

      // Construct a Feature with the image date and the CR mean value
      return ee.Feature(null, {'date': img.date(), 'CR': cr_mean.get('CR')});
      
    };

    // Apply above function over all images inside filtered collection
    var fc = cr_db.map(reduceCR);
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
  return ee.FeatureCollection(mean_cr.flatten());
};

// Apply getSeries with custom season parameters
// Note: Inside the example agricultural holding plots (barley), 
// the season starts in early October and ends in late July.

// Season start years
var years = ee.List.sequence(2014, 2019);
// Season start: month and day
var start = ee.List([10,1]);
// Season end: month and day
var end = ee.List([8,1]);
// End is exclusive in ee.ImageCollection.filterDate(start,end)

// Get FeatureCollection with all images reduced by cultivated fields
var cr_mean = getSeries(sentinel1_des, years, start, end).sort('date');

print('Number of computed images:', cr_mean.size());

// Plot values in a chart
var chart = ui.Chart.feature.byFeature(cr_mean, 'date', 'CR');
print(chart);

// Export FeatureCollection to CSV
Export.table.toDrive({
  collection: cr_mean,
  description: 'Export_CR_mean_over_cultivated_plots',
  folder: 'data',
  fileNamePrefix: 'cr_mean'
});
// Note: system:index property is the image ID from the CR mean is computed.
