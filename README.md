# Extract indices over agricultural plots

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4318779.svg)](https://doi.org/10.5281/zenodo.4318780)

Available languages: [Spanish](https://github.com/CristianICS/sentinel-iv-gee/blob/main/README_esp.md)

Functions written in JavaScript for the Google Earth Engine API that enable the extraction of vegetation indices from Sentinel-1 and Sentinel-2 satellite imagery over agricultural plots, specifically those undergoing crop rotation.

Note: The file c6s_S2_batch.ipynb can be used to correct S2 images in a loop via [Google Colaboratory](https://colab.research.google.com/?utm_source=scs-index).

## Requirements

- Google Earth Engine
- [Shapefile format layer with the cultivated plots](https://code.earthengine.google.com/?asset=users/iranzocristian/explotacion_blcht_buffer20m):
  - Number of columns equivalent to the number of years in the series (column name = year)
  - The values in the rows correspond to the production in each plot, or, if not available, a code indicating whether they were produced or not in that year.

## Vegetation Indices

Available in the Sentinel-1 function:

- CR (*)

And in the Sentinel-2 function:

- NDVI
- NDRE
- IRECI

Any index can be added by defining its equation within the code and calling it in the main function.

(*) The calculation of the CR is performed to achieve a result similar to studies where it has been applied:

Veloso, A., Mermoz, S., Bouvet, A., Le Toan, T., Planells, M., Dejoux, J.-F., & Ceschia, E. (2017). Understanding the temporal behavior of crops using Sentinel-1 and Sentinel-2-like data for agricultural applications. Remote Sensing of Environment, 199, 415-426. https://doi.org/10.1016/j.rse.2017.07.015

Vreugdenhil, M., Wagner, W., Bauer-Marschallinger, B., Pfeil, I., Teubner, I., Rüdiger, C., & Strauss, P. (2018). Sensitivity of Sentinel-1 Backscatter to Vegetation Dynamics: An Austrian Case Study. Remote Sensing, 10(9), 1396. https://doi.org/10.3390/rs10091396

Sonobe, R. (2019). Combining ASNARO-2 XSAR HH and Sentinel-1 C-SAR VH/VV Polarization Data for Improved Crop Mapping. Remote Sensing, 11(16), 1920. https://doi.org/10.3390/rs11161920

## Result

The output of the functions is a graph displaying the average vegetation index (VI) values over the cultivated plots, calculated for each image in the collection.

## Specifications of the Sentinel-2 Function

From June 2015 to March 28, 2017, S2 images in Google Earth Engine (GEE) are only available at treatment level 1C, i.e., without atmospheric correction.

In the Sentinel-2 function's code, the second section allows the integration of images from this period that have been automatically corrected using the 6S model and stored within the GEE cloud. Atmospheric correction is applied using the code created by [Sam Murphy](https://github.com/samsammurphy/gee-atmcorr-S2), modified to correct all the images in the loop. This step can be omitted depending on the objectives of the study.

The file `c6s_S2_batch.ipynb` allows the 6S correction to be applied in a loop from Google Colaboratory. Download the file and upload it to Colab via the `File > Upload Notebook` option.

Important: Registration with the EarthEngine API in Python is done using the ID of a project in Google Cloud. A Google account and a Cloud project are required, which can be created by accessing the [Google Cloud Console](https://console.cloud.google.com/welcome).

![Find the project ID in Goolge Cloud.](https://github.com/CristianICS/sentinel-iv-gee/assets/58115393/90e9975f-6173-4899-b0bb-8f9ce9fa09b7)

The 6S correction can also be performed in a console with the [Py6S](https://py6s.readthedocs.io/en/latest/) and [EarthEngine API](https://developers.google.com/earth-engine/guides/python_install-conda#windows) packages installed (a [conda](https://docs.conda.io/projects/miniconda/en/latest/) workflow is recommended, as the Py6S module can be downloaded with the compiled 6S model). The code is the same as that used in the Colab file, but authentication must be done in the console (`EarthEngine` module installed). Instead of the function `ee.Authenticate()`, the command `earthengine authenticate` should be entered in the console.
