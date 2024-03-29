# sentinel-iv-gee

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4318779.svg)](https://doi.org/10.5281/zenodo.4318780)

Funciones para las constelaciones Sentinel-1 y Sentinel-2 que permiten la extracción de índices de vegetación sobre explotaciones agrícolas sometidas a un régimen de rotación de cultivos.

Se añade el fichero `c6s_S2_batch.ipynb` con el que se pueden corregir en bucle imágenes S2 a través de [Google Colaboratory](https://colab.research.google.com/?utm_source=scs-index).

## Requisitos
- Google Earth Engine
- Capa en formato *shapefile* con los recintos cultivados:
  1. Nº de columnas equivalente al número de años de la serie (nombre de columna = año)
  // Enlace para obtener la capa con la que se desarrolla el SCRIPT https://code.earthengine.google.com/?asset=users/iranzocristian/explotacion_blcht_buffer20m
  2. Los valores de las filas se corresponden con la producción en cada recinto, o en su defecto con un código que indique si han sido producidos o no en dicho año.
  
## Índices de vegetación
Están disponibles, en la función de Sentinel-1
- CR (\*)

y en la función Sentinel-2:
- NDVI
- NDRE
- IRECI

Puede añadirse cualquier índice definiendo su ecuación dentro del código y llamándola en la función principal.

(\*) El cálculo del CR se realiza para obtener un resultado similar a los trabajos donde ha sido aplicado:

Veloso, A., Mermoz, S., Bouvet, A., Le Toan, T., Planells, M., Dejoux, J.-F., & Ceschia, E. (2017). Understanding the temporal behavior of crops using Sentinel-1 and Sentinel-2-like data for agricultural applications. Remote Sensing of Environment, 199, 415-426. https://doi.org/10.1016/j.rse.2017.07.015

Vreugdenhil, M., Wagner, W., Bauer-Marschallinger, B., Pfeil, I., Teubner, I., Rüdiger, C., & Strauss, P. (2018). Sensitivity of Sentinel-1 Backscatter to Vegetation Dynamics: An Austrian Case Study. Remote Sensing, 10(9), 1396. https://doi.org/10.3390/rs10091396

Sonobe, R. (2019). Combining ASNARO-2 XSAR HH and Sentinel-1 C-SAR VH/VV Polarization Data for Improved Crop Mapping. Remote Sensing, 11(16), 1920. https://doi.org/10.3390/rs11161920


## Resultado
La salida de las funciones es un gráfico con los valores medios del IV sobre los recintos cultivados, calculado sobre cada imagen de la colección.

## Especificaciones de la función Sentinel-2
Desde junio de 2015 hasta el 28 de marzo de 2017 las imágenes S2 en GEE solo disponen del nivel de tratamiento 1C, i.e., sin corrección atmosférica. 

En el código de la función para Sentinel-2, el segundo apartado permite integrar imágenes de este período corregidas automáticamente con el modelo 6S y almacenadas dentro de la nube de GEE. La corrección atmosférica se aplica mediante el código creado por [Sam Murphy](https://github.com/samsammurphy/gee-atmcorr-S2), modificado para corregir en bucle todas las imágenes utilizadas. Este paso puede omitirse dependiendo de los objetivos del trabajo.

El fichero `c6s_S2_batch.ipynb` permite realizar la corrección 6S en bucle desde Google Colaboratory. Descargar el fichero y subirlo a Colab mediante la opción `Archivo > Subir cuaderno`.


**Importante:** El registro en la API de EarthEngine en Python se realiza mediante el ID de un proyecto en Google Cloud. Es necesario tener una cuenta en Google y un proyecto en Cloud, el cual puede crearse accediendo a la [consola de Google Cloud](https://console.cloud.google.com/welcome).

![Localización del ID del proyecto en Goolge Cloud.](https://github.com/CristianICS/sentinel-iv-gee/assets/58115393/90e9975f-6173-4899-b0bb-8f9ce9fa09b7)

La corrección 6S también puede realizarse en una consola con los paquetes [Py6S](https://py6s.readthedocs.io/en/latest/) y [EarthEngine API](https://developers.google.com/earth-engine/guides/python_install-conda#windows) instalados (recomendado el flujo de trabajo en [conda](https://docs.conda.io/projects/miniconda/en/latest/), pues el módulo Py6S puede descargarse con el modelo 6S compilado). El código es igual al aplicado en el archivo de Colab, pero el registro se debe realizar en la consola (módulo `EarthEngine` instalado). En lugar de la función `ee.Authenticate()`, debe introducirse en la consola el comando `earthengine authenticate`.

