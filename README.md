# sentinel-iv-gee
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4318779.svg)](https://doi.org/10.5281/zenodo.4318780)

Funciones para las constelaciones Sentinel-1 y Sentinel-2 que permiten la extracción de índices de vegetación sobre explotaciones agrícolas sometidas a un régimen de rotación de cultivos.

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

En el código de la función para Sentinel-2, el segundo apartado permite integrar imágenes de este período corregidas automáticamente con el modelo 6s y almacenadas dentro de la nube de GEE. La corrección atmosférica se aplica mediante el código creado por Sam Murphy (https://github.com/samsammurphy/gee-atmcorr-S2), modificado para corregir en bucle todas las imágenes utilizadas. Este paso puede omitirse dependiendo de los objetivos del trabajo.

### Corrección atmosférica en bucle
A continuación se detalla la modificación del código de Sam Murphy realizada para corregir en bucle todas las imágenes de una colección. Esta se incluye dentro del archivo jupyter notebook que se ejecuta al iniciar el repositorio: `sentinel2_atmospheric_correction.ipynb`.

1) Se incluyen los paquetes necesarios y se inicia GEE.
```py
# Incluir paquetes necesarios
import ee # Modulo de GEE
from Py6S import *
import datetime
import math
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.getcwd()),'bin'))
from atmospheric import Atmospheric
# Inicio de Earth Engine con la API de Python
ee.Initialize()
```
2) Se definen las tres funciones creadas por Sam Murphy
```py
# 1- Extraer la respuesta espectral
def spectralResponseFunction(bandname):

  """
  Respuesta espectral
  de las bandas Sentinel 2
  """
  bandSelect = {
  'B1':PredefinedWavelengths.S2A_MSI_01,
  'B2':PredefinedWavelengths.S2A_MSI_02,
  'B3':PredefinedWavelengths.S2A_MSI_03,
  'B4':PredefinedWavelengths.S2A_MSI_04,
  'B5':PredefinedWavelengths.S2A_MSI_05,
  'B6':PredefinedWavelengths.S2A_MSI_06,
  'B7':PredefinedWavelengths.S2A_MSI_07,
  'B8':PredefinedWavelengths.S2A_MSI_08,
  'B8A':PredefinedWavelengths.S2A_MSI_09,
  'B9':PredefinedWavelengths.S2A_MSI_10,
  'B10':PredefinedWavelengths.S2A_MSI_11,
  'B11':PredefinedWavelengths.S2A_MSI_12,
  'B12':PredefinedWavelengths.S2A_MSI_13,
  }
  return Wavelength(bandSelect[bandname])
  
# 2- Conversion de reflectancia a radiancia3
def toa_to_rad(bandname):
  """
  Reflectividad TOA a radiancia
  """
  # Calcular la irradiancia solar exoatmosferica
  ESUN = info['SOLAR_IRRADIANCE_'+bandname]
  solar_angle_correction = math.cos(math.radians(solar_z))
  # Distancia Tierra-Sol (doy)
  doy = scene_date.timetuple().tm_yday
  d = 1 - 0.01672 * math.cos(0.9856 * (doy-4))
  # http://physics.stackexchange.com/
  # questions/177949/earth-sun-distance-on-a-given-day-of-the-year
  # factor de conversion
  multiplier = ESUN*solar_angle_correction/(math.pi*d**2)
  # calculo de radiancia
  rad = toa.select(bandname).multiply(multiplier)
  return rad
 
# 3- Calculo de reflectividad BOA
def surface_reflectance(bandname):
  """
  Calculo de la reflectividad en superficie a traves de la
  radiancia del sensor (en funcion de la longitud de onda,
  especifica para cada banda)
  """
  # Extraer la respuesta espectral de la banda
  s.wavelength = spectralResponseFunction(bandname)
  # Ejecutar los objetos 6s (definidos en la funcion principal)
  s.run()
  # Extraer las incognitas atmosfericas
  Edir = s.outputs.direct_solar_irradiance # irradiancia solar directa
  Edif = s.outputs.diffuse_solar_irradiance # irradiancia solar difusa
  Lp = s.outputs.atmospheric_intrinsic_radiance # path radiance
  absorb = s.outputs.trans['global_gas'].upward # absorption transmissivity
  scatter = s.outputs.trans['total_scattering']\
  .upward # scattering transmissivity
  tau2 = absorb*scatter # total transmissivity
  
  # Nota: los s.outputs son calculados automaticamente a partir de los
  # objetos 6s definidos en la funcion de conversion principal, "conversion".
  # Transformar los valores de reflectividad TOA a radiancia
  rad = toa_to_rad(bandname)
  # despejar la ecuacion de transferencia radiativa
  ref = rad.subtract(Lp).multiply(math.pi).divide(tau2*(Edir+Edif))
  # Devuelve la reflectividad a BOA de una banda
  return ref
```
3) Definir el AOI sobre el que corregir la imagen
```py
# Incluir el area de estudio:
# Poligono a partir del cual se filtra la coleccion,
# se recorta la imagen final y se calculan los parametros
# atmosfericos necesarios
geom = ee.Geometry.Polygon([[-0.9570796519011493,40.98197275411647],
[-0.5670650034636493,40.98197275411647],
[-0.5670650034636493,41.45919658393617],
[-0.9570796519011493,41.45919658393617],
[-0.9570796519011493,40.98197275411647]])

# Descomentar la siguiente linea si ee.Geometry.Polygon no funciona
# geom = ee.Geometry.Rectangle(-0.996, 41.508, -0.568, 40.992)
# Obtener las coordenadas de geom para recortar las imagenes de la coleccion
region = geom.buffer(1000).bounds().getInfo()['coordinates']
```
4) Creacion de la función para corregir imágenes en bucle
```py
def conversion(img):
  # Incorporar la fecha de la imagen
  date = img.date()
  # Definir las variables globales:
  # Aquellas que pueden ser llamadas fuera del entorno de la funcion.
  global toa
  global info
  global scene_date
  global solar_z
  
  # calcular la reflectividad a TOA
  toa = img.divide(10000)
  
  # Escribir los metadatos de la imagen
  # Recopilar las propiedades
  info = img.getInfo()['properties']
  # Fecha: Python utiliza segundos, EE milisegundos
  scene_date = datetime.datetime\
  .utcfromtimestamp(info['system:time_start']/1000)
  # Angulo cenital solar
  solar_z = info['MEAN_SOLAR_ZENITH_ANGLE']
  # Valores sobre la composicion atmosferica
  # El codigo de las funciones se encuentra dentro
  # del repositorio de samsammurphy (gee-atmcorr-S2/bin/atmospheric.py)
  h2o = Atmospheric.water(geom,date).getInfo()
  o3 = Atmospheric.ozone(geom,date).getInfo()
  # Atmospheric Optical Thickness
  aot = Atmospheric.aerosol(geom,date).getInfo()
  # Altura de la superficie, a partir del MDE de la mision
  # Shuttle Radar Topography mission (STRM) en GEE
  SRTM = ee.Image('CGIAR/SRTM90_V4')
  
  # Calculo de la altura media del ´area de estudio (geom)
  alt = SRTM.reduceRegion(reducer = ee.Reducer.mean(),
  geometry = geom.centroid()).get('elevation').getInfo()
  # Transformar a km, medida utilizada por Py6s
  km = alt/1000
  
  """
  Inicio de los objetos 6s, columna vertebral de Py6s
  A partir de la clase 6s se definen los parametros
  requeridos por la funcion de transferencia radiativa
  Llamar a los objetos 6s
  """
  global s
  s = SixS()
  # Integrar los componentes atmosfericos
  s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o,o3)
  s.aero_profile = AeroProfile.Continental
  s.aot550 = aot
  # Calcular la geometria Earth-Sun-satellite
  s.geometry = Geometry.User()
  s.geometry.view_z = 0 # calculo asumiendo vision en NADIR
  s.geometry.solar_z = solar_z # angulo cenital solar
  s.geometry.month = scene_date.month # mes usado en la distancia Earth-Sun
  s.geometry.day = scene_date.day # dia usado en la distancia Earth-Sun
  s.altitudes\
  .set_sensor_satellite_level() # Altitud del sensor
  s.altitudes\
  .set_target_custom_altitude(km) # Altitud de la superficie
  
  # Aplicar la conversion a cada banda de la imagen
  # 1. Generar el objeto (imagen) a exportar
  output = img.select('QA60')
  # 2. Bucle de correccion: aplica la funcion de correccion a las bandas
  # de la lista
  for band in ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12']:
    print(band)
    # Corregir la banda e incluirla en la imagen a exportar
    output = output.addBands(surface_reflectance(band))
    
  # Exportar la imagen a una carpeta en GEE
  # 1. Definir parametros de la imagen a exportar
  dateString = scene_date.strftime("%Y-%m-%d")
  ref = output.set({'satellite':'Sentinel 2',
    'fileID':info['system:index'],
    'date':dateString,
    'aerosol_optical_thickness':aot,
    'water_vapour':h2o,
    'ozone':o3})
  # 2. Definir la carpeta de destino (dentro de GEE)
  assetID = 'users/iranzocristian/6s_test/S2SR_'+dateString
  # 3. Opciones de la imagen a exportar
  export = ee.batch.Export.image.toAsset(\
    image=ref,
    description='sentinel2_atmcorr_export',
    assetId = assetID,
    region = region,
    crs = 'EPSG:4326',
    scale = 20)
  # 4. Exportar la imagen
  export.start()
  return print("imagen "+assetID+" exportada")
# Final de la funcion de conversion en bucle
```
5) Aplicar la funcion de conversion en bucle a una coleccion GEE. Primero se define la colección y despues se aplica la conversión dentro de un bucle a todas las imágenes que contiene.
```py
# Definir coleccion GEE
S2 = ee.ImageCollection('COPERNICUS/S2')\
  .filterBounds(geom)\
  .filterDate('2015-10-01','2017-04-30')\
  .filterMetadata('MGRS_TILE', 'equals', '30TXL')\
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 20)\
  .sort('system:time_start')\
  .distinct('system:time_start')
  
# Definir en una lista las imagenes a filtrar
features = S2.getInfo()['features']

"""
CORRECCION DE LA COLECCION AUTOMATICAMENTE (bucle for)
1. Recorre cada imagen de la lista anterior
2. Obtiene su id
3. Se llama a la imagen de la coleccion GEE con el id anterior
4. Se aplica la funcion de conversion principal
"""
for i in features:
  id = i['id']
  conversion(ee.Image(id))
# Final del Script
```
