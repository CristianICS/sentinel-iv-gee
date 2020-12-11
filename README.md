# sentinel-iv-gee

Funciones para las constelaciones Sentinel-1 y Sentinel-2 que permiten la extracción de índices de vegetación sobre explotaciones agrícolas sometidas a un régimen de rotación de cultivos.

## Requisitos
- Google Earth Engine
- Capa en formato *shapefile* con los recintos cultivados:
  1. Nº de columnas equivalente al número de años de la serie (nombre de columna = año)
  2. Los valores de las filas se corresponden con la producción en cada recinto, o en su defecto con un código que indique si han sido producidos o no en dicho año.

## Resultado
La salida de las funciones es un gráfico con los valores medios del IV sobre los recintos cultivados, calculado sobre cada imagen de la colección.

## Especificaciones de la función Sentinel-2
Desde junio de 2015 hasta el 28 de marzo de 2017 las imágenes S2 en GEE solo disponen del nivel de tratamiento 1C, i.e., sin corrección atmosférica. 

En el código de la función para Sentinel-2, el segundo apartado permite integrar imágenes de este período corregidas automáticamente con el modelo 6s y almacenadas dentro de la nube de GEE. La corrección atmosférica se aplica mediante el código creado por Sam Murphy (https://github.com/samsammurphy/gee-atmcorr-S2), modificado para corregir en bucle todas las imágenes utilizadas. 

Este paso puede omitirse dependiendo de los objetivos del trabajo.
