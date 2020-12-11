# sentinel-iv-gee

Funciones para las constelaciones Sentinel-1 y Sentinel-2 que permiten la extracción de índices de vegetación sobre explotaciones agrícolas sometidas a un régimen de rotación de cultivos.

## Requisitos
- Google Earth Engine
- Capa en formato *shapefile* con los recintos cultivados:
  1. Nº de columnas equivalente al número de años de la serie (nombre de columna = año)
  2. Los valores de las filas se corresponden con la producción en cada recinto, o en su defecto con un código que indique si han sido producidos o no en dicho año.
