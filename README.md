# sentinel-iv-gee

Funciones para las constelaciones Sentinel-1 y Sentinel-2 que permiten la extracción de índices de vegetación sobre explotaciones agrícolas sometidas a un régimen de rotación de cultivos.

## Requisitos
- Google Earth Engine
- Capa en formato *shapefile* con los recintos cultivados:
  1. Nº de columnas equivalente al número de años de la serie (nombre de columna = año)
  2. Los valores de las filas se corresponden con la producción en cada recinto, o en su defecto con un código que indique si han sido producidos o no en dicho año.

## Función
`S1serietemporal(coleccion, a_inicio, mes_inicio, dia_inicio, a_final, mes_final, dia_final)`

donde
- `coleccion`: Colección Sentinel 1
- `a_inicio`: Año de inicio de la primera temporada de cultivo a analizar
- `mes_inicio`: Mes de inicio de la temporada regular
- `dia_inicio`: Día de inicio "
- `a_final`: Año final de la serie a analizar 
- `mes_final`: Mes final de la temporada regular
- `dia_final`: Día final "

## Desarrollo de la función



**AVISO**: Si la serie de años es muy larga es posible que haya que ejecutar la función varias veces para mostrar el gráfico final.
