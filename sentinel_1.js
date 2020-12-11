/*
=== OBTENCION DE MÉTRICAS SOBRE EL ESTADO DE UN CULTIVO ===
                    SENTINEL 1

El código extrae, en una función aplicada a la colección 
Sentinel 1, los valores medios del cociente VH/VV sobre los 
recintos cultivados. El resultado es un gráfico con todas las
imágenes en las temporadas comprendidas dentro de la serie 
temporal analizada.
  
Las imágenes se recortan sobre los recintos cultivados 
en la temporada durante la que han sido tomadas.

AVISO: Si la serie de años es muy larga es posible que haya 
que ejecutar la función varias veces para mostrar el gráfico
final

------------------|| INPUTS ||----------------------------
Capa con los recintos cultivados:
  - Nº de columnas equivalente al número de años de la serie 
  (nombre de columna = año)
  - Los valores de las filas se corresponden con la producción
  en cada recinto.

Colección Sentinel 1, polarizaciones VH y VV.

------------------|| FUNCIÓN ||---------------------------

S1serietemporal(coleccion, a_inicio,mes_inicio,dia_inicio,
                a_final,mes_final,dia_final);

donde "coleccion": Colección Sentinel 1
"a_inicio": Año de inicio de la primera temporada de
cultivo a analizar
"mes_inicio": Mes de inicio de la temporada regular
"dia_inicio": Día de inicio "
"a_final": Año final de serie a analizar
"mes_final": Mes final de la temporada regular
"dia_final": Día  final "

Nota: En Belchite la cebada se siembra a comienzos de octubre 
y se termina de cosechar a finales de julio.

----------------|| DESARROLLO DE LA FUNCIÓN ||---------------
Se declaran las funciones necesarias para calcular el 
cociente VH/VV:
1) Cálculo de la inversa del logaritmo, obteniendo nuevamente
valores sigma.
2) Aplicar el cociente entre las bandas sin transformar.
3) Volver a convertir el resultado en dB.

Este proceso asegura resultados comparables entre los trabajos
que aplican este índice.

Comienza la función principal:
- Se definen los años inicio y final sobre los que calcular el cociente. 
  El índice se calculará automáticamente entre las temporadas de cultivo 
  presentes dentro de dicho intervalo temporal. 
- Se define una variable que contendrá las colecciones de imágenes con
 las métricas calculadas de todas las temporadas. 
- Se crea la función que recorta las imágenes sobre los recintos cultivados 
  en la temporada de adquisición de la propia imagen.

A continuación se ejecuta un bucle FOR, filtrando la colección de
imágenes por cada temporada de cultivo, recortando las imágenes por 
los cultivos sembrados en la temporada y calculando el valor medio en
cada una de las imágenes. Los pasos se detallan a continuación (se
repiten durante todas las temporadas de la serie)
0) Se define el iterador como el primer año de la serie. A cada iteración 
   completada se sumará un año, y el bucle se interrumpirá al llegar al 
   último año de la serie.
1) Filtrar la capa con los cultivos de la explotación por el año y los
   valores de producción superiores a 0 (son los campos cultivados)
2) Seleccionar las imágenes dentro de la temporada de cultivo
3) Recortarlas por la capa espacial del paso 1)
4) Aplicar las funciones para crear el cociente VH/VV
5) Condicional IF: se vuelve al primer año de la temporada de cultivo del
   bucle y
   a) Si el año coincide con el primer año de la serie la colección de
   imágenes se integra a la colección creada en la primera parte de la función
   b) Si es un año distinto, la colección se suma a las colecciones integradas
   con anterioridad utilizando la función merge()

Al finalizar el bucle se muestra en un gráfico la serie temporal con todos los 
valores medios de CR durante las temporadas de cultivo analizadas. Se exporta a CSV 
manualmente a través de la consola de GEE.
---------------------------------------------------------------
*/

// ----------- CARGAR LA CAPA CON LOS RECINTOS ------------- \\
var AOI = ee.FeatureCollection('users/iranzocristian/explotacion_blcht_buffer20m');


// ----------------   FILTRAR COLECCIÓN   ------------------ \\

// Cargar la colección Sentinel-1
var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  // Obtener imágenes con polarizacion VV y VH
  .filter(ee.Filter.listContains(
    'transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains(
    'transmitterReceiverPolarisation', 'VH'))
                
  // Filtrar por las imágenes obtenidas en modo 
  // 'Interferometric Wide Swath'
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
                
  // Filtrar imágenes de alta resolucion
  .filter(ee.Filter.eq('resolution', 'H'));

// Dividir la coleccion en función de la pasada

// Se selecciona la trayectoria descendente, presenta mejores ángulos 
// de incidencia para analizar cultivos
var s1Descending = sentinel1.filter(
  ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

// ············| Cálculo del cociente VH/VV |·············\\

// Funciones para aplicar la transformación de la inversa 
// del logaritmo (de dB a valores sigma lineales)
  
// Transformación de la banda VH
  
var s1vh_invlog = function(img){
    
  // Expresión
  var vh_log = img.expression(
  'pow(10,(vh/10))',
  {
      vh: img.select('VH')  // Banda con polarizacion VH
  }).toDouble();
    
  // Resultado
  return vh_log.select([0],['vh'])
  // Seleccionar las propiedades a incluir en la nueva imagen
        .copyProperties(img, ['system: index', 'resolution_meters', 
          'totalSlices', 'productType', 'orbitProperties_pass'
          , 'system:time_start']);
};
    
// Transformación de la banda VV
var s1vv_invlog = function(img){
    
  // Expresión
  var vv_log = img.expression(
    'pow(10,(vv/10))',
  {
      vv: img.select('VV')  // Banda con polarizacion VV
  }).toDouble();
    
  // Resultado
  return vv_log.select([0],['vv'])
  // Seleccionar las propiedades a incluir en la nueva imagen
        .copyProperties(img, ['system: index', 'resolution_meters', 
          'totalSlices', 'productType', 'orbitProperties_pass'
          , 'system:time_start']);
};
  
// Funcion que calcula el cociente VH/VV 
// Las bandas seleccionadas son las creadas con las funciones anteriores
var s1Index = function(img){
    
   // Expresion
   var index = img.expression(
	'vh / vv',
   {
      vh: img.select('vh'),    // Banda con polarizacion VH
      vv: img.select('vv'),    // Banda con polarizacion VV
   }).toDouble();
  
   // Resultado
   return index.select([0], ['vh/vv'])
          .copyProperties(img, ['system: index', 'resolution_meters', 
          'totalSlices', 'productType', 'orbitProperties_pass'
          , 'system:time_start']);
};
  
// Funcion que transforma el resultado del cociente a dB

var s1Trans = function(img){
// Expresion
var index = img.expression(
  '10 * log10(cr)',
  {
      cr: img.select('vh/vv'),    // Banda con el cociente
  }).toDouble();
  
// Resultado
return index.select([0], ['vh/vv_db'])
          .copyProperties(img, ['system: index', 
          'resolution_meters', 'totalSlices', 'productType', 
          'orbitProperties_pass', 'system:time_start']);
};


// --------------- INICIO DE LA FUNCIÓN ----------------- \\
var S1serietemporal = function(coleccion,
                                a_inicio,mes_inicio,dia_inicio,
                                a_final,mes_final,dia_final){

  // Establecer los años de la serie
  var yearrangeStart = a_inicio;
  var yearrangeStop = a_final;
  
  // Generar la coleccion con todas las imágenes
  var col;

  // Funcion para recortar la coleccion
  var clip = function(image){return image.clip(cultivo)};
  // La variable "cultivo" se define más adelante,
  // integra las parcelas cultivadas en la temporada seleccionada
  
// ··········| Cálculo de las series temporales |············\\

  // Loop a través de los años 
  // Crear una serie de cada temporada de cultivo
  for(var loopYear = yearrangeStart; loopYear < yearrangeStop; loopYear +=1){
    // Seleccionar la temporada de producción con la que seleccionar recintos
    var prod_temp = loopYear += 1;
    // Se suma un año, pues los valores de producción utilizados en el 
    // filtro de parcelas están en el año siguiente al comienzo de la temporada
    
    // Se vuelve al año inicial de la temporada
    loopYear -= 1;
    
    // Seleccionar los campos sembrados con cebada en la temporada
    var cultivo = AOI.filter(ee.Filter.gt(prod_temp.toString(), 1));
    
    // Seleccionar las fechas de filtrado (de temporada en temporada)
    var start = ee.Date.fromYMD(loopYear, mes_inicio, dia_inicio);
    var end = ee.Date.fromYMD(loopYear +=1, mes_final, dia_final);

    // Filtrar la colección por fecha y recortar la colección con 
    // los recintos anteriores
    var imgClip = coleccion.filterDate(start, end).map(clip);
    
    // CALCULO DEL COCIENTE VH/VV
    // Aplicar la inversa del logaritmo a cada banda
    var imgClip_vh = imgClip.map(s1vh_invlog);
    var imgClip_vv = imgClip.map(s1vv_invlog);
    
    // Combinar las dos colecciones anteriores en una sola
    var combine = imgClip_vh.combine(imgClip_vv);
    
    // Calcular el cociente VH/VV
    var imgCInv = combine.map(s1Index);
    
    // Transformar a dB
    var imgCIndB = imgCInv.map(s1Trans);
    
    // Volver al año de inicio de la temporada
    loopYear -= 1;
    // Condicional IF:
    // Si el bucle trabaja con el primer año de la serie crea la colección
    if(loopYear == a_inicio){
      col = imgCIndB;
    // El resto de años incluirá las imágenes en la coleccion anterior
    } else if(loopYear > a_inicio && loopYear < a_final){ 
      col = col.merge(imgCIndB);
    }
  } // Fin del bucle FOR

  // Crear el grafico
  var chart = ui.Chart.image.series({
      imageCollection: col,
      region: AOI,
      // Calcular la media de cada imagen
      reducer: ee.Reducer.median(),
      scale: 10
  });
  // Nombrarlo en función de los años de la serie temporal
  var filename = ("des_VHVV_sum_").concat(a_inicio.toString().concat("-")
      .concat(a_final.toString()));
  print(chart, filename);
};

// Aplicar la función anterior a la colección
S1serietemporal(s1Descending,2015,10,1,2020,7,31);