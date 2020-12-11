/*
=== OBTENCION DE MÉTRICAS SOBRE EL ESTADO DE UN CULTIVO ===
                    SENTINEL 2

El código extrae, en una función aplicada a la
colección Sentinel 2, los valores medios de varios índices
de vegetación (IV) en cada una de las imágenes de la colección. 
El resultado es un gráfico que muestra el valor medio de las 
imágenes durante el total de temporadas de cultivo, puede ser 
exportado a CSV. 

Los IV disponibles son:
- NDVI
- NDRE
- IRECI
  
El código recorta las imágenes con los recintos
cultivadas dentro de la temporada de cultivo.

AVISO: Si la serie de años es muy larga es posible que haya 
que ejecutar la función varias veces para mostrar el gráfico
final (ejecutar otra vez siempre que salga el error de tiempo
de procesamiento agotado)

------------------|| INPUTS ||----------------------------
Capa con los recintos cultivados:
  - Nº de columnas equivalente al número de años de la serie 
  (nombre de columna = año)
  - Los valores de las filas se corresponden con la producción
  en cada recinto.

Colección de Sentinel 2. 

Nota: La colección incluida en el este script es la unión de 
dos colecciones. La primera (2015-10-01 hasta 2017-03-31)
esta compuesta por imágenes corregidas a reflectividad BOA
a partir de los productos Sentinel 2 nivel 1C (reflectividad
TOA), mediante el modelo 6s. La segunda es la colección de 
Sentinel 2 nivel 2A (2017-04-01 hasta 2020-07-31). 
El proceso de unión de las dos colecciones es el detallado 
en el apartado "Unión de las colecciones".

------------------|| FUNCIÓN ||---------------------------
S2serietemporal(coleccion, IV, a_inicio,mes_inicio,dia_inicio,
                a_final,mes_final,dia_final);

donde "coleccion": Colección Sentinel 2
"IV": Índice de vegetación a calcular. Elegir uno de los disponibles
"a_inicio": Año de inicio de la primera temporada de
cultivo de la serie
"mes_inicio": Mes de inicio de la temporada
"dia_inicio": Día de inicio de la temporada
"a_final": Año final de la temporada contenida en la 
serie
"mes_final": Mes de finalización de la temporada
"dia_final": Día final de la temporada

Nota: En la zona de estudio la cebada se siembra a comienzos 
de octubre y se termina de cosechar a finales de julio.


---------------|| DESARROLLO DEL CÓDIGO ||------------------
1)  Definir las funciones necesarias en la aplicación de
  la máscara de sombras, modificada de (Schmitt, 2019). 

  - Parámetros generales:
    - Umbrales de las máscaras
    - Alturas de las nubes
  
  - Función "dilatedErossion": Extender la máscara de sombras sobre 
    los píxeles adyacentes

  - Función "projectShadows": Crea una nueva banda en la imagen
    con la probabilidad de sombra (0-1) en cada píxel

2) Integrar los inputs de la función:
  - Capa con los recintos
  - Coleccion de Sentinel 2:
    - Es el resultado de la unión entre las dos colecciones,
      integrando en cada imagen la información sobre la
      probabilidad de nubes procedente de la colección de 
      imágenes 'COPERNICUS/S2_CLOUD_PROBABILITY'

3) Aplicar las máscaras de nubes y sombras

4) Crear las funciones para calcular los IV 

5) Función principal: cálculo de series temporales de IV

  - Se definen los años inicio y final sobre los que calcular el IV.
  - Se define una variable que contendrá las colecciones de imágenes con
    las métricas calculadas de todas las temporadas. 
  - Se crea la función que recorta las imágenes sobre los recintos cultivados 
    en la temporada de adquisición de la propia imagen.

  Bucle FOR: filtra la colección de imágenes por cada temporada de cultivo, 
  recortando las imágenes por los cultivos sembrados en la temporada y calculando 
  el valor medio en cada una de las imágenes. Los pasos se detallan a continuación (se
  repiten durante todas las temporadas de la serie):

0) Se define el iterador como el primer año de la serie. A cada iteración 
   completada se sumará un año, y el bucle se interrumpirá al llegar al 
   último año de la serie.
1) Filtrar la capa con los cultivos de la explotación por el año y los
   valores de producción superiores a 0 (son los campos cultivados)
2) Seleccionar las imágenes dentro de la temporada de cultivo
3) Recortarlas por la capa espacial del paso 1)
4) Aplicar la función para calcular el IV
5) Condicional IF: se vuelve al primer año de la temporada de cultivo del
   bucle y
   a) Si el año coincide con el primer año de la serie la colección de
   imágenes se integra a la colección creada en la primera parte de la función
   b) Si es un año distinto, la colección se suma a las colecciones integradas
   con anterioridad utilizando la función merge()

Al finalizar el bucle se muestra en un gráfico la serie temporal con todos los 
valores medios del indicador durante las temporadas de cultivo filtradas.

6) Llamar a la función principal

Referencias:
Schmitt, M., Hughes, L.H., Qiu, C., Zhu, X.X., 2019. Aggregating 
cloud-free Sentinel-2 images with Google Earth Engine, in: 
ISPRS Annals of Photogrammetry, Remote Sensing and Spatial 
Information Sciences. Presented at the ISPRS ICWG II/III
PIA19+MRSS19 - Photogrammetric Image Analysis & Munich Remote Sensing 
Symposium: Joint ISPRS conference (Volume IV-2/W7) - Copernicus GmbH,
pp. 145–152. https://doi.org/10.5194/isprs-annals-IV-2-W7-145-2019


Gao, B., 1996. NDWI – a normalized difference water index for
remote sensing of vegetation liquid water from space. Remote
Sensing of Environment,58(3), 257 – 266.

Hall, D.K., Riggs, G. A., 2011. Normalized-Difference Snow
Index (NDSI). Springer Netherlands, Dordrecht, 779–780.
---------------------------------------------------------------
*/

// ---------------- MÁSCARA DE NUBES Y SOMBRAS-------------------- \\
// Parámetros generales
var max_probN = 65; // umbral a partir del cuál un píxel es considerado nube
var max_probS = 0.02; // umbral a partir del cuál se considera sombra (tantos por uno)
var ndviThresh = -0.1;  // Umbral de sombras en función del NDVI
var irSumThresh = 0.3;  // Umbral de sombras en función de las bandas infrarrojas

// Parámetros de extensión de la máscara sobre píxeles adyacentes
var erodePixels = 1.5;
var dilationPixels = 3;

// Alturas medias de las nubes
// Utilizadas en la proyección de sus sombras
var cloudHeights = ee.List.sequence(200,10000,250);


function dilatedErossion(score) {
// Se aplica a la capa con la probabilidad de sombras
  score = score
          .reproject('EPSG:4326', null, 20)
          .focal_min({radius: erodePixels, 
            kernelType: 'circle', iterations:3})
          .focal_max({radius: dilationPixels, 
            kernelType: 'circle', iterations:3})
          .reproject('EPSG:4326', null, 20);
                
  return(score);
}
  
// ··········| Cálculo de la probabilidad de sombra |············\\
// ··········|          dentro de un píxel          |············\\

// Función que proyecta las sombras de las nubes
function projectShadows(image){
  // Ángulos de iluminación
  var meanAzimuth = image.get('MEAN_SOLAR_AZIMUTH_ANGLE');
  var meanZenith = image.get('MEAN_SOLAR_ZENITH_ANGLE');
  
  // Banda con la probabilidad de nubes, procedente de la colección
  // "COPERNICUS/S2_CLOUD_PROBABILITY", utilizada para localizar
  // los píxeles con nubes
  var cloudMask = ee.Image(image.get('cloud_mask'))
                  .select('probability').gt(max_probN);
  
    
  // Localizar píxeles oscuros (región infrarroja)
  var darkPixelsImg = image.select(['B8','B11','B12'])
                        .reduce(ee.Reducer.sum());

  // Localizar píxeles oscuros asociados con agua
  var ndvi = image.normalizedDifference(['B8','B4']);
  var waterMask = ndvi.lt(ndviThresh);
  
  // Obtener píxeles de sombra
  var darkPixels = darkPixelsImg.lt(irSumThresh);
    
  // Crear la máscara de sombras excluyendo el agua
  var darkPixelMask = darkPixels.and(waterMask.not());
  darkPixelMask = darkPixelMask.and(cloudMask.not());
     
  // Localizar las sombras de las nubes, basado en la 
  // geometría de iluminación (convertida en radianes)
  var azR = ee.Number(meanAzimuth).add(180).multiply(Math.PI)
            .divide(180.0);
  var zenR = ee.Number(meanZenith).multiply(Math.PI)
            .divide(180.0);
      
  // Localizar las sombras de las nubes
  var shadows = cloudHeights.map(function(cloudHeight){
      cloudHeight = ee.Number(cloudHeight);
     
        var shadowCastedDistance = zenR.tan()
               .multiply(cloudHeight);  // Clasificación según altura
        var x = azR.sin().multiply(shadowCastedDistance)
               .multiply(-1); // distancia de sombras, coord. X
        var y = azR.cos().multiply(shadowCastedDistance)
               .multiply(-1); // distancia de sombras, coord. Y

        // Aplica una máscara sobre los píxeles de sombra
        // desplazando los valores de la máscara de nubes sobre sus sombras
        return ee.Image(image.get('cloud_mask'))
                 .select('probability').displace(ee.Image.constant(x)
        .addBands(ee.Image.constant(y)));
  });
    
  // Agrupar las máscaras de sombras creadas en la función anterior
  // dentro de una colección
  var shadowMasks = ee.ImageCollection.fromImages(shadows);
  var shadowMask = shadowMasks.mean();
    
  // Crear la máscara de sombras
  shadowMask = dilatedErossion(shadowMask.multiply(darkPixelMask));
   
  // Crear la nueba banda con los píxeles sobre sombras
  var shadowScore = shadowMask.reduceNeighborhood(
      {
          reducer: ee.Reducer.max(),
          kernel: ee.Kernel.square(1)
      });
    
  image = image.addBands(shadowScore.rename(['shadowScore']));
    
  return image;
} 

// ---------------------- INPUTS ---------------------- \\
// ···········| Recintos de la explotación |·············\\
var AOI = ee.FeatureCollection(
  'users/iranzocristian/explotacion_blcht_buffer20m');

// ·············|       Colección       |·················\\

// UNIR COLECCIÓN DE GEE Y COLECCIÓN CORREGIDA CON 6S

// Cargar la colección Sentinel 2 (nivel 2A) de GEE
var sentinel2 = ee.ImageCollection("COPERNICUS/S2_SR")
        // Filtrar por area de estudio
        .filterBounds(AOI)

        // Filtrar por fecha (comienzo al acabar la serie de 
        // imágenes corregidas con 6S - S2_6s)
        .filterDate('2017-03-21', '2020-07-31')

        // Primer filtro de nubes
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40))

        // Convertir valores a reflectividad, conservando
        // las propiedades de la imagen original
        .map(function(img){
          var properties = img.propertyNames();
          return(img.divide(10000)
                  .copyProperties(img, properties)) })

        // Ordenar según la fecha de adquisición
        .sort('system:time_start');

// AÑADIR LAS IMÁGENES CON EL PORCENTAJE DE NUBES
// Cargar colección con el porcentaje de nubes
var inicio = ee.Date('2015-10-01');
var final = ee.Date('2020-07-31');

var S2Clouds = ee.ImageCollection(
               'COPERNICUS/S2_CLOUD_PROBABILITY')
               .filterBounds(AOI).filterDate(inicio,final);

// Unir la imagen con el porcentaje de nubes a cada imagen de la
// colección S2 nivel 2A
var sentinel2n = ee.Join.saveFirst('cloud_mask').apply({
    primary: sentinel2,
    secondary: S2Clouds,
    condition:
        ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
});       

// Cargar la colección Sentinel-2 corregida 
// Ubicación de las imágenes
var asset = 'users/iranzocristian/s2_6s/';

// Matriz con el número de imágenes
var image = ee.data.listAssets(asset);

// Seleccionar las imágenes del elemento anterior
var imageList = image['assets'];

// Obtener el número de imágenes
var num_img = imageList.length;
// Restarle uno, el iterador comenzará en el 0
num_img = num_img -= 1;

// Crear una nueva lista con el id de las imágenes
// 1. Definir una lista en blanco
var imageIDs = [];
// 2. Crear iterador
var i;

// 3. Bucle for: integrar en la lista vacía cada uno
// de los id's de las imágenes contenidas en la carpeta,
// variable imageList

for(i=0;i <= num_img; i++){
    // guardar el id
    var imageid = imageList[i].id;
    // incluirlo en la lista
    imageIDs.push(ee.Image(imageid.toString()));
}

// Crear la colección de imágenes a partir de la lista anterior
var S2_6s = ee.ImageCollection.fromImages(imageIDs);

// AÑADIR EL PORCENTAJE DE NUBES 
// Unir la imagen con el porcentaje de nubes a cada imagen de la
// coleccion
var S2_6sn = ee.Join.saveFirst('cloud_mask').apply({
  primary: S2_6s,
  secondary: S2Clouds,
  condition:
      ee.Filter.equals({leftField: 'fileID', rightField: 'system:index'})
});

// UNIR LAS DOS COLECCIONES
var S2_SR = ee.ImageCollection(S2_6sn.merge(sentinel2n));

// ·············| APLICAR LA MÁSCARA DE NUBES |·················\\
// Función que crea la máscara de nubes, en función del umbral
// establecido previamente {max_probN}
var aplicar_maskN = function(img){
  // Crea la máscara con la banda del porcentaje de nubes dentro de la
  // imagen unida a ambas colecciones en los pasos previos
  var mask = ee.Image(img.get('cloud_mask'))
             .select('probability').gt(max_probN).eq(0);
  // Devuelve la imagen con la máscara aplicada
  return img.updateMask(mask);
};

// Aplicar máscara de nubes a la colección
var S2_filtroNubes = S2_SR.map(aplicar_maskN);

// ·············| APLICAR LA MÁSCARA DE SOMBRAS |·················\\
// Obtener la probabilidad de sombra en cada píxel
var probSombras = S2_filtroNubes.map(projectShadows);

// Obtener máscara de sombras, umbral definido previamente {max_probS}
var aplicar_masks = function(img){
var mask = img.select('shadowScore').gt(max_probS).eq(0);
return img.updateMask(mask);
};

// Aplicar máscara de sombras
var S2_maskNS = probSombras.map(aplicar_masks);

// --------------- ÍNDICES DE VEGETACIÓN ----------------- \\

// NDVI
var NDVI = function(image){
  var ndvi = image.normalizedDifference(['B8', 'B4'])
  // Ajustar los valores evitando píxeles saturados y con suelo
  .clamp(0.1, 0.8);
  return ndvi.select([0], ['NDVI'])
         .copyProperties(image, ['PRODUCT_ID',
         'system:index', 'system:time_start']);
};

// NDRE
var NDRE = function(image){
  var ndre = image.normalizedDifference(['B7', 'B4'])
  // Ajustar los valores evitando píxeles saturados y con suelo
  //.clamp(0.1, 0.8);
  return ndre.select([0], ['NDRE'])
         .copyProperties(image, ['PRODUCT_ID',
         'system:index', 'system:time_start']);
};

// IRECI
var IRECI = function(image){
  var ireci = image.expression(
    '(re3-r)/(re1/re2)',
    {
      re3: image.select('B7'),
      r: image.select('B4'),
      re1: image.select('B5'),
      re2: image.select('B6')
    });
  return ireci.select([0], ['IRECI'])
         .copyProperties(image, ['PRODUCT_ID',
         'system:index', 'system:time_start']);
};

// --------------- INICIO DE LA FUNCIÓN ----------------- \\

var S2serietemporal = function(coleccion, IV,
                              a_inicio,mes_inicio,dia_inicio,
                              a_final,mes_final,dia_final){
  
  // Establecer el rango de años de la serie
  var yearrangeStart = a_inicio;
  var yearrangeStop = a_final;
  
  // Generar la coleccion con todas las imágenes
  var col;
     
  // Funcion de recorte, en función de los recintos cultivados
  var clip = function(image){return image.clip(cultivo)};
  // La variable "cultivo" se define más adelante,
  // integra las parcelas cultivadas en el ciclo seleccionado
  
// ··········| Cálculo de las series temporales |············\\
  
  // Loop a través de los años (el código siguiente se ejecuta para cada año)
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
    // La temporada de la cebada comienza en octubre y finaliza en julio
    var start = ee.Date.fromYMD(loopYear, mes_inicio, dia_inicio);
    var end = ee.Date.fromYMD(loopYear +=1, mes_final, dia_final);

    // Filtrar la colección por fecha y recortar la colección con 
    // los recintos anteriores
    var imgClip = coleccion.filterDate(start, end).map(clip);
    
    // Aplicar la función del IV llamado en la función
    var imgIV = imgClip.map(IV);

    // Volver al año de inicio
    loopYear -= 1;

    // Condicional IF:
    // Si el bucle trabaja con el primer año de la serie crea la colección
    // El primer año crea la colección
    if(loopYear == a_inicio){
      col = imgIV;
    // El resto de años incluirá las imágenes en la coleccion anterior
    } else if(loopYear > a_inicio && loopYear < a_final){ 
      col = col.merge(imgIV);
    }
  } // Fin del bucle FOR

// Crear el gráfico
  var chart = ui.Chart.image.series({
    imageCollection: col,
      region: AOI,
      // Calcular la media de cada imagen
      reducer: ee.Reducer.mean(), //opciones= mean, median, stdDev, sum
      scale: 20
    });

// Definir el nombre del gráfico en función de los años de la serie temporal
  var filename = "".concat(a_inicio.toString())
      .concat('-').concat(a_final.toString());
  print(chart);
};// Fin de la función

// --------     LLAMAR A LA FUNCIÓN PRINCIPAL     -------------- \\
S2serietemporal(S2_maskNS, IRECI, 2015,10,1,2020,7,31);
