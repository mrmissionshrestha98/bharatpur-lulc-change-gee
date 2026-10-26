// ============================================================
// MULTI-YEAR LAND USE LAND COVER CLASSIFICATION
// Study Area : Bharatpur Metropolitan City, Nepal
// Years : 2005 (Landsat 5), 2015 (Landsat 8), 2025 (Landsat 9)
// Classifier : Random Forest Supervised Classification
// Classes : Forest, Water Body, Built-up, Agriculture, Other Land
//
// Outputs:
// - LULC maps
// - Accuracy assessment
// - Area by class
// - Transition / change tables
// - Forest gain, loss, persistence
// - Rate of change
// - CSV exports
// - Change maps with proper color coding and legend
// ============================================================
// -----------------------------------------------------------
// 1. STUDY AREA
// -----------------------------------------------------------
69
print(table.first());
var bharatpur = table.filter(
 ee.Filter.eq('GaPa_NaPa', 'Bharatpur')
);
var studyArea = bharatpur;
Map.centerObject(bharatpur, 11);
Map.addLayer(bharatpur, {color: 'red'}, 'Bharatpur Boundary');
var area = bharatpur.geometry().area().divide(1e6);
print('Area of Bharatpur (km²):', area);
// -----------------------------------------------------------
// 2. YEAR-SPECIFIC TRAINING SAMPLES
// -----------------------------------------------------------
// Class codes
// 0 = Forest
// 1 = Water Body
// 2 = Built-up
// 3 = Agriculture
// 4 = Other Land
70
// -- 2005 --
var landcover_2005 = ee.FeatureCollection(forest_2005).map(function(f) {
 return f.set('landcover', 0);
 })
 .merge(ee.FeatureCollection(water_body_2005).map(function(f) {
 return f.set('landcover', 1);
 }))
 .merge(ee.FeatureCollection(urban_2005).map(function(f) {
 return f.set('landcover', 2);
 }))
 .merge(ee.FeatureCollection(agriculture_land_2005).map(function(f) {
 return f.set('landcover', 3);
 }))
 .merge(ee.FeatureCollection(otherLand_2005).map(function(f) {
 return f.set('landcover', 4);
 }));
// -- 2015 --
var landcover_2015 = ee.FeatureCollection(forest_2015).map(function(f) {
 return f.set('landcover', 0);
 })
 .merge(ee.FeatureCollection(water_body_2015).map(function(f) {
71
 return f.set('landcover', 1);
 }))
 .merge(ee.FeatureCollection(urban_2015).map(function(f) {
 return f.set('landcover', 2);
 }))
 .merge(ee.FeatureCollection(agriculture_land_2015).map(function(f) {
 return f.set('landcover', 3);
 }))
 .merge(ee.FeatureCollection(otherLand_2015).map(function(f) {
 return f.set('landcover', 4);
 }));
// -- 2025 --
var landcover_2025 = ee.FeatureCollection(forest_2025).map(function(f) {
 return f.set('landcover', 0);
 })
 .merge(ee.FeatureCollection(water_body_2025).map(function(f) {
 return f.set('landcover', 1);
 }))
 .merge(ee.FeatureCollection(urban_2025).map(function(f) {
 return f.set('landcover', 2);
 }))
 .merge(ee.FeatureCollection(agriculture_land_2025).map(function(f) {
72
 return f.set('landcover', 3);
 }))
 .merge(ee.FeatureCollection(otherLand_2025).map(function(f) {
 return f.set('landcover', 4);
 }));
// -----------------------------------------------------------
// 3. CONFIGURATION
// -----------------------------------------------------------
var classProperty = 'landcover';
var commonBands = [
 'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2',
 'NDVI', 'NDBI', 'MNDWI', 'NDMI', 'BSI', 'IBI', 'SAVI'
];
var palette = [
 '008080', // Forest
 '7B68EE', // Water Body
 'FF0000', // Built-up
 'FFFF00', // Agriculture
 '808080' // Other Land
73
];
var classNames = ['Forest', 'Water Body', 'Built-up', 'Agriculture', 'Other Land'];
var classNamesEE = ee.List(classNames);
var classOrder = ee.List([0, 1, 2, 3, 4]);
var vizParams = {min: 0, max: 4, palette: palette};
var rgbViz = {bands: ['Red', 'Green', 'Blue'], min: 0, max: 0.3, gamma: 1.4};
// -----------------------------------------------------------
// 4. HELPER FUNCTIONS
// -----------------------------------------------------------
// ---------- Cloud mask for Landsat Collection 2 Level-2 ----------
function maskLandsatC2L2(img) {
 var qa = img.select('QA_PIXEL');
 var cloud = qa.bitwiseAnd(1 << 3).neq(0);
 var shadow = qa.bitwiseAnd(1 << 4).neq(0);
 var snow = qa.bitwiseAnd(1 << 5).neq(0);
 var cirrus = qa.bitwiseAnd(1 << 7).neq(0);
74
 var mask = cloud.or(shadow).or(snow).or(cirrus).not();
 var sat = img.select('QA_RADSAT').eq(0);
 return img.updateMask(mask).updateMask(sat);
}
// ---------- Build SR composite ----------
function getLandsatSRComposite(collectionId, srcBands, startDate, endDate) {
 var col = ee.ImageCollection(collectionId)
 .filterBounds(studyArea)
 .filterDate(startDate, endDate)
 .map(maskLandsatC2L2);
 return col.median()
 .select(srcBands, ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
 .multiply(0.0000275)
 .add(-0.2)
 .clip(studyArea);
}
// ---------- Add spectral indices ----------
function addIndices(img) {
 var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
75
 var ndbi = img.normalizedDifference(['SWIR1', 'NIR']).rename('NDBI');
 var mndwi = img.normalizedDifference(['Green', 'SWIR1']).rename('MNDWI');
 var ndmi = img.normalizedDifference(['NIR', 'SWIR1']).rename('NDMI');
 var bsi = img.expression(
 '((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue))', {
 'SWIR1': img.select('SWIR1'),
 'Red' : img.select('Red'),
 'NIR' : img.select('NIR'),
 'Blue' : img.select('Blue')
 }).rename('BSI');
 var ibi = img.expression(
 '(ndbi - ((ndvi + mndwi) / 2)) / (ndbi + ((ndvi + mndwi) / 2))', {
 'ndbi': ndbi,
 'ndvi': ndvi,
 'mndwi': mndwi
 }).rename('IBI');
 var savi = img.expression(
 '((NIR - Red) / (NIR + Red + 0.5)) * 1.5', {
 'NIR': img.select('NIR'),
 'Red': img.select('Red')
76
 }).rename('SAVI');
 return img.addBands([ndvi, ndbi, mndwi, ndmi, bsi, ibi, savi]);
}
// ---------- Water refinement ----------
function refineWater(classified, image) {
 var waterMask = image.select('MNDWI').gt(0.15)
 .and(image.select('NDVI').lt(0.1));
 return classified.where(waterMask, 1);
}
// ---------- Built-up refinement ----------
function refineBuiltUp(classified, image) {
 var urbanMask = image.select('NDBI').gt(0.2)
 .and(image.select('IBI').gt(0.05))
 .and(image.select('NDVI').lt(0.25))
 .and(image.select('MNDWI').lt(0));
 return classified.where(urbanMask, 2);
}
// ---------- Strong urban refinement ----------
function strongUrbanRefinement(classified, image) {
77
 var urbanMask = image.select('NDBI').gt(0.25)
 .and(image.select('IBI').gt(0.1))
 .and(image.select('NDVI').lt(0.2))
 .and(image.select('BSI').lt(0.35));
 return classified.where(urbanMask, 2);
}
// ---------- Smoothing ----------
function smoothClassified(img) {
 return img.focal_mode({radius: 1, units: 'pixels'});
}
// ---------- Per-class accuracy ----------
function getPerClassAcc(validatedSamples, yearLabel) {
 var features = classOrder.map(function(i) {
 i = ee.Number(i);
 var refClass = validatedSamples.filter(ee.Filter.eq(classProperty, i));
 var predClass = validatedSamples.filter(ee.Filter.eq('class', i));
 var correctClass = validatedSamples
 .filter(ee.Filter.eq(classProperty, i))
 .filter(ee.Filter.eq('class', i));
78
 var reference_total = refClass.size();
 var predicted_total = predClass.size();
 var correct = correctClass.size();
 var producers_accuracy = ee.Algorithms.If(
 ee.Number(reference_total).gt(0),
 ee.Number(correct).divide(reference_total),
 null
 );
 var users_accuracy = ee.Algorithms.If(
 ee.Number(predicted_total).gt(0),
 ee.Number(correct).divide(predicted_total),
 null
 );
 return ee.Feature(null, {
 year : yearLabel,
 class_id : i,
 class_name : classNamesEE.get(i),
 correct : correct,
 reference_total : reference_total,
 predicted_total : predicted_total,
79
 producers_accuracy: producers_accuracy,
 users_accuracy : users_accuracy
 });
 });
 var fc = ee.FeatureCollection(features);
 print('Per-class accuracy table — ' + yearLabel, fc);
 return fc;
}
// ---------- Area by class ----------
function getAreaByClass(classifiedImage, yearLabel) {
 var areaImage = ee.Image.pixelArea().divide(1e6).rename('area')
 .addBands(classifiedImage.rename('class'));
 var areaByClass = areaImage.reduceRegion({
 reducer : ee.Reducer.sum().group({ groupField: 1, groupName: 'class' }),
 geometry : studyArea.geometry(),
 scale : 30,
 maxPixels: 1e13
 });
 var groups = ee.List(areaByClass.get('groups'));
80
 var rows = groups.map(function(g) {
 g = ee.Dictionary(g);
 var c = ee.Number(g.get('class')).toInt();
 return ee.Feature(null, {
 year : yearLabel,
 class_id : c,
 class_name: classNamesEE.get(c),
 area_km2 : ee.Number(g.get('sum'))
 });
 });
 var fc = ee.FeatureCollection(rows).sort('class_id');
 print('Area by class (km²) — ' + yearLabel, fc);
 return fc;
}
// -----------------------------------------------------------
// 5. TRAIN, CLASSIFY & VALIDATE
// -----------------------------------------------------------
function trainClassifyAndValidate(image, trainingPoints, yearLabel) {
 print('==============================');
81
 print('Classification Year: ' + yearLabel);
 print('==============================');
 var samples = image.sampleRegions({
 collection: trainingPoints,
 properties: [classProperty],
 scale : 30,
 geometries: true,
 tileScale : 8
 });
 var needed = commonBands.concat([classProperty]);
 var samplesClean = samples.filter(ee.Filter.notNull(needed));
 print('Sampled class histogram (' + yearLabel + '):',
samplesClean.aggregate_histogram(classProperty));
 var withRandom = samplesClean.randomColumn('random', 42);
 var trainSet = withRandom.filter(ee.Filter.lt('random', 0.7));
 var testSet = withRandom.filter(ee.Filter.gte('random', 0.7));
 print('Total samples :', samplesClean.size());
 print('Train samples :', trainSet.size());
 print('Test samples :', testSet.size());
82
 print('Train class hist:', trainSet.aggregate_histogram(classProperty));
 print('Test class hist :', testSet.aggregate_histogram(classProperty));
 var classifier = ee.Classifier.smileRandomForest({
 numberOfTrees : 200,
 bagFraction : 0.6,
 minLeafPopulation: 2,
 seed : 42
 }).train({
 features : trainSet,
 classProperty : classProperty,
 inputProperties: commonBands
 });
 var trainCM = classifier.confusionMatrix();
 print('-- Training Accuracy --');
 print('Confusion Matrix :', trainCM);
 print('Overall Accuracy :', trainCM.accuracy());
 print('Kappa :', trainCM.kappa());
 var classified = image.classify(classifier).rename('class');
 classified = refineWater(classified, image);
83
 classified = refineBuiltUp(classified, image);
 classified = strongUrbanRefinement(classified, image);
 classified = smoothClassified(classified).rename('class').toByte();
 var validatedFinal = classified.sampleRegions({
 collection: testSet,
 properties: [classProperty],
 scale : 30,
 geometries: false,
 tileScale : 8
 });
 var testCM = validatedFinal.errorMatrix(classProperty, 'class', classOrder);
 print('-- Validation Accuracy (final refined map) --');
 print('Confusion Matrix :', testCM);
 print('Overall Accuracy :', testCM.accuracy());
 print('Kappa :', testCM.kappa());
 print('Producers Acc :', testCM.producersAccuracy());
 print('Users Acc :', testCM.consumersAccuracy());
 var perClassAcc = getPerClassAcc(validatedFinal, yearLabel);
 print('Order of classes :', classNames);
84
 return {
 classified : classified,
 accuracyTable: perClassAcc
 };
}
// -----------------------------------------------------------
// 6. LOAD IMAGERY
// -----------------------------------------------------------
var image2005 = getLandsatSRComposite(
 'LANDSAT/LT05/C02/T1_L2',
 ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
 '2005-09-01', '2005-12-31'
);
var image2015 = getLandsatSRComposite(
 'LANDSAT/LC08/C02/T1_L2',
 ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
 '2015-09-01', '2015-12-31'
);
85
var image2025 = getLandsatSRComposite(
 'LANDSAT/LC09/C02/T1_L2',
 ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
 '2025-09-01', '2025-12-31'
);
image2005 = addIndices(image2005);
image2015 = addIndices(image2015);
image2025 = addIndices(image2025);
Map.addLayer(image2005, rgbViz, 'RGB 2005');
Map.addLayer(image2015, rgbViz, 'RGB 2015');
Map.addLayer(image2025, rgbViz, 'RGB 2025');
print('2005 image bands:', image2005.bandNames());
print('2015 image bands:', image2015.bandNames());
print('2025 image bands:', image2025.bandNames());
// -----------------------------------------------------------
// 7. CLASSIFY EACH YEAR
// -----------------------------------------------------------
var result2005 = trainClassifyAndValidate(image2005, landcover_2005, '2005');
86
var classified2005 = result2005.classified.toByte();
var accuracy2005 = result2005.accuracyTable;
Map.addLayer(classified2005, vizParams, 'LULC 2005');
var area2005 = getAreaByClass(classified2005, '2005');
var result2015 = trainClassifyAndValidate(image2015, landcover_2015, '2015');
var classified2015 = result2015.classified.toByte();
var accuracy2015 = result2015.accuracyTable;
Map.addLayer(classified2015, vizParams, 'LULC 2015');
var area2015 = getAreaByClass(classified2015, '2015');
var result2025 = trainClassifyAndValidate(image2025, landcover_2025, '2025');
var classified2025 = result2025.classified.toByte();
var accuracy2025 = result2025.accuracyTable;
Map.addLayer(classified2025, vizParams, 'LULC 2025');
var area2025 = getAreaByClass(classified2025, '2025');
// -----------------------------------------------------------
// 8. TRANSITION MATRIX / CHANGE TABLE
// -----------------------------------------------------------
function transitionMatrixArea(fromImage, toImage, yearLabel) {
 var transition = fromImage.toInt16().multiply(10)
87
 .add(toImage.toInt16())
 .toInt16()
 .rename('transition');
 var areaImage = ee.Image.pixelArea().divide(1e6).rename('area')
 .addBands(transition);
 var transitionStats = areaImage.reduceRegion({
 reducer : ee.Reducer.sum().group({ groupField: 1, groupName: 'transition' }),
 geometry : studyArea.geometry(),
 scale : 30,
 maxPixels: 1e13,
 tileScale: 16
 });
 var groups = ee.List(transitionStats.get('groups'));
 var rows = groups.map(function(g) {
 g = ee.Dictionary(g);
 var code = ee.Number(g.get('transition')).toInt();
 var fromClass = code.divide(10).floor();
 var toClass = code.mod(10);
88
 return ee.Feature(null, {
 period : yearLabel,
 transition_code : code,
 from_class_id : fromClass,
 from_class_name : classNamesEE.get(fromClass),
 to_class_id : toClass,
 to_class_name : classNamesEE.get(toClass),
 area_km2 : ee.Number(g.get('sum'))
 });
 });
 var fc = ee.FeatureCollection(rows).sort('transition_code');
 print('Transition Matrix / Change Table (km²) — ' + yearLabel, fc);
 return fc;
}
var transition2005_2015 = transitionMatrixArea(classified2005, classified2015, '2005-
2015');
var transition2015_2025 = transitionMatrixArea(classified2015, classified2025, '2015-
2025');
var transition2005_2025 = transitionMatrixArea(classified2005, classified2025, '2005-
2025');
89
// -----------------------------------------------------------
// 9. FOREST GAIN / LOSS / PERSISTENCE
// -----------------------------------------------------------
function forestChangeAnalysis(img1, img2, yearLabel) {
 var forest1 = img1.eq(0);
 var forest2 = img2.eq(0);
 var forestLoss = forest1.and(forest2.not()).rename('forest_loss').toByte();
 var forestGain = forest1.not().and(forest2).rename('forest_gain').toByte();
 var forestStable = forest1.and(forest2).rename('forest_persistence').toByte();
 var areaImage = ee.Image.pixelArea().divide(1e6);
 var lossArea = forestLoss.multiply(areaImage).reduceRegion({
 reducer: ee.Reducer.sum(), geometry: studyArea.geometry(), scale: 30, maxPixels: 1e13
 });
 var gainArea = forestGain.multiply(areaImage).reduceRegion({
 reducer: ee.Reducer.sum(), geometry: studyArea.geometry(), scale: 30, maxPixels: 1e13
 });
 var stableArea = forestStable.multiply(areaImage).reduceRegion({
 reducer: ee.Reducer.sum(), geometry: studyArea.geometry(), scale: 30, maxPixels: 1e13
 });
90
 var forestSummary = ee.FeatureCollection([
 ee.Feature(null, { period: yearLabel, metric: 'Forest Loss', area_km2:
ee.Number(lossArea.get('forest_loss')) }),
 ee.Feature(null, { period: yearLabel, metric: 'Forest Gain', area_km2:
ee.Number(gainArea.get('forest_gain')) }),
 ee.Feature(null, { period: yearLabel, metric: 'Forest Persistence', area_km2:
ee.Number(stableArea.get('forest_persistence')) })
 ]);
 print('Forest Change Summary (km²) — ' + yearLabel, forestSummary);
 Map.addLayer(forestLoss.selfMask(), {palette: ['FF0000']}, 'Forest Loss ' + yearLabel);
 Map.addLayer(forestGain.selfMask(), {palette: ['00FF00']}, 'Forest Gain ' + yearLabel);
 Map.addLayer(forestStable.selfMask(), {palette: ['006400']}, 'Forest Persistence ' +
yearLabel);
 return { loss: forestLoss, gain: forestGain, stable: forestStable, summary: forestSummary
};
}
var forest2005_2015 = forestChangeAnalysis(classified2005, classified2015, '2005-2015');
var forest2015_2025 = forestChangeAnalysis(classified2015, classified2025, '2015-2025');
var forest2005_2025 = forestChangeAnalysis(classified2005, classified2025, '2005-2025');
91
// -----------------------------------------------------------
// 10. RATE OF CHANGE BY CLASS
// -----------------------------------------------------------
function compareAreaChange(areaFC1, areaFC2, year1, year2, intervalYears) {
 var joined = classOrder.map(function(i) {
 i = ee.Number(i);
 var f1 = areaFC1.filter(ee.Filter.eq('class_id', i)).first();
 var f2 = areaFC2.filter(ee.Filter.eq('class_id', i)).first();
 var area1 = ee.Number(ee.Algorithms.If(f1, ee.Feature(f1).get('area_km2'), 0));
 var area2 = ee.Number(ee.Algorithms.If(f2, ee.Feature(f2).get('area_km2'), 0));
 var netChange = area2.subtract(area1);
 var annualRatePercent = ee.Algorithms.If(
 area1.neq(0),
 netChange.divide(area1).divide(intervalYears).multiply(100),
 null
 );
 return ee.Feature(null, {
 from_year : year1,
 to_year : year2,
92
 class_id : i,
 class_name : classNamesEE.get(i),
 area_start_km2 : area1,
 area_end_km2 : area2,
 net_change_km2 : netChange,
 annual_rate_percent: annualRatePercent
 });
 });
 var fc = ee.FeatureCollection(joined);
 print('Rate of Change by Class — ' + year1 + '-' + year2, fc);
 return fc;
}
var rate2005_2015 = compareAreaChange(area2005, area2015, '2005', '2015', 10);
var rate2015_2025 = compareAreaChange(area2015, area2025, '2015', '2025', 10);
var rate2005_2025 = compareAreaChange(area2005, area2025, '2005', '2025', 20);
// -----------------------------------------------------------
// 11. CHANGE DETECTION MAPS (proper color coding)
// -----------------------------------------------------------
93
// Transition code = fromClass * 10 + toClass
// All 25 possible transitions across 5 classes
var change_2005_2015 = classified2005.toInt16()
 .multiply(10).add(classified2015.toInt16()).toInt16().rename('change');
var change_2015_2025 = classified2015.toInt16()
 .multiply(10).add(classified2025.toInt16()).toInt16().rename('change');
var change_2005_2025 = classified2005.toInt16()
 .multiply(10).add(classified2025.toInt16()).toInt16().rename('change');
// Full list of all 25 transition codes (fromClass * 10 + toClass)
var transitionCodes = [
 0, 1, 2, 3, 4, // Forest → [Forest, Water, Built-up, Agriculture, Other]
 10, 11, 12, 13, 14, // Water → [Forest, Water, Built-up, Agriculture, Other]
 20, 21, 22, 23, 24, // Built-up → [Forest, Water, Built-up, Agriculture, Other]
 30, 31, 32, 33, 34, // Agriculture → [Forest, Water, Built-up, Agriculture, Other]
 40, 41, 42, 43, 44 // Other Land → [Forest, Water, Built-up, Agriculture, Other]
];
// Remap to sequential 0-24 for palette indexing
var remappedValues = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
94
// Color scheme logic:
// - Diagonal (stable classes): deep/saturated version of original class color
// - Forest loss (row 0, off-diagonal): warm reds/oranges — ecologically alarming
// - Forest gain (column 0, off-diagonal): cool greens — positive change
// - Urban expansion (any → Built-up, off-diagonal): red tones
// - Other transitions: muted/neutral tones
var transitionPalette = [
 // Row 0: Forest as FROM
 '006400', // 00: Forest → Forest (stable dark green)
 '4169E1', // 01: Forest → Water (royal blue)
 'FF0000', // 02: Forest → Built-up (red — deforestation/urbanization)
 'FFA500', // 03: Forest → Agriculture (orange — deforestation/agriculture)
 'A9A9A9', // 04: Forest → Other (grey — degradation)
 // Row 1: Water as FROM
 '90EE90', // 10: Water → Forest (light green)
 '7B68EE', // 11: Water → Water (stable medium purple)
 'FF69B4', // 12: Water → Built-up (hot pink — water loss to urban)
 'FFD700', // 13: Water → Agriculture (gold)
 'C0C0C0', // 14: Water → Other (silver)
 // Row 2: Built-up as FROM
 '00C800', // 20: Built-up → Forest (bright green — urban greening)
95
 '00BFFF', // 21: Built-up → Water (deep sky blue)
 '8B0000', // 22: Built-up → Built-up (stable dark red)
 'DAA520', // 23: Built-up → Agriculture (goldenrod)
 'B8B8B8', // 24: Built-up → Other (light grey)
 // Row 3: Agriculture as FROM
 '32CD32', // 30: Agriculture → Forest (lime green — reforestation)
 '1E90FF', // 31: Agriculture → Water (dodger blue)
 'DC143C', // 32: Agriculture → Built-up (crimson — agri loss to urban)
 'FFFF00', // 33: Agriculture → Agriculture (stable yellow)
 'D3D3D3', // 34: Agriculture → Other (light grey)
 // Row 4: Other Land as FROM
 '228B22', // 40: Other → Forest (forest green)
 '6495ED', // 41: Other → Water (cornflower blue)
 'CD5C5C', // 42: Other → Built-up (indian red)
 'F0E68C', // 43: Other → Agriculture (khaki)
 '808080' // 44: Other → Other (stable grey)
];
var changeViz = {
 min : 0,
 max : 24,
96
 palette: transitionPalette
};
// Remap each change image for display
function colorChangeMap(changeImg) {
 return changeImg.remap(transitionCodes,
remappedValues).rename('change_remapped');
}
var changeRemapped_2005_2015 = colorChangeMap(change_2005_2015);
var changeRemapped_2015_2025 = colorChangeMap(change_2015_2025);
var changeRemapped_2005_2025 = colorChangeMap(change_2005_2025);
Map.addLayer(changeRemapped_2005_2015, changeViz, 'Change Map 2005-2015');
Map.addLayer(changeRemapped_2015_2025, changeViz, 'Change Map 2015-2025');
Map.addLayer(changeRemapped_2005_2025, changeViz, 'Change Map 2005-2025');
// -----------------------------------------------------------
// 12. LEGENDS
// -----------------------------------------------------------
// ---- LULC Legend ----
var lulcLegend = ui.Panel({
97
 style: {position: 'bottom-left', padding: '8px 15px'}
});
lulcLegend.add(ui.Label({
 value: 'Land Cover Classes',
 style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 6px 0'}
}));
for (var i = 0; i < classNames.length; i++) {
 var colorBox = ui.Label({
 style: {
 backgroundColor: '#' + palette[i],
 padding : '8px',
 margin : '0 4px 4px 0'
 }
 });
 var label = ui.Label({
 value: classNames[i],
 style: {margin: '0 0 4px 6px'}
 });
 lulcLegend.add(ui.Panel([colorBox, label], ui.Panel.Layout.Flow('horizontal')));
}
98
Map.add(lulcLegend);
// ---- Change Map Legend (key transitions only) ----
var changeLegend = ui.Panel({
 style: {position: 'bottom-right', padding: '8px 15px'}
});
changeLegend.add(ui.Label({
 value: 'Key LULC Transitions',
 style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 6px 0'}
}));
var keyTransitions = [
 {label: 'Forest → Forest (stable)', color: '006400'},
 {label: 'Forest → Built-up', color: 'FF0000'},
 {label: 'Forest → Agriculture', color: 'FFA500'},
 {label: 'Forest → Other Land', color: 'A9A9A9'},
 {label: 'Water → Water (stable)', color: '7B68EE'},
 {label: 'Built-up → Built-up (stable)', color: '8B0000'},
 {label: 'Built-up → Forest', color: '00C800'},
 {label: 'Agriculture → Built-up', color: 'DC143C'},
 {label: 'Agriculture → Agriculture (stable)', color: 'FFFF00'},
 {label: 'Agriculture → Forest', color: '32CD32'},
99
 {label: 'Other → Other (stable)', color: '808080'},
];
keyTransitions.forEach(function(item) {
 var row = ui.Panel([
 ui.Label({style: {backgroundColor: '#' + item.color, padding: '8px', margin: '0 4px 4px 0'}}),
 ui.Label({value: item.label, style: {margin: '0 0 4px 6px', fontSize: '11px'}})
 ], ui.Panel.Layout.Flow('horizontal'));
 changeLegend.add(row);
});
Map.add(changeLegend);
// -----------------------------------------------------------
// 13. EXPORT LULC MAPS
// -----------------------------------------------------------
Export.image.toDrive({
 image : classified2005.toByte(),
 description : 'LULC_Bharatpur_2005_raw',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'LULC_Bharatpur_2005_raw',
 scale : 30,
100
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : classified2015.toByte(),
 description : 'LULC_Bharatpur_2015_raw',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'LULC_Bharatpur_2015_raw',
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : classified2025.toByte(),
 description : 'LULC_Bharatpur_2025_raw',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'LULC_Bharatpur_2025_raw',
 scale : 30,
 region : studyArea.geometry(),
101
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : classified2005.visualize({min: 0, max: 4, palette: palette}).toByte(),
 description : 'LULC_Bharatpur_2005_color',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'LULC_Bharatpur_2005_color',
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : classified2015.visualize({min: 0, max: 4, palette: palette}).toByte(),
 description : 'LULC_Bharatpur_2015_color',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'LULC_Bharatpur_2015_color',
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
102
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : classified2025.visualize({min: 0, max: 4, palette: palette}).toByte(),
 description : 'LULC_Bharatpur_2025_color',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'LULC_Bharatpur_2025_color',
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
// -----------------------------------------------------------
// 14. EXPORT CHANGE MAPS (color-coded, no randomVisualizer)
// -----------------------------------------------------------
Export.image.toDrive({
 image : change_2005_2015.toInt16(),
 description : 'ChangeMap_2005_2015_raw',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'ChangeMap_2005_2015_raw',
103
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : change_2015_2025.toInt16(),
 description : 'ChangeMap_2015_2025_raw',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'ChangeMap_2015_2025_raw',
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : change_2005_2025.toInt16(),
 description : 'ChangeMap_2005_2025_raw',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'ChangeMap_2005_2025_raw',
 scale : 30,
104
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : changeRemapped_2005_2015.visualize(changeViz).toByte(),
 description : 'ChangeMap_2005_2015_color',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'ChangeMap_2005_2015_color',
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : changeRemapped_2015_2025.visualize(changeViz).toByte(),
 description : 'ChangeMap_2015_2025_color',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'ChangeMap_2015_2025_color',
 scale : 30,
 region : studyArea.geometry(),
105
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
Export.image.toDrive({
 image : changeRemapped_2005_2025.visualize(changeViz).toByte(),
 description : 'ChangeMap_2005_2025_color',
 folder : 'GEE_LULC_Bharatpur',
 fileNamePrefix: 'ChangeMap_2005_2025_color',
 scale : 30,
 region : studyArea.geometry(),
 maxPixels : 1e13,
 crs : 'EPSG:32644'
});
// -----------------------------------------------------------
// 15. EXPORT FOREST CHANGE MAPS
// -----------------------------------------------------------
Export.image.toDrive({
 image: forest2005_2015.loss.toByte(), description: 'ForestLoss_2005_2015',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestLoss_2005_2015',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
106
});
Export.image.toDrive({
 image: forest2005_2015.gain.toByte(), description: 'ForestGain_2005_2015',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestGain_2005_2015',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
Export.image.toDrive({
 image: forest2005_2015.stable.toByte(), description: 'ForestPersistence_2005_2015',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestPersistence_2005_2015',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
Export.image.toDrive({
 image: forest2015_2025.loss.toByte(), description: 'ForestLoss_2015_2025',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestLoss_2015_2025',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
Export.image.toDrive({
 image: forest2015_2025.gain.toByte(), description: 'ForestGain_2015_2025',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestGain_2015_2025',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
Export.image.toDrive({
107
 image: forest2015_2025.stable.toByte(), description: 'ForestPersistence_2015_2025',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestPersistence_2015_2025',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
Export.image.toDrive({
 image: forest2005_2025.loss.toByte(), description: 'ForestLoss_2005_2025',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestLoss_2005_2025',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
Export.image.toDrive({
 image: forest2005_2025.gain.toByte(), description: 'ForestGain_2005_2025',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestGain_2005_2025',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
Export.image.toDrive({
 image: forest2005_2025.stable.toByte(), description: 'ForestPersistence_2005_2025',
 folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'ForestPersistence_2005_2025',
 scale: 30, region: studyArea.geometry(), maxPixels: 1e13, crs: 'EPSG:32644'
});
// -----------------------------------------------------------
108
// 16. EXPORT TABLES AS CSV
// -----------------------------------------------------------
Export.table.toDrive({ collection: area2005, description: 'AreaByClass_2005_csv',
folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'AreaByClass_2005', fileFormat:
'CSV' });
Export.table.toDrive({ collection: area2015, description: 'AreaByClass_2015_csv',
folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'AreaByClass_2015', fileFormat:
'CSV' });
Export.table.toDrive({ collection: area2025, description: 'AreaByClass_2025_csv',
folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'AreaByClass_2025', fileFormat:
'CSV' });
Export.table.toDrive({ collection: accuracy2005, description: 'AccuracyTable_2005_csv',
folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'AccuracyTable_2005', fileFormat:
'CSV' });
Export.table.toDrive({ collection: accuracy2015, description: 'AccuracyTable_2015_csv',
folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'AccuracyTable_2015', fileFormat:
'CSV' });
Export.table.toDrive({ collection: accuracy2025, description: 'AccuracyTable_2025_csv',
folder: 'GEE_LULC_Bharatpur', fileNamePrefix: 'AccuracyTable_2025', fileFormat:
'CSV' });
Export.table.toDrive({ collection: transition2005_2015, description:
'TransitionMatrix_2005_2015_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'TransitionMatrix_2005_2015', fileFormat: 'CSV' });
Export.table.toDrive({ collection: transition2015_2025, description:
'TransitionMatrix_2015_2025_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'TransitionMatrix_2015_2025', fileFormat: 'CSV' });
109
Export.table.toDrive({ collection: transition2005_2025, description:
'TransitionMatrix_2005_2025_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'TransitionMatrix_2005_2025', fileFormat: 'CSV' });
Export.table.toDrive({ collection: forest2005_2015.summary, description:
'ForestChangeSummary_2005_2015_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'ForestChangeSummary_2005_2015', fileFormat: 'CSV' });
Export.table.toDrive({ collection: forest2015_2025.summary, description:
'ForestChangeSummary_2015_2025_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'ForestChangeSummary_2015_2025', fileFormat: 'CSV' });
Export.table.toDrive({ collection: forest2005_2025.summary, description:
'ForestChangeSummary_2005_2025_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'ForestChangeSummary_2005_2025', fileFormat: 'CSV' });
Export.table.toDrive({ collection: rate2005_2015, description:
'RateOfChange_2005_2015_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'RateOfChange_2005_2015', fileFormat: 'CSV' });
Export.table.toDrive({ collection: rate2015_2025, description:
'RateOfChange_2015_2025_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'RateOfChange_2015_2025', fileFormat: 'CSV' });
Export.table.toDrive({ collection: rate2005_2025, description:
'RateOfChange_2005_2025_csv', folder: 'GEE_LULC_Bharatpur', fileNamePrefix:
'RateOfChange_2005_2025', fileFormat: 'CSV' });
