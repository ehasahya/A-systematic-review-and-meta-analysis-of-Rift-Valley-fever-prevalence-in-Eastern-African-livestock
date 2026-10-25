/******************************************************************************
 * Ecological classification of sampling districts 
 * A systematic review and meta-analysis of Rift Valley fever prevalence in 
 * Eastern African livestock: associations with ecological 
 * and epidemiological factors
 * ----------------------------------------------------------------------------
 * WHAT THIS SCRIPT DOES
 *   For each administrative unit (district) in which livestock were sampled in
 *   the included studies, this script derives four ecological variables from
 *   public geospatial datasets in Google Earth Engine and classifies the
 *   district as "yes" or "no" for each. Districts are selected by their name
 *   in the FAO GAUL 2015 boundaries. The output contains only district-level 
 *   ecological information. 
 *
 * VARIABLES AND CLASSIFICATION RULES (Supplementary Table 3)
 *
 *   Flood-prone
 *     Dataset : Global Flood Database, MODIS 250 m flood events 2000-2018
 *               (GLOBAL_FLOOD_DB/MODIS_EVENTS/V1); permanent water masked
 *               with JRC Global Surface Water (occurrence >= 50%).
 *     "yes"   : at least one flood event recorded within the district polygon
 *               (flood_max_events >= 1).
 *
 *   Urban setting
 *     Dataset : GHSL Degree of Urbanisation settlement model, GHS-SMOD,
 *               epoch 2020 (JRC/GHSL/P2023A/GHS_SMOD).
 *     "urban" : SMOD classes 23 (dense urban cluster) and 30 (urban centre)
 *               together cover more than 10% of the district's land area
 *               (urban_23_30 > 10); otherwise it is classifed as "rural".
 *
 *   Near national park (wildlife protected area)
 *     Dataset : World Database on Protected Areas (WCMC/WDPA/current/polygons).
 *     "yes"   : the district polygon intersects, or lies within 6 km of, a
 *               wildlife protected area (wild_dist_km <= 6).
 *
 *   Near international border
 *     Dataset : US Department of State Large-Scale International Boundaries,
 *               LSIB 2017 (USDOS/LSIB_SIMPLE/2017); land boundaries only,
 *               coastlines excluded.
 *     "yes"   : the district polygon touches an international land boundary
 *               (border_dist_km = 0).
 *
 *   The thresholds are set in the "classification thresholds" block below and
 *   can be changed to test other cut-offs. Prior RVF outbreaks (FAO EMPRES-i)
 *   were classified separately and are not part of this script.
 *
 * HOW TO RUN
 *   1. Sign in to the Earth Engine Code Editor (code.earthengine.google.com);
 *      a free Earth Engine account is required.
 *   2. Paste this script into a new script and click Run.
 *   3. Check the Console: "Units NOT matched to exactly one polygon" must be
 *      empty. 
 *   4. Open the Tasks tab and click Run next to 'RVF_eco_v11_units' to export
 *      the results as a CSV file to your Google Drive.
 *
 * INPUT (UNITS list below)
 *   One row per district:
 *     unit_id      : district identifier
 *     country      : country name as used in the data workbook
 *     gaul_country : country name as spelt in GAUL (ADM0_NAME)
 *     level        : 1 = GAUL level 1 (Uganda districts); 2 = GAUL level 2
 *                    (districts in all other countries)
 *     name         : district name exactly as spelt in GAUL
 *
 * OUTPUT (one row per district)
 *   unit_id, country, admin_level, admin_unit : district identification
 *   n_matches        : number of GAUL polygons matching the name; must be 1
 *                      (0 = name not found, 2 or more = ambiguous name)
 *   flood_max_events : maximum number of flood events in any 250 m land pixel
 *   flood            : flood-prone, yes / no
 *   urban_23_30      : % of land area in GHS-SMOD classes 23 + 30
 *   setting          : urban / rural
 *   wild_dist_km     : distance to the nearest wildlife protected area (km)
 *   national_park    : near national park (wildlife protected area), yes / no
 *   border_dist_km   : distance to the nearest international land boundary (km)
 *   border           : near international border, yes / no
 *   Distances are measured from the district boundary; 0 means the district
 *   intersects or touches the feature, and 999 means none was found within
 *   250 km (protected areas) or 800 km (borders).
 ******************************************************************************/

// ---- classification thresholds ---------------------------------------------
var FLOOD_MIN = 1;    // flood = yes if flood_max_events >= FLOOD_MIN
var URBAN_MIN = 10;   // setting = urban if urban_23_30 > URBAN_MIN (%)
var PARK_KM   = 6;    // national_park = yes if wild_dist_km <= PARK_KM (km)
var BORDER_KM = 0;    // border = yes if border_dist_km <= BORDER_KM (km)
var SIZE_MIN  = 50;   // minimum size of a wildlife protected area (km2)

// ---- administrative units (GAUL names) ------------------------------------
var UNITS = [
  {unit_id:"U01", country:"Ethiopia", gaul_country:"Ethiopia", level:2, name:"Nuer"},
  {unit_id:"U02", country:"Ethiopia", gaul_country:"Ethiopia", level:2, name:"Shabelle"},
  {unit_id:"U03", country:"Ethiopia", gaul_country:"Ethiopia", level:2, name:"South Omo"},
  {unit_id:"U04", country:"Ethiopia", gaul_country:"Ethiopia", level:2, name:"Zone 3"},
  {unit_id:"U05", country:"Kenya", gaul_country:"Kenya", level:2, name:"Baringo"},
  {unit_id:"U06", country:"Kenya", gaul_country:"Kenya", level:2, name:"Busia"},
  {unit_id:"U07", country:"Kenya", gaul_country:"Kenya", level:2, name:"Garissa"},
  {unit_id:"U08", country:"Kenya", gaul_country:"Kenya", level:2, name:"Isiolo"},
  {unit_id:"U09", country:"Kenya", gaul_country:"Kenya", level:2, name:"Kiambu"},
  {unit_id:"U10", country:"Kenya", gaul_country:"Kenya", level:2, name:"Laikipia"},
  {unit_id:"U11", country:"Kenya", gaul_country:"Kenya", level:2, name:"Machakos"},
  {unit_id:"U12", country:"Kenya", gaul_country:"Kenya", level:2, name:"Makueni"},
  {unit_id:"U13", country:"Kenya", gaul_country:"Kenya", level:2, name:"Migori"},
  {unit_id:"U14", country:"Kenya", gaul_country:"Kenya", level:2, name:"Mombasa"},
  {unit_id:"U15", country:"Kenya", gaul_country:"Kenya", level:2, name:"Nakuru"},
  {unit_id:"U16", country:"Kenya", gaul_country:"Kenya", level:2, name:"Nyeri"},
  {unit_id:"U17", country:"Kenya", gaul_country:"Kenya", level:2, name:"Taita Taveta"},
  {unit_id:"U18", country:"Kenya", gaul_country:"Kenya", level:2, name:"Tana River"},
  {unit_id:"U19", country:"Kenya", gaul_country:"Kenya", level:2, name:"Trans Nzoia"},
  {unit_id:"U20", country:"Madagascar", gaul_country:"Madagascar", level:2, name:"Mananjary"},
  {unit_id:"U21", country:"Rwanda", gaul_country:"Rwanda", level:2, name:"Bugesera"},
  {unit_id:"U22", country:"Rwanda", gaul_country:"Rwanda", level:2, name:"Gakenke"},
  {unit_id:"U23", country:"Rwanda", gaul_country:"Rwanda", level:2, name:"Kamonyi"},
  {unit_id:"U24", country:"Rwanda", gaul_country:"Rwanda", level:2, name:"Kirehe"},
  {unit_id:"U25", country:"Rwanda", gaul_country:"Rwanda", level:2, name:"Ngoma"},
  {unit_id:"U26", country:"Rwanda", gaul_country:"Rwanda", level:2, name:"Nyagatare"},
  {unit_id:"U27", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Bahi"},
  {unit_id:"U28", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Ikungi"},
  {unit_id:"U29", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Ilala"},
  {unit_id:"U30", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Kibondo"},
  {unit_id:"U31", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Kyela"},
  {unit_id:"U32", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Momba"},
  {unit_id:"U33", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Morogoro"},
  {unit_id:"U34", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Rufiji"},
  {unit_id:"U35", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Sengerema"},
  {unit_id:"U36", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Simanjiro"},
  {unit_id:"U37", country:"Tanzania", gaul_country:"United Republic of Tanzania", level:2, name:"Tunduma"},
  {unit_id:"U38", country:"Uganda", gaul_country:"Uganda", level:1, name:"Hoima"},
  {unit_id:"U39", country:"Uganda", gaul_country:"Uganda", level:1, name:"Kabale"},
  {unit_id:"U40", country:"Uganda", gaul_country:"Uganda", level:1, name:"Kanungu"},
  {unit_id:"U41", country:"Uganda", gaul_country:"Uganda", level:1, name:"Kasese"},
  {unit_id:"U42", country:"Uganda", gaul_country:"Uganda", level:1, name:"Kibaale"},
  {unit_id:"U43", country:"Uganda", gaul_country:"Uganda", level:1, name:"Kisoro"},
  {unit_id:"U44", country:"Uganda", gaul_country:"Uganda", level:1, name:"Masindi"},
  {unit_id:"U45", country:"Uganda", gaul_country:"Uganda", level:1, name:"Rubirizi"},
  {unit_id:"U46", country:"Uganda", gaul_country:"Uganda", level:1, name:"Ssembabule"},
  {unit_id:"U47", country:"Zambia", gaul_country:"Zambia", level:2, name:"Chadiza"},
  {unit_id:"U48", country:"Zambia", gaul_country:"Zambia", level:2, name:"Chavuma"},
  {unit_id:"U49", country:"Zambia", gaul_country:"Zambia", level:2, name:"Chibombo"},
  {unit_id:"U50", country:"Zambia", gaul_country:"Zambia", level:2, name:"Chililabombwe"},
  {unit_id:"U51", country:"Zambia", gaul_country:"Zambia", level:2, name:"Kalomo"},
  {unit_id:"U52", country:"Zambia", gaul_country:"Zambia", level:2, name:"Livingstone"},
  {unit_id:"U53", country:"Zambia", gaul_country:"Zambia", level:2, name:"Mazabuka"},
  {unit_id:"U54", country:"Zambia", gaul_country:"Zambia", level:2, name:"Mbala"},
  {unit_id:"U55", country:"Zambia", gaul_country:"Zambia", level:2, name:"Monze"},
  {unit_id:"U56", country:"Zambia", gaul_country:"Zambia", level:2, name:"Mumbwa"},
  {unit_id:"U57", country:"Zambia", gaul_country:"Zambia", level:2, name:"Mungwi"},
  {unit_id:"U58", country:"Zambia", gaul_country:"Zambia", level:2, name:"Namwala"},
  {unit_id:"U59", country:"Zambia", gaul_country:"Zambia", level:2, name:"Siavonga"},
  {unit_id:"U60", country:"Zimbabwe", gaul_country:"Zimbabwe", level:2, name:"Chiredzi"},
];

// ---- datasets ----------------------------------------------------------------
var ME = 1000, PA_SEARCH = 250000, BRD_SEARCH = 800000;
var REGION = ee.Geometry.Rectangle([20, -27, 52, 16]);
var gaul1 = ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level1').filterBounds(REGION);
var gaul2 = ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level2').filterBounds(REGION);
var world = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017');
var px = ee.Image.pixelArea();

// flooding: Global Flood Database, permanent water masked with JRC GSW
var gfd  = ee.ImageCollection('GLOBAL_FLOOD_DB/MODIS_EVENTS/V1');
var gsw  = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0);
var land = gsw.lt(50);
var freqLand = gfd.select('flooded').sum().unmask(0).multiply(land).rename('freq');

// settlement: GHS-SMOD 2020, land area in all classes and in classes 23 + 30
var smodCol = ee.ImageCollection('JRC/GHSL/P2023A/GHS_SMOD');
var spick = smodCol.filterDate('2020-01-01', '2020-12-31');
spick = ee.ImageCollection(ee.Algorithms.If(spick.size().gt(0), spick,
          smodCol.filter(ee.Filter.stringContains('system:index', '2020'))));
var smod = ee.Image(spick.first());
smod = smod.select([ee.String(smod.bandNames().get(0))]).rename('smod');
var smodImg = smod.gte(11).multiply(px).rename('s_land')
  .addBands(smod.eq(23).or(smod.eq(30)).multiply(px).rename('s_urban'));

// wildlife protected areas: WDPA
var wdpa = ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(REGION);
var paBase = wdpa.filter(ee.Filter.neq('MARINE', '2'))
  .filter(ee.Filter.inList('STATUS', ['Designated', 'Inscribed', 'Established', 'Adopted']));
var wildDesig = ee.Filter.or(
  ee.Filter.stringContains('DESIG_ENG', 'National Park'), ee.Filter.stringContains('DESIG_ENG', 'Game Reserve'),
  ee.Filter.stringContains('DESIG_ENG', 'National Reserve'), ee.Filter.stringContains('DESIG_ENG', 'Wildlife'),
  ee.Filter.stringContains('DESIG_ENG', 'Conservancy'), ee.Filter.stringContains('DESIG_ENG', 'Conservation Area'),
  ee.Filter.stringContains('DESIG_ENG', 'Game Controlled'), ee.Filter.stringContains('DESIG_ENG', 'Ramsar'));
var wildSet = paBase.filter(ee.Filter.or(wildDesig, ee.Filter.inList('IUCN_CAT', ['II', 'IV', 'VI'])))
  .filter(ee.Filter.or(ee.Filter.gte('REP_AREA', SIZE_MIN), ee.Filter.gte('GIS_AREA', SIZE_MIN)));

// ---- helpers -----------------------------------------------------------------
function r3(x) { return ee.Number(x).multiply(1000).round().divide(1000); }
function gnum(d, k) { return ee.Number(ee.Algorithms.If(d.get(k), d.get(k), 0)); }
function yesno(cond) { return ee.String(ee.Algorithms.If(cond, 'yes', 'no')); }
// distance (km) from a polygon to the nearest feature of a collection
function nearestKm(geom, fc, search) {
  var near = fc.filterBounds(geom.buffer(search));
  return ee.Number(ee.Algorithms.If(near.size().gt(0),
    ee.Number(near.iterate(function (ft, prev) {
      return ee.Number(prev).min(geom.distance(ee.Feature(ft).geometry(), ME));
    }, 1e13)).divide(1000), 999));
}

// ---- match units to GAUL polygons by name -------------------------------------
var unitsFC = ee.FeatureCollection(UNITS.map(function (u) {
  return ee.Feature(null, {unit_id: u.unit_id, country: u.country, admin_level: 'GAUL' + u.level,
    admin_unit: u.name, key: u.gaul_country + '|' + u.level + '|' + u.name});
}));
var gaulKeyed = gaul1.map(function (f) {
    return f.set('key', ee.String(f.get('ADM0_NAME')).cat('|1|').cat(f.get('ADM1_NAME')));
  }).merge(gaul2.map(function (f) {
    return f.set('key', ee.String(f.get('ADM0_NAME')).cat('|2|').cat(f.get('ADM2_NAME')));
  }));
var byKey = ee.Filter.equals({leftField: 'key', rightField: 'key'});
var matched   = ee.Join.saveAll({matchesKey: 'polys'}).apply(unitsFC, gaulKeyed, byKey);
var unmatched = ee.Join.inverted().apply(unitsFC, gaulKeyed, byKey).map(function (f) {
  return f.set('n_matches', 0);
});

// ---- per unit ----------------------------------------------------------------
var computed = matched.map(function (f) {
  var polys = ee.List(f.get('polys'));
  var poly = ee.Feature(polys.get(0)).geometry();

  // flood: maximum number of flood events in any land pixel
  var mxD = freqLand.reduceRegion({reducer: ee.Reducer.max(), geometry: poly, scale: 250,
    maxPixels: 1e10, bestEffort: true, tileScale: 4}).get('freq');
  var mx = ee.Number(ee.Algorithms.If(mxD, mxD, 0));

  // urban: % of land area in SMOD classes 23 + 30
  var su = smodImg.reduceRegion({reducer: ee.Reducer.sum(), geometry: poly, scale: 1000,
    maxPixels: 1e10, bestEffort: true, tileScale: 4});
  var urbanPct = gnum(su, 's_urban').divide(gnum(su, 's_land').max(1)).multiply(100);

  // wildlife protected area: distance from the district boundary
  var wildKm = nearestKm(poly, wildSet, PA_SEARCH);

  // international border: distance to the nearest other country
  var own = world.filterBounds(poly.centroid(100));
  var ownName = ee.Algorithms.If(own.size().gt(0), ee.Feature(own.first()).get('country_na'), 'NA');
  var others = world.filterBounds(poly.buffer(BRD_SEARCH)).filter(ee.Filter.neq('country_na', ownName));
  var borderKm = ee.Number(ee.Algorithms.If(others.size().gt(0),
    ee.Number(others.iterate(function (ft, prev) {
      return ee.Number(prev).min(poly.distance(ee.Feature(ft).geometry(), ME));
    }, 1e13)).divide(1000), 999));

  return ee.Feature(null, {
    unit_id: f.get('unit_id'), country: f.get('country'), admin_level: f.get('admin_level'),
    admin_unit: f.get('admin_unit'), n_matches: polys.size(),
    flood_max_events: mx, flood: yesno(mx.gte(FLOOD_MIN)),
    urban_23_30: r3(urbanPct), setting: ee.String(ee.Algorithms.If(urbanPct.gt(URBAN_MIN), 'urban', 'rural')),
    wild_dist_km: r3(wildKm), national_park: yesno(wildKm.lte(PARK_KM)),
    border_dist_km: r3(borderKm), border: yesno(borderKm.lte(BORDER_KM))
  });
});
var results = computed.merge(unmatched).sort('unit_id');

// ---- checks -------------------------------------------------------------------
print('Units in list:', UNITS.length);
print('Units NOT matched to exactly one polygon (fix these names):',
      results.filter(ee.Filter.neq('n_matches', 1)));
print('First 5 results:', results.limit(5));

// ---- export -------------------------------------------------------------------
Export.table.toDrive({
  collection: results, description: 'RVF_eco_v11_units', fileFormat: 'CSV',
  selectors: ['unit_id', 'country', 'admin_level', 'admin_unit', 'n_matches',
              'flood_max_events', 'flood', 'urban_23_30', 'setting',
              'wild_dist_km', 'national_park', 'border_dist_km', 'border']
});
