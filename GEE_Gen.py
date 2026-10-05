import ee
try:
    ee.Initialize(project='apt-trainer-483518-q8')
except Exception:
    ee.Authenticate()
    ee.Initialize(project='apt-trainer-483518-q8')
from functools import reduce

from datetime import datetime
from dateutil.relativedelta import relativedelta   # add these two to the imports at the top

########################### 
#Imports
########################### 

LFMC = ee.FeatureCollection("projects/apt-trainer-483518-q8/assets/LFMC_WUS")
SAR1C = ee.ImageCollection("COPERNICUS/S1_GRD")
STATE = ee.FeatureCollection("TIGER/2018/States")
S2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
DAYMT = ee.ImageCollection("NASA/ORNL/DAYMET_V4")
GRDMT = ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")
NLCD = ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/NLCD")
DW = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
PURE = ee.FeatureCollection("projects/apt-trainer-483518-q8/assets/PURESITES")


########################### 
# Date and Location 
########################### 

WUS = STATE.filter(ee.Filter.inList('NAME', [
    'Washington', 'Oregon', 'California', 'Nevada', 'Idaho', 'Montana',
    'Utah', 'Arizona', 'Wyoming', 'Colorado', 'New Mexico', 'Texas'
]))

bound = WUS

#puresites = PURE.aggregate_array('﻿Puresites') #Invisible character, be aware

########################### 
# Collections
########################### 

def MaskS2(image):
    QA60 = image.select('QA60')
    return image.updateMask(QA60.bitwiseAnd(1 << 10).eq(0)
                            .And(QA60.bitwiseAnd(1 << 11).eq(0)))

COLLECTIONS = [
    {
        'name': 'MSI',  # Sentinel 2
        'collection': S2,
        'bands': ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'],
        'newnames': ['B', 'G', 'R', 'NIR', 'SW1', 'SW2'],
        'scale': 10,
        'maskFn': MaskS2,
    },
    {
        'name': 'SAR1C',  # Sentinel 1
        'collection': SAR1C,
        'bands': ['VV', 'VH', 'angle'],
        'newnames': ['VV', 'VH', 'angle'],
        'scale': 10,
        'maskFn': None,
    },
]
########################### 
# Sites
########################### 


def add_time_and_buffer(f):
    """Adding the ee date and the buffer around points where samples are taken"""
    date_parsed = ee.Date.parse('yyyy-MM-dd', ee.String(f.get('DATE')))
    return f.set('system:time_start', date_parsed.millis()).buffer(80)


def add_lat_lon(f):
    # lat/lon can help match sites with static parameters later
    coords = f.geometry().centroid(1).coordinates()
    return f.set({'lon': coords.get(0), 'lat': coords.get(1)})

def add_year_month(f):
    d = ee.Date(f.get('system:time_start'))
    return f.set({'year': d.get('year'), 'month': d.get('month')})




########################### 
# Match and Reduce
########################### 


TEMPORAL_SPATIAL_FILTER = ee.Filter.And(
    ee.Filter.maxDifference(
        difference= 10 * 24 * 60 * 60 * 1000,
        leftField=  'system:time_start',
        rightField= 'system:time_start'
  ),
  ee.Filter.intersects(
        leftField=  '.geo',
        rightField= '.geo'
  )
) 

def match_and_reduce(sites, image_col, scale):
    def imgfind(feature):
        closest_image = ee.Image(feature.get('closestImage'))
        medians = closest_image.reduceRegion(
            reducer=ee.Reducer.median(),
            geometry=feature.geometry(),
            scale=scale,
            maxPixels=1e9,
            bestEffort=True
        )
        new_keys = medians.keys().map(lambda k: ee.String(k).cat('_0'))
        renamed = ee.Dictionary.fromLists(new_keys, medians.values())
        return (feature.set(renamed)
                .set('matched_img_date', closest_image.date().format('YYYY-MM-dd')))

    join = ee.Join.saveBest(
        matchKey='closestImage',
        measureKey='timeDiff'
    )
    # back to using the shared filter
    joined = join.apply(sites, image_col, TEMPORAL_SPATIAL_FILTER)
    return joined.map(imgfind)
    
def matchAndReduce_at(sites, imageCol, scale, days, suffix):
    def imgfind(feature):
        img = ee.Image(feature.get('closestImage'))
        values = img.reduceRegion(
            reducer=    ee.Reducer.median(),
            geometry=   feature.geometry(),
            scale=      scale,
            maxPixels=  1e9,
            bestEffort= True
        )
        # renaming each band with the suffix
        new_keys = values.keys().map(lambda k: ee.String(k).cat(suffix))
        renamed = ee.Dictionary.fromLists(new_keys, values.values())
        return (feature.set(renamed)
                .set('dateDiff' + suffix, feature.get('dateDiff')))
            
    FILTER = ee.Filter.And(
        ee.Filter.maxDifference(
            difference= days * 24 * 60 * 60 * 1000,
            leftField=  'system:time_start',
            rightField= 'system:time_start'
        ),
        ee.Filter.greaterThanOrEquals(
            leftField=  'system:time_start',
            rightField= 'system:time_start'
        ),
        ee.Filter.intersects(
            leftField= '.geo',
            rightField= '.geo'
        )
    )
    
    join = ee.Join.saveBest(
        matchKey= 'closestImage',
        measureKey= 'dateDiff'
    )
    
    joined = join.apply(sites,imageCol, FILTER)
    return joined.map(imgfind)


def matchAndReduce30(sites, imageCol, scale):
    return matchAndReduce_at(sites, imageCol, scale, 30, '_30')

def matchAndReduce60(sites, imageCol, scale):
    return matchAndReduce_at(sites, imageCol, scale, 60, '_60')

def matchAndReduce90(sites, imageCol, scale):
    return matchAndReduce_at(sites, imageCol, scale, 90, '_90')

def matchAndReduce120(sites, imageCol, scale):
    return matchAndReduce_at(sites, imageCol, scale, 120, '_120')



########################### 
# Automatic Import (work in progress)
###########################                                    



########################### 
# Preprocess
###########################    



def preprocess(config, image_start, end_date):
    filtered = (config['collection']
                .filterBounds(bound)
                .filterDate(image_start,end_date))

    if config['name'] == 'SAR1C':
        filtered = (filtered
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')))

    mask_fn = config['maskFn']
    if mask_fn:
        filtered = filtered.map(
            lambda img: mask_fn(img).set('system:time_start', img.get('system:time_start')))

    if config['bands']:
        if config['newnames']:
            filtered = filtered.select(config['bands'], config['newnames'])
        else:
            filtered = filtered.select(config['bands'])

    return filtered

########################### 
# FMport
###########################    

def FMport(start_date, end_date, pf_type, folder='GEE_Py_Data'):
    image_start = ee.Date(start_date).advance(-120, 'day')
    start_year, end_year = start_date[:4], end_date[:4]    
    sites = (LFMC
        .filterBounds(bound)
        .map(add_time_and_buffer)
        .filterDate(start_date, end_date)
        #.filter(ee.Filter.inList('SITENAME', puresites))
        .filter(ee.Filter.eq('FUNCTIONALTYPE', pf_type)))
    sites = sites.map(add_lat_lon).map(add_year_month)

    def combine_each(current, config):
        processed = preprocess(config, image_start, end_date)
        output0 = match_and_reduce(current, processed, config['scale'])
        output30 = matchAndReduce30(output0, processed, config['scale'])
        output60 = matchAndReduce60(output30, processed, config['scale'])
        # output90 = matchAndReduce90(output60,processed, config['scale'])
        # output120 = matchAndReduce120(output90,processed, config['scale'])
        # i commented the last 2 out for simplicity but i can uncomment if i decide i want to include them 
        return output60
        
    combined_sites = reduce(combine_each, COLLECTIONS, sites)

    task = ee.batch.Export.table.toDrive(
        collection=combined_sites,
        description=f'WUS{end_year}{pf_type}_export',
        folder=folder,
        fileNamePrefix=f'{start_date}_{end_date}_{pf_type.upper()}',
        fileFormat='CSV'
    )
    task.start()
    return task

########################### 
# FMport batched
###########################   

def FMbatch(start_date, end_date, pf_type, folder='GEE_Py_Data', chunk_months=12):
    fmt = '%Y-%m-%d'
    window_start = datetime.strptime(start_date, fmt)
    final_end = datetime.strptime(end_date, fmt)

    tasks = []
    while window_start < final_end:
        window_end = min(window_start + relativedelta(months=chunk_months), final_end)
        task = FMport(window_start.strftime(fmt), window_end.strftime(fmt), pf_type, folder)
        tasks.append(task)
        window_start = window_end
    return tasks







