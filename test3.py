from math import floor, radians, log, cos, pi, tan, sqrt, atan, exp, degrees
import geopandas as gpd
from shapely.geometry import Polygon
from osgeo import gdal
from collections import Counter
import numpy as np
import csv
import gzip

import time

import sys

# first argument of script is the zoom level.
# second argument of script is the filename.
#
# Example z=19 and filename=clipped-esa.tif:
#
# python3 test3.py 19 clipped-esa.tif
#

z = int(sys.argv[1])
print('z =', z)

filename = sys.argv[2]
print('filename =', filename)

# map from pixel value to index for summing, color, and name
classes = {
    10: { 'index': 0, 'color': '0,100,0,255', 'name': 'Tree cover' },
    20: { 'index': 1, 'color': '255,187,34,255', 'name': 'Shrubland' },
    30: { 'index': 2, 'color': '255,255,76,255', 'name': 'Grassland' },
    40: { 'index': 3, 'color': '240,150,255,255', 'name': 'Cropland' },
    50: { 'index': 4, 'color': '250,0,0,255', 'name': 'Built-up' },
    60: { 'index': 5, 'color': '180,180,180,255', 'name': 'Bare / sparse vegetation' },
    70: { 'index': 6, 'color': '240,240,240,255', 'name': 'Snow and Ice' },
    80: { 'index': 7, 'color': '0,100,200,255', 'name': 'Permanent water bodies' },
    90: { 'index': 8, 'color': '0,150,160,255', 'name': 'Herbaceous wetland' },
    95: { 'index': 9, 'color': '0,207,117,255', 'name': 'Mangroves' },
    100: { 'index': 10, 'color': '250,230,160,255', 'name': 'Moss and lichen' },
}

classes_index_to_value = {
}
for c in classes:
    classes_index_to_value[classes[c]['index']] = c

# for c in classes:
#     print(f'{c},')
#     print(f'"rgba({classes[c]["color"]})"')

def get_x_max(z):
    return 2 ** z

def get_triangle(x, y, z):
    x_max = get_x_max(z)
    i = floor(x)
    j = floor(y)

    x_triangle = None
    y_triangle = None

    if j % 2 == 0: # j is even
        if i % 2 == 1: # i is odd
            if (x - i) > (y - j):
                x_triangle = i + 1
            else:
                x_triangle = i
        else: # i is even
            if (x - i) < 1 - (y - j):
                x_triangle = i
            else:
                x_triangle = i + 1
    else: # j is odd
        if i % 2 == 1: # i is odd
            if (x - i) < 1 - (y - j):
                x_triangle = i
            else:
                x_triangle = i + 1
        else: # i is even
            if (x - i) > (y - j):
                x_triangle = i + 1
            else:
                x_triangle = i

    x_triangle %= x_max
    y_triangle = j

    return x_triangle, y_triangle

def get_x_y(lat, lon, z):
    x_max = get_x_max(z)

    latRad = radians(lat)
    lonRad = radians(lon)
    
    x_norm = (lonRad + pi) / (2 * pi) # between 0 and 1
    y_norm = (1 - log(tan(latRad) + 1 / cos(latRad)) / pi) / 2 # between 0 and 1 for inputs between ~-85 deg and ~85 deg

    x = x_norm * x_max
    y = y_norm * x_max / sqrt(3)

    return x, y


def get_lat_lon(x, y, z):
    x_max = get_x_max(z)
    x_norm = x / x_max
    y_norm = y / (x_max / sqrt(3))
    
    lonRad = (x_norm * 2 * pi) - pi
    latRad = (2 * atan(exp((1 - 2 * y_norm) * pi)) - pi / 2)
    
    lat = degrees(latRad)
    lon = degrees(lonRad)
    
    return lat, lon

def get_vertices(x_triangle, y_triangle):
    i = x_triangle
    j = y_triangle

    p0 = None
    p1 = None
    p2 = None

    if j % 2 == 1: # j is odd
        if i % 2 == 0: # i is even
            p0 = (i, j)
            p1 = (i+1, j+1)
            p2 = (i-1, j+1)
        else: # i is odd
            p0 = (i, j+1)
            p1 = (i-1, j)
            p2 = (i+1, j)
    else:
        if i % 2 == 0: # i is even
            p0 = (i, j+1)
            p1 = (i-1, j)
            p2 = (i+1, j)
        else: # i is odd
            p0 = (i, j)
            p1 = (i+1, j+1)
            p2 = (i-1, j+1)
    
    return [p0, p1, p2]
        
def get_polygon(x_triangle, y_triangle, z):
    vertices = get_vertices(x_triangle, y_triangle)
    vertices_lon_lat = []
    for vertex in vertices:
        lat, lon = get_lat_lon(vertex[0], vertex[1], z)
        vertices_lon_lat.append([lon, lat])
    return Polygon(vertices_lon_lat)


def get_hexagon(x_triangle, y_triangle):
    # (x_hexagon, y_hexagon) is the (x_triangle, y_triangle) of the top left triangle of the hexagon
    x_hexagon = None
    y_hexagon = None
    if x_triangle % 6 in [1, 2, 3]:
        if y_triangle % 2 == 0:
            x_hexagon, y_hexagon = (int(x_triangle / 6) * 6 + 1, y_triangle)
        else:
            x_hexagon, y_hexagon = (int(x_triangle / 6) * 6 + 1, y_triangle - 1)
    if (x_triangle - 3) % 6 in [1, 2, 3]:
        if y_triangle % 2 == 1:
            x_hexagon, y_hexagon = (int((x_triangle - 3) / 6) * 6 + 4, y_triangle)
        else:
            x_hexagon, y_hexagon = (int((x_triangle - 3) / 6) * 6 + 4, y_triangle - 1)
    return x_hexagon, y_hexagon

def get_triangles(x_hexagon, y_hexagon):
    result = []
    for dx in range(3):
        for dy in range(2):
            result.append((x_hexagon + dx, y_hexagon + dy))
    return result

# test get_triangle
# z = 2
# x = 3.1
# y = 1.5
# print(get_triangle(x, y, z))

# # test get_x_y
# lat = 0
# lon = 0
# z = 2
# x, y = get_x_y(lat, lon, z)
# print(get_triangle(x, y, z))

# # test get_lat_lon
# z = 2
# x = 4
# y = 4
# print(get_lat_lon(x, y, z))


# # test get_vertices
# x_triangle = 0
# y_triangle = 0
# x_max = 128


# data = {
#     'x_triangle': [],
#     'y_triangle': [],
#     'value': [],
#     'geometry': []
# }
# for x_triangle in range(128):
#     for y_triangle in range(64):
#         data['x_triangle'].append(x_triangle)
#         data['y_triangle'].append(y_triangle)
#         value = ''
#         if x_triangle % 6 in [1, 2, 3]:
#             if y_triangle % 2 == 0:
#                 value = f'({int(x_triangle / 6) * 6 + 1},{y_triangle})'
#             else:
#                 value = f'({int(x_triangle / 6) * 6 + 1},{y_triangle - 1})'
#         if (x_triangle - 3) % 6 in [1, 2, 3]:
#             if y_triangle % 2 == 1:
#                 value = f'({int((x_triangle - 3) / 6) * 6 + 4},{y_triangle})'
#             else:
#                 value = f'({int((x_triangle - 3) / 6) * 6 + 4},{y_triangle - 1})'
#         data['value'].append(value)
#         data['geometry'].append(get_polygon(x_triangle, y_triangle, x_max))

# gdf = gpd.GeoDataFrame(data, crs='epsg:4326')
# gdf = gdf.dissolve('value')
# print(gdf)
# gdf.to_file(filename='polygon.gpkg', driver="GPKG")




# # test get_vertices
# x_triangle = 0
# y_triangle = 0
# x_max = 6 * 20


# data = {
#     'value': [],
#     'geometry': []
# }
# for i in range(1, 21, 6):
#     for j in range(0, 10, 2):
#         value = f'({i},{j})'
#         y_triangle = j
#         for x_triangle in range(i, i + 3):
#             data['value'].append(value)
#             data['geometry'].append(get_polygon(x_triangle, y_triangle, x_max))
#             data['value'].append(value)
#             data['geometry'].append(get_polygon(x_triangle, y_triangle + 1, x_max))

# for i in range(1, 21, 6):
#     for j in range(0, 10, 2):
#         value = f'({i},{j})'
#         y_triangle = j
#         for x_triangle in range(i, i + 3):
#             data['value'].append(value)
#             data['geometry'].append(get_polygon(x_triangle, y_triangle, x_max))
#             data['value'].append(value)
#             data['geometry'].append(get_polygon(x_triangle, y_triangle + 1, x_max))

# gdf = gpd.GeoDataFrame(data, crs='epsg:4326')
# gdf = gdf.dissolve('value')
# print(gdf)
# gdf.to_file(filename='polygon.gpkg', driver="GPKG")




def process_image(filename, z):

    # Open the geotiff file
    dataset = gdal.Open(filename)
    
    # Get the pixel values and the geotransform information
    band = dataset.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    geotransform = dataset.GetGeoTransform()
    pixel_values = band.ReadAsArray().astype(float)
    
    hexagons = {}

    tic = time.time()
    # Loop through each pixel and calculate the covered area
    # print(pixel_values.shape[0])

    # ## NumPy variant
    # nodata_mask = (pixel_values != nodata)

    # # Create arrays for j and i values
    # j_values, i_values = np.meshgrid(np.arange(pixel_values.shape[1]), np.arange(pixel_values.shape[0]))

    # # Calculate lon and lat for valid pixels in a vectorized manner
    # lon = geotransform[0] + (j_values[nodata_mask] + 0.5) * geotransform[1] + (i_values[nodata_mask] + 0.5) * geotransform[2]
    # lat = geotransform[3] + (j_values[nodata_mask] + 0.5) * geotransform[4] + (i_values[nodata_mask] + 0.5) * geotransform[5]
    # x, y = get_x_y_np(lat, lon, z)
    # x_triangle, y_triangle = get_triangle_np(x, y, z)
    # x_hexagon, y_hexagon = get_hexagon_np(x_triangle, y_triangle)

    for i in range(pixel_values.shape[0]):
        for j in range(pixel_values.shape[1]):
            # Skip nodata pixels
            if pixel_values[i][j] == nodata:
                continue
                        
            # Convert the pixel coordinates to map coordinates
            lon = geotransform[0] + (j + 0.5) * geotransform[1] + (i + 0.5) * geotransform[2]
            lat = geotransform[3] + (j + 0.5) * geotransform[4] + (i + 0.5) * geotransform[5]

            x, y = get_x_y(lat, lon, z)
            x_triangle, y_triangle = get_triangle(x, y, z)
            x_hexagon, y_hexagon = get_hexagon(x_triangle, y_triangle)
            
            if (x_hexagon, y_hexagon) in hexagons:
                if int(pixel_values[i][j]) in hexagons[(x_hexagon, y_hexagon)]:
                    hexagons[(x_hexagon, y_hexagon)][int(pixel_values[i][j])] += 1
                else:
                    hexagons[(x_hexagon, y_hexagon)][int(pixel_values[i][j])] = 1
            else:
                hexagons[(x_hexagon, y_hexagon)] = {
                    int(pixel_values[i][j]): 1
                }

            # Print the progress
            pixels_processed = (i * pixel_values.shape[1]) + j + 1
            total_pixels = pixel_values.shape[0] * pixel_values.shape[1]
            if pixels_processed % 10000 == 0:
                progress = pixels_processed / total_pixels
                print(f"{filename}: {progress:.2%} complete", end="\r")
    
    print('total time', time.time() - tic)
    hexagons_majority = {}
    for key in hexagons:
        max_index = None
        max_value = 0
        for landcover_index in hexagons[key]:
            if hexagons[key][landcover_index] > max_value:
                max_index = landcover_index
                max_value = hexagons[key][landcover_index]
        hexagons_majority[key] = max_index
    
    # print(hexagons_majority)

    data = {
        'value': [],
        'geometry': []
    }
    for key in hexagons_majority:
        triangles = get_triangles(key[0], key[1])
        for triangle in triangles:
            data['value'].append(hexagons_majority[key])
            data['geometry'].append(get_polygon(triangle[0], triangle[1], z))
    
    gdf = gpd.GeoDataFrame(data, crs='epsg:4326')
    gdf = gdf.dissolve('value')
    print(gdf)
    gdf.to_file(filename='polygon.gpkg', driver="GPKG")
            

def step_1(filename, z):

    # compute the x_triangle, y_triangle values for each pixel and write them to a csv file with the pixel value
    # produces a csv file which looks like this:
    # 
    # x_triangle, y_triangle, value
    # 51886, 19110, 70
    # 51886, 19110, 70
    # ...

    tic = time.time()

    # Open the geotiff file
    dataset = gdal.Open(f'downloads/{filename}')
    
    # Get the pixel values and the geotransform information
    band = dataset.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    geotransform = dataset.GetGeoTransform()
    # pixel_values = band.ReadAsArray().astype(float)
    pixel_values = band.ReadAsArray().astype(np.uint8)
    
    longitudes = [geotransform[0] + (j + 0.5) * geotransform[1] + (0 + 0.5) * geotransform[2] for j in range(pixel_values.shape[1])]
    latitudes = [geotransform[3] + (0 + 0.5) * geotransform[4] + (i + 0.5) * geotransform[5] for i in range(pixel_values.shape[0])]

    # print('longitudes', longitudes[0], longitudes[-1])
    # print('latitudes', latitudes[0], latitudes[-1])

    xs = []
    for longitude in longitudes:
        x, _ = get_x_y(lat=0, lon=longitude, z=z)
        xs.append(x)
    
    ys = []
    for latitude in latitudes:
        _, y = get_x_y(lat=latitude, lon=0, z=z)
        ys.append(y)
    
    # print('xs', xs[0], xs[-1])
    # print('ys', ys[0], ys[-1])

    # with gzip.open(f'output_step_1/{filename}-triangles-z{z}.csv.gz', 'wt') as f:
    with open(f'output_step_1/{filename}-triangles-z{z}.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['x_triangle', 'y_triangle', 'value'])
            
        for i in range(pixel_values.shape[0]):
            for j in range(pixel_values.shape[1]):
                # Skip nodata pixels
                if pixel_values[i][j] == nodata:
                    continue

                x, y = xs[j], ys[i]
                x_triangle, y_triangle = get_triangle(x, y, z)

                writer.writerow([x_triangle, y_triangle, int(pixel_values[i][j])])

                # Print the progress
                pixels_processed = (i * pixel_values.shape[1]) + j + 1
                total_pixels = pixel_values.shape[0] * pixel_values.shape[1]
                if pixels_processed % 1000000 == 0:
                    progress = pixels_processed / total_pixels
                    print(f"step_1 {progress:.2%} complete")

    # print(x_triangles)
    print('step_1', time.time() - tic)


def step_2(filename, z):

    # sum the triangle values from step 1 and creates a triangle counts csv file
    # where triangles only appear once
    # example csv file:
    # 
    # x_triangle,y_triangle,count_index_0,count_index_1,count_index_2,count_index_3,count_index_4,count_index_5,count_index_6,count_index_7,count_index_8,count_index_9,count_index_10
    # 51849,19102,0,0,0,0,0,0,0,501,0,0,0
    # 51850,19102,0,0,0,0,0,0,0,386,0,0,0
    # 51851,19102,0,0,0,0,0,0,0,667,0,0,0
    # 51852,19102,0,0,0,0,0,0,0,387,0,0,0
    # ...

    tic = time.time()

    triangles = {}

    # with gzip.open(f'output_step_1/{filename}-triangles-z{z}.csv.gz', 'rt') as f:
    with open(f'output_step_1/{filename}-triangles-z{z}.csv', 'r') as f:
        reader = csv.reader(f)
        i = 0
        next(reader) # skip header line
        for row in reader:
            x_triangle, y_triangle, value = [int(element) for element in row]
            
            if (x_triangle, y_triangle) not in triangles:            
                triangles[(x_triangle, y_triangle)] = [0 for _ in range(11)]
            
            index = classes[value]['index']
            triangles[(x_triangle, y_triangle)][index] += 1

    # with gzip.open(f'output_step_2/{filename}-triangle-counts-z{z}.csv.gz', 'wt') as f:
    with open(f'output_step_2/{filename}-triangle-counts-z{z}.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow([
            'x_triangle',
            'y_triangle',
            'count_index_0',
            'count_index_1',
            'count_index_2',
            'count_index_3', 
            'count_index_4',
            'count_index_5',
            'count_index_6',
            'count_index_7',
            'count_index_8',
            'count_index_9',
            'count_index_10'
        ])
        for triangle in triangles:
            x_triangle, y_triangle = triangle
            row = [x_triangle, y_triangle] + triangles[(x_triangle, y_triangle)]
            writer.writerow(row)

    print('step_2', time.time() - tic)

def step_3(filename, z):

    # turns the triangles from step 2 into hexagons with a single value

    tic = time.time()

    hexagons = {}

    # with gzip.open(f'output_step_2/{filename}-triangle-counts-z{z}.csv.gz', 'rt') as f:
    with open(f'output_step_2/{filename}-triangle-counts-z{z}.csv', 'r') as f:
        reader = csv.reader(f)
        next(reader) # skip header line
        for row in reader:
            elements = [int(element) for element in row]
            x_triangle, y_triangle = elements[0], elements[1]
            counts = elements[2:]
            
            x_hexagon, y_hexagon = get_hexagon(x_triangle, y_triangle)

            if (x_hexagon, y_hexagon) not in hexagons:            
                hexagons[(x_hexagon, y_hexagon)] = counts
            else:
                for i in range(len(counts)):
                    hexagons[(x_hexagon, y_hexagon)][i] += counts[i]
    
    hexagons_majority = {}
    for key in hexagons:
        max_index = hexagons[key].index(max(hexagons[key]))
        hexagons_majority[key] = classes_index_to_value[max_index]
    
    data = {
        'value': [],
        'geometry': []
    }

    for key in hexagons_majority:
        triangles = get_triangles(key[0], key[1])
        for triangle in triangles:
            data['value'].append(hexagons_majority[key])
            data['geometry'].append(get_polygon(triangle[0], triangle[1], z))
    
    gdf = gpd.GeoDataFrame(data, crs='epsg:4326')
    gdf = gdf.dissolve('value')
    gdf.to_file(filename=f'output_step_3/{filename}-polygon-z{z}.gpkg', driver="GPKG")

    print('step_3', time.time() - tic)

def all_in_one(filename, z_min, z_max):
    # from and including z_min to and including z_max

    tic = time.time()

    # Open the geotiff file
    dataset = gdal.Open(f'downloads/{filename}')
    
    # Get the pixel values and the geotransform information
    band = dataset.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    geotransform = dataset.GetGeoTransform()
    # pixel_values = band.ReadAsArray().astype(float)
    pixel_values = band.ReadAsArray().astype(np.uint8)
    
    longitudes = [geotransform[0] + (j + 0.5) * geotransform[1] + (0 + 0.5) * geotransform[2] for j in range(pixel_values.shape[1])]
    latitudes = [geotransform[3] + (0 + 0.5) * geotransform[4] + (i + 0.5) * geotransform[5] for i in range(pixel_values.shape[0])]

    # print('longitudes', longitudes[0], longitudes[-1])
    # print('latitudes', latitudes[0], latitudes[-1])


    data = {
        # z : {
        #     'xs': list of x values
        #     'ys': list of y values
        #     'has_first_or_last_line': boolean if upper or lower line has first or last line
        #     'inner_queue_part': number 
        #     'triangles_current_line': {
        #         'y': number
        #         'counts': {
        #             (x_triangle, y_triangle): list of counts
        #             ...
        #         }
        #     }
        #     'triangles_upper_line': ...
        #     'triangles_lower_line': ...
        #     'hexagons_outer_queue': {
        #        (x_hexagon, y_hexagon): list of counts
        #        ...
        #     }
        #     'hexagons_inner_queue': ...
        # }
    }
    for z in range(z_min, z_max+1):
        data[z] = {
            'has_first_or_last_line': True,
            'inner_queue_part': 0,
            'triangles_current_line': {
                'y': None,
                'counts': {}
            },
            'triangles_lower_line': {
                'y': None,
                'counts': {}
            },
            'triangles_upper_line': {
                'y': None,
                'counts': {}
            },
            'hexagons_outer_queue': {},
            'hexagons_inner_queue': {}
        }

    for z in data:
        data[z]['xs'] = []
        for longitude in longitudes:
            x, _ = get_x_y(lat=0, lon=longitude, z=z)
            data[z]['xs'].append(x)
        
        data[z]['ys'] = []
        for latitude in latitudes:
            _, y = get_x_y(lat=latitude, lon=0, z=z)
            data[z]['ys'].append(y)
    

    for i in range(pixel_values.shape[0]):
        for j in range(pixel_values.shape[1]):
            # Skip nodata pixels
            if pixel_values[i][j] == nodata:
                continue

            value = pixel_values[i][j]

            for z in data:
                x, y = data[z]['xs'][j], data[z]['ys'][i]
                x_triangle, y_triangle = get_triangle(x, y, z)

                if y_triangle != data[z]['triangles_current_line']['y']:
                    if data[z]['triangles_upper_line']['y'] != None:
                        data[z]['has_first_or_last_line'] = False
                    data[z]['triangles_upper_line'] = data[z]['triangles_lower_line']
                    data[z]['triangles_lower_line'] = data[z]['triangles_current_line']
                    data[z]['triangles_current_line'] = {
                        'y': y_triangle,
                        'counts': {}
                    }
                    generate_hexagons(data, z)

                    if len(data[z]['hexagons_inner_queue']) > 1e6:
                        clear_inner_queue(data, z)

                counts = data[z]['triangles_current_line']['counts']
                if (x_triangle, y_triangle) not in counts:            
                    counts[(x_triangle, y_triangle)] = [0 for _ in range(11)]
        
                index = classes[value]['index']
                counts[(x_triangle, y_triangle)][index] += 1


                # # Print the progress
                # pixels_processed = (i * pixel_values.shape[1]) + j + 1
                # total_pixels = pixel_values.shape[0] * pixel_values.shape[1]
                # if pixels_processed % 1000000 == 0:
                #     progress = pixels_processed / total_pixels
                #     print(f"step_1 {progress:.2%} complete")

    for z in data:
        data[z]['has_first_or_last_line'] = True

        data[z]['triangles_upper_line'] = data[z]['triangles_lower_line']
        data[z]['triangles_lower_line'] = data[z]['triangles_current_line']
        data[z]['triangles_current_line'] = {
            'y': data[z]['triangles_lower_line']['y'] + 1,
            'counts': {}
        }
        generate_hexagons(data, z)

        data[z]['triangles_upper_line'] = data[z]['triangles_lower_line']
        data[z]['triangles_lower_line'] = data[z]['triangles_current_line']
        data[z]['triangles_current_line'] = {
            'y': data[z]['triangles_lower_line']['y'] + 1,
            'counts': {}
        }
        generate_hexagons(data, z)

    for z in data:
        clear_inner_queue(data=data, z=z)
        clear_outer_queue(data=data, z=z)

    print('total duration', time.time() - tic)

def generate_hexagons(data, z):
    counts = {
        # (x_hexagon, y_hexagon): list of counts
        # ...
    }

    # sum lower triangles to hexagons
    for x_y_triangle in data[z]['triangles_lower_line']['counts']:
        x_triangle, y_triangle = x_y_triangle
        x_hexagon, y_hexagon = get_hexagon(x_triangle, y_triangle)
        if y_hexagon == y_triangle:
            # the triangle is part of the next line of hexagons
            continue

        if (x_hexagon, y_hexagon) not in counts:
            counts[(x_hexagon, y_hexagon)] = [0 for _ in range(11)]

        for i in range(11):
            counts[(x_hexagon, y_hexagon)][i] += data[z]['triangles_lower_line']['counts'][x_y_triangle][i]

    # sum upper triangles to hexagons
    for x_y_triangle in data[z]['triangles_upper_line']['counts']:
        x_triangle, y_triangle = x_y_triangle
        x_hexagon, y_hexagon = get_hexagon(x_triangle, y_triangle)
        if y_hexagon != y_triangle:
            # the triangle is part of the last line of hexagons
            continue

        if (x_hexagon, y_hexagon) not in counts:            
            counts[(x_hexagon, y_hexagon)] = [0 for _ in range(11)]

        for i in range(11):
            counts[(x_hexagon, y_hexagon)][i] += data[z]['triangles_upper_line']['counts'][x_y_triangle][i]

    # stop here if there were no hexagons generated
    if len(counts) == 0:
        return
    
    # get first and last hexagon, put to outer queue
    x_values = []
    for x_y_hexagon in counts:
        x_hexagon, _ = x_y_hexagon
        x_values.append(x_hexagon)
    x_min = min(x_values)
    x_max = max(x_values)
    y = data[z]['triangles_lower_line']['y'] - 1
    data[z]['hexagons_outer_queue'][(x_min, y)] = counts[(x_min, y)]
    counts.pop((x_min, y))
    data[z]['hexagons_outer_queue'][(x_max, y)] = counts[(x_max, y)]
    counts.pop((x_max, y))
    
    if data[z]['has_first_or_last_line']:
        data[z]['hexagons_outer_queue'].update(counts)
    else:
        data[z]['hexagons_inner_queue'].update(counts)

def clear_inner_queue(data, z):
    hexagons_majority = {}
    for key in data[z]['hexagons_inner_queue']:
        max_index = data[z]['hexagons_inner_queue'][key].index(max(data[z]['hexagons_inner_queue'][key]))
        hexagons_majority[key] = classes_index_to_value[max_index]
    
    data[z]['hexagons_inner_queue'].clear()

    write_data = {
        'value': [],
        'geometry': []
    }

    for key in hexagons_majority:
        triangles = get_triangles(key[0], key[1])
        for triangle in triangles:
            write_data['value'].append(hexagons_majority[key])
            write_data['geometry'].append(get_polygon(triangle[0], triangle[1], z))
    
    gdf = gpd.GeoDataFrame(write_data, crs='epsg:4326')
    gdf = gdf.dissolve('value')

    part = data[z]['inner_queue_part']
    data[z]['inner_queue_part'] += 1
    gdf.to_file(filename=f'output_step_3/{filename}-polygon-inner-part-{part}-z{z}.gpkg', driver="GPKG")



def clear_outer_queue(data, z):
    # writes polygons of the outer hexagons for debugging
    # later should write hexagon counts instead for combination with neighboring outer hexagons

    hexagons_majority = {}
    for key in data[z]['hexagons_outer_queue']:
        max_index = data[z]['hexagons_outer_queue'][key].index(max(data[z]['hexagons_outer_queue'][key]))
        hexagons_majority[key] = classes_index_to_value[max_index]
    
    data[z]['hexagons_outer_queue'].clear()

    write_data = {
        'value': [],
        'geometry': []
    }

    for key in hexagons_majority:
        triangles = get_triangles(key[0], key[1])
        for triangle in triangles:
            write_data['value'].append(hexagons_majority[key])
            write_data['geometry'].append(get_polygon(triangle[0], triangle[1], z))
    
    gdf = gpd.GeoDataFrame(write_data, crs='epsg:4326')
    gdf = gdf.dissolve('value')

    gdf.to_file(filename=f'output_step_3/{filename}-polygon-outer-z{z}.gpkg', driver="GPKG")




# process_image(filename='clipped-esa.tif', z=17)

#step_1(filename='ESA_WorldCover_10m_2020_v100_N66W039_Map.tif', z=17)

# filename = 'clipped-esa.tif'
# filename = 'ESA_WorldCover_10m_2020_v100_N45E006_Map.tif'
# step_1(filename=filename, z=z)
# step_2(filename=filename, z=z)
# step_3(filename=filename, z=z)

# process_image_np(filename='ESA_WorldCover_10m_2020_v100_N66W039_Map.tif', z=17)
# process_image_np(filename='clipped-esa.tif', z=17)

all_in_one(filename=filename, z_min=z, z_max=z)
# print(get_hexagon(x_triangle=5, y_triangle=1))
