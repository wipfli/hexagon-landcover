from math import floor, radians, log, cos, pi, tan, sqrt, atan, exp, degrees
import geopandas as gpd
from shapely.geometry import Polygon
from osgeo import gdal
from collections import Counter
import numpy as np
import csv
import gzip

import time


def get_x_max(z):
    return 2 ** z

def get_triangle(x, y, z):
    x_max = 2 ** z
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
            if (x - i) > (y -j):
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

            # # Print the progress
            # pixels_processed = (i * pixel_values.shape[1]) + j + 1
            # total_pixels = pixel_values.shape[0] * pixel_values.shape[1]
            # if pixels_processed % 10000 == 0:
            #     progress = pixels_processed / total_pixels
            #     print(f"{filename}: {progress:.2%} complete", end="\r")
    
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
            

def process_image_csv(filename, z):

    # Open the geotiff file
    dataset = gdal.Open(filename)
    
    # Get the pixel values and the geotransform information
    band = dataset.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    geotransform = dataset.GetGeoTransform()
    pixel_values = band.ReadAsArray().astype(float)
    

    tic = time.time()

    with gzip.open('output.csv.gz', 'wt') as f:
        writer = csv.writer(f)
        writer.writerow(['x_triangle', 'y_triangle', 'value'])
            
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
                writer.writerow([x_triangle, y_triangle, int(pixel_values[i][j])])
                # Print the progress
                pixels_processed = (i * pixel_values.shape[1]) + j + 1
                total_pixels = pixel_values.shape[0] * pixel_values.shape[1]
                if pixels_processed % 10000 == 0:
                    progress = pixels_processed / total_pixels
                    print(f"{filename}: {progress:.2%} complete", end="\r")

    print('total time', time.time() - tic)
         

process_image_csv(filename='ESA_WorldCover_10m_2020_v100_N66W039_Map.tif', z=17)
# process_image_csv(filename='clipped-esa.tif', z=17)
