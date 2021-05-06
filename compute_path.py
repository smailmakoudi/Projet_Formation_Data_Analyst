import networkx as nx
import numpy as np
import pandas as pd
import json		
def geocalc(lat0, lon0, lat1, lon1):
    EARTH_R = 6378137
    """Return the distance (in m) between two points
    in geographical coordinates."""
    lat0 = np.radians(lat0)
    lon0 = np.radians(lon0)
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    dlon = lon0 - lon1
    y = np.sqrt((np.cos(lat1) * np.sin(dlon)) ** 2 +
        (np.cos(lat0) * np.sin(lat1) - np.sin(lat0) *
         np.cos(lat1) * np.cos(dlon)) ** 2)
    x = np.sin(lat0) * np.sin(lat1) + \
        np.cos(lat0) * np.cos(lat1) * np.cos(dlon)
    c = np.arctan2(y, x)
    return EARTH_R * c
def wgs84_to_web_mercator(lat, lon):
    '''
    from https://stackoverflow.com/questions/57178783/how-to-plot-latitude-and-longitude-in-bokeh
    Cette fonction calcule les coordonnées
    UTM (mercator) à partir de la longitude et de la latitude
    df : le dataframe où se trouvent les données à convertir
    lon : nom de la colonne pour les longitudes
    lat : nom de la colonne pour les latitudes
    return : retourne le dataframe en y ajoutant les colonnes 
    pour les coordonnées x_utm et y_utm
    '''
    k = 6378137 # rayon de la Terre en mètres
    x = lon * (k * np.pi/180.0)
    y = np.log(np.tan((90 + lat) * np.pi/360.0)) * k
    return x, y


def MetersToLatLon(mx, my):
    "Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum"
    "Source : https://www.maptiler.com/google-maps-coordinates-tile-bounds-projection/"
    originShift = 2 * np.pi * 6378137 / 2.0
    lon = (mx / originShift) * 180.0
    lat = (my / originShift) * 180.0
    lat = 180 / np.pi * (2 * np.arctan( np.exp( lat * np.pi / 180.0)) - np.pi / 2.0)
    return lat, lon
def distance_euclidean(lat0, lon0, lat1, lon1):
    x0, y0 = wgs84_to_web_mercator(lat0, lon0)
    x1, y1 = wgs84_to_web_mercator(lat1, lon1)
    return np.sqrt( (x0-x1)**2 + (y0-y1)**2)
	
def get_path_length(path):
    return np.sum(geocalc(path[1:, 1], path[1:, 0],
                          path[:-1, 1], path[:-1, 0]))
def get_path_length_euclidean(path):
    return np.sum(distance_euclidean(path[1:, 1], path[1:, 0],
                          path[:-1, 1], path[:-1, 0]))

def compute_nodes(sgs, dist_type='euclidean'):
    i = np.argmax([len(sg) for sg in sgs])
    sg = sgs[i]
    # Compute the length of the road segments.
    for n0, n1 in sg.edges:
        path = np.array(json.loads(sg[n0][n1]['Json'])['coordinates'])
        if dist_type == 'euclidean':
            distance = get_path_length_euclidean(path)
        elif dist_type == 'haversine':
            distance = get_path_length(path)
        else:
            print('dist_type doit être euclidean ou haversine')
        sg.edges[n0, n1]['distance'] = distance
    nodes = np.array(sg.nodes())
    return nodes, sg
	
def get_linepath(sg, path):
    """Return the positions along a path."""
    p_list = []
    p_inter = []
    curp = None
    for i in range(len(path) - 1):
        p = np.array(json.loads(sg[path[i]][path[i + 1]]['Json'])['coordinates'])
        if curp is None:
            curp = p
        if (np.sum((p[0] - curp) ** 2) >
                np.sum((p[-1] - curp) ** 2)):
            p = p[::-1, :]
        p_list.append(p)
        for x, y in sg.edges(path[i]):
            if y != path[i + 1]:
               v = np.array(json.loads(sg[x][y]['Json'])['coordinates'])
               p_inter.append(v)
        curp = p[-1]
    return np.vstack(p_list), np.vstack(p_inter)
	
def compute_path(nodes, sg, o_long, o_lat, d_long, d_lat):
    pos0 = (o_lat, o_long)
    pos1 = (d_lat, d_long)
    # Get the closest nodes in the graph.
    pos0_i = np.argmin(
        np.sum((nodes[:, ::-1] - pos0)**2, axis=1))
    pos1_i = np.argmin(
        np.sum((nodes[:, ::-1] - pos1)**2, axis=1))
    # Compute the shortest path.
    path = nx.shortest_path(
        sg,
        source=tuple(nodes[pos0_i]),
        target=tuple(nodes[pos1_i]),
        weight='distance')
    length = nx.shortest_path_length(
        sg,
        source=tuple(nodes[pos0_i]),
        target=tuple(nodes[pos1_i]),
        weight='distance')
    return path, length