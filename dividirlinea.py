import folium
from folium import plugins
import geopandas as gpd
import shapely.ops as shp
import shapely
from shapely.ops import unary_union, substring
from shapely.geometry import Point, Polygon, MultiLineString, MultiPoint, LineString
import pandas as pd
import contextily as cx
import numpy as np
import networkx as nx
import itertools
from fiona.drvsupport import supported_drivers


# df = SHP.set_geometry(
#     SHP.geometry.map(
#         lambda linestring: shapely.ops.transform(lambda x, y, *_: (x, y), linestring)
#     ))

# Para luego, Dividir Linea por puntos --> https://shapely.readthedocs.io/en/stable/manual.html#splitting


def tramos_basico(shpinterno, longtitudtramo):
    l = len(shpinterno)
    nombres = shpinterno["Linea"].unique()
    A = []
    for i in range(l):
        linea = shpinterno.at[i, "geometry"]
        distances = np.arange(0, linea.length, longtitudtramo)
        L = []

        for w in distances:
            b = substring(linea, w, w + longtitudtramo)
            L.append(b)
        aGDF = gpd.GeoDataFrame(
            list(range(len(L))), geometry=L, crs=22185)
        aGDF.columns = ['index', 'geometry']
        aGDF["index"] = aGDF["index"].astype(int)
        aGDF['centroide'] = aGDF['geometry'].centroid
        aGDF['mitad'] = aGDF['geometry'].interpolate((aGDF['geometry'].length) / 2)
        aGDF['mitad'] = aGDF['mitad'].to_crs(4326)
        for c in range(len(aGDF)):
            aGDF.at[c, "inicio"] = Point(np.array(aGDF["geometry"][c].coords)[0])
        aGDF.at[len(aGDF), "inicio"] = Point(np.array(aGDF.at[len(aGDF)-1, "geometry"].coords)[1])
        aGDF.at[len(aGDF)-1, "index"] = aGDF.at[len(aGDF)-2, "index"] + 1
        aGDF["Linea"] = nombres[i]
        aGDF['IDtramo'] = (aGDF["Linea"].astype(str) + "-" + aGDF["index"].astype(str)).str.split('.').str[0]
        auxiliar = aGDF[["inicio", "IDtramo"]]
        auxiliar2 = gpd.GeoDataFrame(auxiliar, geometry="inicio", crs=22185)
        aGDF = aGDF.drop("inicio", axis=1)
        aGDF = pd.merge(aGDF, auxiliar2, how="inner", on="IDtramo" )
        aGDF['inicio'] = aGDF['inicio'].to_crs(4326)
        A.append(aGDF)
    B = pd.concat(A)
    B = B.rename(columns={'geometry': 'tramos'})
    B.drop(["centroide"], axis=1, inplace=True)
    B = B.reset_index()
    B.drop(["level_0", "index"], axis=1, inplace=True)
    return B



