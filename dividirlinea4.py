import geopandas as gpd
import shapely
from shapely.ops import substring
from shapely.geometry import Point, Polygon, MultiLineString, MultiPoint, LineString
import pandas as pd
import numpy as np
from fiona.drvsupport import supported_drivers
supported_drivers['LIBKML'] = 'rw'




# # # #borrar

# SHP = gpd.read_file("SHP Files/Nuevas Lineas/Nuevas LineasSR.kml") #cargo el shp de las lineas
# SHP=SHP.rename(columns={'Name':'Lineas'})
# SHP = SHP.to_crs(epsg=22185)
# Lineas = SHP["Lineas"].unique()
# # SHP["Lineas"] = SHP["Linea"]
# # SHP.drop([8,9], inplace=True)
#
# df = SHP.set_geometry(
#     SHP.geometry.map(
#         lambda linestring: shapely.ops.transform(lambda x, y, *_: (x, y), linestring)
#     ))




def tramos(shpinterno, longitudtramo, direcciondfconseparaciones, dseparadas):
    separaciones = gpd.read_file(dseparadas, driver="LIBKML")
    separaciones = separaciones[["geometry", "Name"]]
    separaciones = separaciones.rename(columns={"Name": "ZonaParadas"})
    dfconseparaciones = pd.read_excel(direcciondfconseparaciones)
    l = len(shpinterno)
    nombres = shpinterno["Lineas"].unique()
    A = []
    for i in range(l):
        linea = shpinterno.at[i, "geometry"]
        distances = np.arange(0, linea.length, longitudtramo)
        L = []
        for w in distances:
            b = substring(linea, w, w + longitudtramo)
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
        aGDF["Lineas"] = nombres[i]
        aGDF['IDtramo'] = (aGDF["Lineas"].astype(str) + "-" + aGDF["index"].astype(str)).str.split('.').str[0]
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
    B.set_geometry("inicio", inplace=True)
    B = gpd.sjoin(B, separaciones, how="left", op="intersects")
    B["ZonaParadas"] = B["ZonaParadas"].fillna("Normal")
    B = B[["Lineas", "IDtramo", "inicio", "ZonaParadas"]]
    zonas = B["ZonaParadas"].unique()
    C = pd.DataFrame()
    for n in nombres:
        auxinterno = B[B["Lineas"] == n]
        D = pd.DataFrame()
        for z in zonas:
            auxinternointerno = auxinterno[auxinterno["ZonaParadas"] == z]
            if len(auxinternointerno) > 0:
                factorx = dfconseparaciones[dfconseparaciones["Zona"] == z]["Factor"].values[0]
                auxinternoav = auxinternointerno[auxinternointerno["ZonaParadas"] == z].reset_index()
                df_pares = auxinternoav.iloc[::factorx].reset_index()
                df_pares["Distancia"] = longitudtramo * factorx
                df_pares["factorx"] = factorx
                D = pd.concat([D, df_pares]).drop(["index", "level_0"], axis=1)
            else:
                continue
        D["indice"] = D["IDtramo"].str.split('-').str[1].astype(int)
        D = D.sort_values(by=["indice"]).drop(["indice", "IDtramo"], axis=1).reset_index().reset_index()
        D["IDtramo"] = (D["Lineas"].astype(str) + "-" + D["level_0"].astype(str)).str.split('.').str[0]
        D = D.drop(["level_0", "index"], axis=1)
        C = pd.concat([C, D])
    return C