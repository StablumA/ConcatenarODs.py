import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as cx
import folium
from folium import plugins
import pandas as pd
from shapely.geometry import Point
import pygeos
import numpy
from fiona.drvsupport import supported_drivers
supported_drivers['LIBKML'] = 'rw'

#SHP = gpd.read_file("SHP Files/Recorridos/Recorridos2021/Recorridos octubre 2021.kml", driver = "LIBKML")
SHP = gpd.read_file("SHP Files/TUP 2023/TUP - 2023.shp", driver = "LIBKML")
#SHP = SHP.drop([32,33,34,35])
#lista = [16,13,8,1,5,18,14,4,9,15,132,3,10,102,11,112,2,21]

#arr= SHP["Name"].unique()
#SHP = SHP.replace(arr, lista)
SHP= SHP.rename(columns={'NOMBRE': 'Linea'})
SHP = SHP[["Linea", "geometry"]]
SHP["Linea"].replace(["LINEA 1A/1B","LINEA 2A","LINEA 2B","LINEA 3A/3B","LINEA 4A","LINEA 4B"],["A","BA","BB","C","DA","DB"],inplace=True)
lineas = SHP["Linea"].unique()
for i in lineas:
    a = SHP[SHP['Linea'] == i]
    a.to_file("Preprocessdata/Recorridos/Recorrido_"+ str(i) + ".shp")

SHP.to_file("SHP Files/RecorridosF.shp")

#visualizacion de recorridos
mapa = folium.Map([-33.63238910, -61.69945910], tiles="OpenStreetMap", zoom_start=10 )
style1 = {'fillColor': '#5b5b5f', 'color': '#5b5b5f'}

for i in lineas:
    a = gpd.read_file(
        "Preprocessdata/Recorridos/Recorrido_" + str(i) + ".shp")  # Leer el archivo con el recorrido
    fg = folium.FeatureGroup("LINEA" + " " + str(i))
    folium.GeoJson(data=a["geometry"], style_function=lambda x: style1).add_to(fg)
    fg.add_to(mapa)
folium.LayerControl().add_to(mapa)
loc = "Recorridos Abril 2023"
title_html = '''
                    <h3 align="center" style="font-size:16px"><b>{}</b></h3>
                   '''.format(loc)
mapa.get_root().html.add_child(folium.Element(title_html))
mapa.save("Output/Recorridos" + ".html")  # Guardo.
print("Mapa Generado")
