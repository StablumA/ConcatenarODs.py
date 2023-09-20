import pandas as pd
import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
import shapely
import numpy as np
from shapely.ops import transform
from fiona.drvsupport import supported_drivers
from dividirlinea4 import tramos
from dividirlinea import tramos_basico
import warnings
import folium
from folium.plugins import GroupedLayerControl
from branca.colormap import linear
import matplotlib.pyplot as plt

supported_drivers['kml'] = 'rw' # enable KML support which is disabled by default
supported_drivers['KML'] = 'rw'
warnings.filterwarnings('ignore')
supported_drivers['LIBKML'] = 'rw'


# preprocesamos el SHP también
SHP = gpd.read_file(direccionlineassalida)
SHP = SHP.to_crs(epsg=4326)
SHP.rename(columns={"NOMBRE":"Name"},inplace=True)
SHP.drop("id",axis=1,inplace=True)
SHP["Name"] = SHP["Name"].str.split(" ",expand=True)[1]
aux = SHP["Name"].str.split("/",expand=True)[0].reset_index()
aux["geometry"] = SHP["geometry"]
aux.rename(columns={0:"Name"},inplace=True)
aux2=SHP["Name"].str.split("/",expand=True)[1]
aux2= pd.DataFrame(aux2[aux2.values!=None])
aux2["geometry"]=0
for i in aux2.index:
    aux2.at[i,"geometry"] = SHP.at[i,"geometry"]
aux2=aux2.reset_index()
aux2.rename(columns={1:"Name"},inplace=True)
SHP=pd.concat([aux,aux2],axis=0).reset_index()
SHP.drop(["level_0","index"],inplace=True,axis=1)
SHP.sort_values(by=["Name"],inplace=True)
df = SHP.set_geometry(
     SHP.geometry.map(
         lambda linestring: shapely.ops.transform(lambda x, y, *_: (x, y), linestring)
     ))
SHP = gpd.GeoDataFrame(df,geometry="geometry")
SHP.crs = 'EPSG:4326'
SHP["Lineas"] = SHP["Name"]
SHP.to_file(direccionlineassalida)

Lineas = SHP["Name"].unique()
nuevaslineas = ['A', 'A', 'BA', 'BB', 'C', 'C', 'DA', 'DB']
SHP["Name"] = SHP["Name"].replace(Lineas, nuevaslineas)
SHP["Linea"] = SHP["Name"]
SHP["Lineas"] = SHP["Linea"]
SHP= SHP[["Name","Linea","geometry","Lineas"]]
for l in SHP["Lineas"].unique():
    shpinterno = SHP[SHP["Linea"] == l]
    shpinterno.crs= 'EPSG:4326'
    shpinterno.to_file("Preprocessdata/Recorridos/Recorrido_"+ str(l) + ".shp")






# INTERPROCESAMIENTO
OD1 = pd.read_csv("Output/ODfinalmayo.csv")
OD2= pd.read_csv("Output/ODfinalSeptiembre.csv")
OD3 = pd.read_csv("Output/ODfinalNoviembre.csv")
ODaux = pd.concat([OD1,OD2],axis=0)
ODFINAL = pd.concat([ODaux,OD3],axis=0)
ODFINAL.to_csv("Output/OD_MaySepNovFINAL.csv")

ExL1 = pd.read_excel("Output/estadisticosodporlineaMayo.xlsx")
ExL2 = pd.read_excel("Output/estadisticosodporlineaSeptiembre.xlsx")
ExL3 = pd.read_excel("Output/estadisticosodporlineaNoviembre.xlsx")

ExLF = pd.DataFrame()
ExLF["idlinea"] = ExL1["idlinea"]
ExLF["CuentaOD"] = ExL1["CuentaOD"] + ExL2["CuentaOD"] + ExL3["CuentaOD"]
ExLF["CuentaDF"] = ExL1["CuentaDF"] + ExL2["CuentaDF"] + ExL3["CuentaDF"]
ExLF["% por linea"] =(  ExL1["% por linea"] + ExL2["% por linea"] + ExL3["% por linea"]   )/3
ExLF["Factor"] =(  ExL1["Factor"] + ExL2["Factor"] + ExL3["Factor"]   )/3
ExLF.to_excel("Output/estadisticosodporlineaMSN.xlsx")


RC=gpd.read_file("SHP Files/Radios Censales VT/RC_VT_Cortado.shp")
gdfRC= gpd.GeoDataFrame(RC,geometry="geometry")
destinomapaRC = "Output/mapaRC.html"
mapaRC = folium.Map(location=[-31.6324, -60.6991], zoom_start=12)
gdfRC["geometry"] = gdfRC['geometry'].to_crs(epsg=3857)        # aparto a la geometría
folium.GeoJson(gdfRC).add_to(mapaRC)  # Grafico los polígonos

# importo od
mapaOD = folium.Map(location=[-31.6324, -60.6991], zoom_start=12)
odexp = pd.read_csv("Output/ODexpVT.csv")
#odexp = pd.read_csv("Output/OD_MaySepNovFINAL.csv")
fgod = folium.FeatureGroup("puntosOD")
for i in range(0,len(odexp),6):
    lat,long = odexp.at[i,"latitude_Origen"],odexp.at[i,"longitude_Origen"]
    fgod.add_child(folium.CircleMarker(
            location=(lat,long),
            radius=0.2,
            fill=True,
            color="#12bc8e" ))
fgod.add_to(mapaOD)
mapaOD.save("Output/MapaODexp.html")
print('mapa listo en '+ "Output/MapaODexp.html")





direccionlineas = "SHP Files/RecorridosF.shp"  # direccion de las lineas
direccionod = "Output/ODexpVT.csv"  # direccion de la OD
direccionparadas = "Preprocessdata/paradasVF.xlsx"  # direccion de las paradas
direcciondfseparaciones = "Rawdata/SeparacionParadas.xlsx"  # direccion del df de separaciones  # modificado # puse todas a 200 # antes era 6 3 2 respectivamente
direccionseparacionkml = "SHP Files/Velocidades/ZonaParadasVT.kml"
## MAPA
import folium
from dividirlinea4 import tramos
coordenadasmapa = [-33.732235653522516, -61.9855704545104]
longitudtramo=100
directoriosalida = "Preprocessdata/"


shpnodos = tramos(df, longitudtramo, direcciondfseparaciones, direccionseparacionkml)

mapa1 = folium.Map(coordenadasmapa, tiles="OpenStreetMap", zoom_start=14)

for lv in Lineas:
    fg = folium.FeatureGroup("LINEA" + " " + lv, show=False)
    shpnodosinterno = shpnodos[shpnodos["Lineas"] == lv].reset_index()
    for j, nmg in enumerate(shpnodosinterno["inicio"]):
        fg.add_child(folium.CircleMarker(
            location=(nmg.y, nmg.x),
            radius=4,
            fill=True,
            color="#12bc8e", popup=shpnodosinterno.at[j, "IDtramo"]
        ))
    fg.add_to(mapa1)
folium.LayerControl().add_to(mapa1)
loc = "Paradas Magicas"
title_html = '''
                <h3 align="center" style="font-size:16px"><b>{}</b></h3>
               '''.format(loc)
mapa1.get_root().html.add_child(folium.Element(title_html))

mapa1.save(directoriosalida + "ParadasMagicas.html")  # Guardo.

