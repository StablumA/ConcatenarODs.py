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


warnings.filterwarnings('ignore')
supported_drivers['LIBKML'] = 'rw'



# Funciones
def generate_color_scale(min_value, max_value):
    colormap = linear.RdYlGn_11.scale(min_value, max_value)
    return colormap

def calcular_distancia(row):
    if pd.notnull(row['next_latitude']) and pd.notnull(row['next_longitude']):
        return (geodesic((row['latitude'], row['longitude']),
                        (row['next_latitude'], row['next_longitude'])).km)*1.27
    else:
        return None

# Parametros
direcciongps ="Venado/Rawdata/Registros GPS 1º SEM 2023 .xlsx"
direccionlineasnuevas = "Venado/SHP Files/Lineas/LineasNuevas/LineasNuevas1F.shp"
direccionlineasgps = "Venado/SHP Files/Lineas/RecorridosF.shp"
direcciondfseparaciones = "Venado/Rawdata/SeparacionParadas.xlsx"  # direccion del df de separaciones
direccionseparacionkml  = "Venado/SHP Files/Velocidades/ZonasParadasVT.kml"
direcciondestino = "Venado/Output/Pruebas Nuevos Sistemas/Alternativa1Septiembre/NodosVelocidades.xlsx"
longitudtramo = 100
tipodedia = "Habil"
diasferiados = []
listahoras = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]

# Cargar datos
gps = pd.read_excel(direcciongps)
gps = gps[["Interno", "Ramal", "Fecha completa", "Longitud", "Latitud"]]
gps = gps.dropna()
gps = gps.drop_duplicates(subset=["Interno", "Ramal", "Fecha completa", "Longitud", "Latitud"])
gps.rename(columns={"Longitud": "longitude", "Latitud": "latitude", "Fecha completa": "date_time", "Interno" : "interno", "Ramal" : "idlinea"}, inplace=True)
gps = gps[gps["idlinea"] != 0]
gps["date_time"] = pd.to_datetime(gps["date_time"])
gps["DIA"] = gps["date_time"].dt.day
gps["HORA"] = gps["date_time"].dt.hour
gps['nDIA']=gps['date_time'].dt.day_name() #le agrego los dias por nombre


if tipodedia == "Habil":
    gps = gps[gps['nDIA'].isin(['Wednesday', 'Thursday', 'Friday', 'Monday',
           'Tuesday'])]
    gps = gps[~gps["DIA"].isin(diasferiados)]
elif tipodedia == "Sabado":
    gps = gps[gps['nDIA'].isin(['Saturday'])]
    gps = gps[~gps["DIA"].isin(diasferiados)]
else:
    aux1 = gps[gps['nDIA'].isin(['Sunday'])]
    aux2 = gps[gps['nDIA'].isin(['Wednesday', 'Thursday', 'Friday', 'Monday',
                                     'Tuesday', "Saturday"])]
    aux2 = aux2[aux2["DIA"].isin(diasferiados)]
    gps = pd.concat([aux1, aux2])
SHP = gpd.read_file(direccionlineasnuevas)  # cargo el shp de las lineas
SHP = SHP.rename(columns={'Name': 'Lineas'})
SHP = SHP.to_crs(epsg=22185)
Lineas = gps.idlinea.unique()
df = SHP.set_geometry(
    SHP.geometry.map(
        lambda linestring: shapely.ops.transform(lambda x, y, *_: (x, y), linestring)
    ))
shpnodos = tramos(df, longitudtramo, direcciondfseparaciones, direccionseparacionkml)
shpnodos = shpnodos.set_geometry("inicio")
shpnodos = shpnodos.to_crs(epsg=4326)


from geopy.distance import geodesic



# Calcular la distancia entre filas consecutivas
A = pd.DataFrame(columns=["idlinea", "longitude", "latitude", "HORA", "Velocidad"])
for l in Lineas:
    gpsinterno = gps[gps["idlinea"] == l]
    for d in gps.DIA.unique():
        gpsinternodia = gpsinterno[gpsinterno["DIA"] == d]
        for i in gpsinternodia.interno.unique():
            gpsinternodiainterno = gpsinternodia[gpsinternodia["interno"] == i]
            gpsinternodiainterno = gpsinternodiainterno.sort_values(by=["date_time"])
            gpsinternodiainterno['next_latitude'] = gpsinternodiainterno['latitude'].shift(-1)
            gpsinternodiainterno['next_longitude'] = gpsinternodiainterno['longitude'].shift(-1)
            gpsinternodiainterno['distancia'] = gpsinternodiainterno.apply(calcular_distancia, axis=1)
            gpsinternodiainterno['tiempo_minutos'] = (gpsinternodiainterno['date_time'].shift(-1) - gpsinternodiainterno['date_time']).dt.total_seconds() / 60
            gpsinternodiainterno = gpsinternodiainterno[(gpsinternodiainterno["tiempo_minutos"] > 0) & (gpsinternodiainterno["tiempo_minutos"] < 10)]
            gpsinternodiainterno["Velocidad"] = gpsinternodiainterno["distancia"] / gpsinternodiainterno["tiempo_minutos"] * 60
            gpsinternodiainterno = gpsinternodiainterno[gpsinternodiainterno["Velocidad"] > 8]
            gpsinternodiainterno = gpsinternodiainterno[gpsinternodiainterno["Velocidad"] < 50]
            gpsinternodiainterno = gpsinternodiainterno[["idlinea", "longitude", "latitude", "HORA", "Velocidad"]]
            A = pd.concat([A, gpsinternodiainterno], axis=0)


A = A.reset_index(drop=True)
A = gpd.GeoDataFrame(A, geometry=gpd.points_from_xy(A.longitude, A.latitude), crs="EPSG:4326")
A = A.set_crs(epsg=4326)

velocidadesporpuntos = gpd.sjoin_nearest(A, shpnodos, how="left")
velocidadpromedio = velocidadesporpuntos.groupby(["Lineas", "IDtramo", "HORA"])["Velocidad"].mean().reset_index()
desvio = velocidadesporpuntos.groupby(["Lineas", "IDtramo", "HORA"])["Velocidad"].std().reset_index(name='desvio')
conteo = velocidadesporpuntos.groupby(["Lineas", "IDtramo", "HORA"]).size().reset_index(name='conteo')
velocidadpromedio = velocidadpromedio.merge(conteo, on=["Lineas", "IDtramo", "HORA"])
velocidadpromedio = velocidadpromedio.merge(desvio, on=["Lineas", "IDtramo", "HORA"])
velocidadpromedio[["NADA", "ORDEN"]] = velocidadpromedio["IDtramo"].str.split('-', -1, expand=True)
velocidadpromedio = velocidadpromedio.drop(columns=["NADA"])
velocidadpromedio["ORDEN"] = velocidadpromedio["ORDEN"].astype(int)
velocidadpromedio = velocidadpromedio[velocidadpromedio["HORA"].isin(listahoras)]
velocidadpromedio = velocidadpromedio.sort_values(by=["Lineas", "ORDEN", "HORA"])


shpnodosh = shpnodos[['Lineas', 'IDtramo', 'inicio', "ZonaParadas", "factorx", "Distancia"]]
A = []
for h in listahoras:
    dfinterno = shpnodosh.copy(deep=True)
    dfinterno["HORA"] = h
    A.append(dfinterno)
shpnodosh = pd.concat(A, axis=0)

velocidadpromediogeo = pd.merge(velocidadpromedio, shpnodosh, on=["Lineas", "IDtramo", "HORA"], how="outer")
promedios = velocidadpromediogeo.groupby('Lineas')['Velocidad'].mean()
velocidadpromediogeo['Velocidad'] = velocidadpromediogeo.apply(lambda row: promedios[row['Lineas']] if np.isnan(row['Velocidad']) else row['Velocidad'], axis=1)
velocidadpromediogeo = gpd.GeoDataFrame(velocidadpromediogeo, geometry=velocidadpromediogeo["inicio"], crs="EPSG:4326")
velocidadpromediogeo["latitude"] = velocidadpromediogeo["inicio"].y
velocidadpromediogeo["longitude"] = velocidadpromediogeo["inicio"].x
velocidadpromediogeo[["NADA", "ORDEN"]] = velocidadpromediogeo["IDtramo"].str.split('-', -1, expand=True)
velocidadpromediogeo = velocidadpromediogeo.drop(columns=["NADA"])
velocidadpromediogeo["ORDEN"] = velocidadpromediogeo["ORDEN"].astype(int)
velocidadpromediogeo = velocidadpromediogeo.fillna(0)
velocidadpromediogeo = velocidadpromediogeo.sort_values(by=["Lineas", "ORDEN", "HORA"])
velocidadpromediogeo = velocidadpromediogeo[['Lineas', 'IDtramo', 'ORDEN', 'HORA', 'Velocidad', 'desvio', 'conteo', 'latitude', 'longitude', 'geometry', "ZonaParadas", "factorx", "Distancia"]]
velocidadpromediogeo.rename(columns={"Lineas" : "idlinea"}, inplace=True)
velocidadpromediogeo.to_excel(direcciondestino, index=False)


mapa3 = folium.Map([-33.74595770280119, -61.96865467712429], tiles="OpenStreetMap", zoom_start=14 )

listadelistasdefg = []
Lineas = velocidadpromediogeo["idlinea"].unique()
for l in Lineas:
        listafg = []
        velocidadpromedioint = velocidadpromediogeo[velocidadpromediogeo["idlinea"] == l].reset_index(drop=True)
        velocidadpromedioint = velocidadpromedioint.sort_values(by=["ORDEN"]).reset_index()
        velocidadpromedioint['VelocidadReversa'] = velocidadpromedioint['Velocidad'].rank(ascending=True)
        min = velocidadpromedioint["VelocidadReversa"].min()
        max = velocidadpromedioint["VelocidadReversa"].max()
        colormap = generate_color_scale(min, max)
        horas = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
        for h in horas:
            FG = folium.FeatureGroup(name="Hora: " + str(h), control=True, show=False)
            velocidadpromediointg = velocidadpromedioint[velocidadpromedioint["HORA"] == h].reset_index(drop=True)
            velocidadpromediointg = velocidadpromediointg.sort_values(by=["ORDEN"]).reset_index()
            # velocidadpromediointg['VelocidadReversa'] = velocidadpromediointg['Velocidad'].rank(ascending=False)
            # min = velocidadpromediointg["VelocidadReversa"].min()
            # max = velocidadpromediointg["VelocidadReversa"].max()
            # colormap = generate_color_scale(min, max)
            for i, (lo, la) in enumerate(zip(velocidadpromediointg["longitude"], velocidadpromediointg["latitude"])):
                color = colormap(velocidadpromediointg["VelocidadReversa"][i])
                try:
                    folium.PolyLine(locations=[[la, lo], [velocidadpromediointg["latitude"][i + 1], velocidadpromediointg["longitude"][i + 1]]], color=color, weight=3.5, opacity=0.8, popup=velocidadpromediointg.at[i, "Velocidad"]).add_to(FG)
                except KeyError:
                    pass
            FG.add_to(mapa3)
            listafg.append(FG)
        listadelistasdefg.append(listafg)


GroupedLayerControl(
    groups={l: listadelistasdefg[i] for i, l in enumerate(Lineas)},
    exclusive_groups=False,
    collapsed=True
).add_to(mapa3)
loc = "Mapa de calor de las velocidades por hora"
title_html = '''
                <h3 align="center" style="font-size:16px"><b>{}</b></h3>
               '''.format(loc)
mapa3.get_root().html.add_child(folium.Element(title_html))

mapa3.save("Venado/Output/HeatVelocidadesAlternativa1.html") # Guardo.
print("Mapa guardado")
