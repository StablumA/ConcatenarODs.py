import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
import pandas as pd
import ast
import warnings
import time
import numpy as np
from shapely import Polygon
import h3
from fiona.drvsupport import supported_drivers
from openpyxl import load_workbook
from geojson import Feature, Point, FeatureCollection
import json
import folium
import matplotlib.pyplot as plt
import matplotlib
supported_drivers['LIBKML'] = 'rw'
warnings.filterwarnings('ignore')



# Funciones

def get_color(custom_cm, val, vmin, vmax):
    return matplotlib.colors.to_hex(custom_cm((val - vmin) / (vmax - vmin)))


def h3pol(h3_id):
    strings_separados = h3_id.split("-")
    h3code = strings_separados[0]
    return Polygon(h3.h3_to_geo_boundary(h3code, geo_json=True))

def hexagons_dataframe_to_geojson(df_hex, column_name="value"):
    """
    Produce the GeoJSON for a dataframe, constructing the geometry from the "hex_id" column
    and with a property matching the one in column_name
    """
    list_features = []

    for i, row in df_hex.iterrows():
        try:
            geometry_for_row = {"type": "Polygon",
                                "coordinates": [h3.h3_to_geo_boundary(h=row["h3_code"], geo_json=True)]}
            feature = Feature(geometry=geometry_for_row, id=row["h3_code"], properties={column_name: row[column_name]})
            list_features.append(feature)
        except:
            print("An exception occurred for hex " + row["h3_code"])

    feat_collection = FeatureCollection(list_features)
    geojson_result = json.dumps(feat_collection)
    return geojson_result

def choropleth_map(df_aggreg, column_name="value", border_color='black', fill_opacity=0.4, color_map_name="RdYlGn",
                   initial_map=None):
    """
    Creates choropleth maps given the aggregated data. initial_map can be an existing map to draw on top of.
    """
    # colormap
    min_value = df_aggreg[column_name].min()
    max_value = df_aggreg[column_name].max()
    mean_value = df_aggreg[column_name].mean()
    print(f"Colour column min value {min_value}, max value {max_value}, mean value {mean_value}")
    print(f"Hexagon cell count: {df_aggreg['h3_code'].nunique()}")

    # the name of the layer just needs to be unique, put something silly there for now:
    name_layer = "Choropleth " + str(df_aggreg)

    if initial_map is None:
        initial_map = folium.Map(location=[47, 4], zoom_start=5.5, tiles="cartodbpositron")

    # create geojson data from dataframe
    geojson_data = hexagons_dataframe_to_geojson(df_hex=df_aggreg, column_name=column_name)

    # color_map_name 'Blues' for now, many more at https://matplotlib.org/stable/tutorials/colors/colormaps.html to choose from!
    custom_cm = matplotlib.cm.get_cmap(color_map_name)

    folium.GeoJson(
        geojson_data,
        style_function=lambda feature: {
            'fillColor': get_color(custom_cm, feature['properties'][column_name], vmin=max_value, vmax=min_value),
            'color': border_color,
            'weight': 1,
            'fillOpacity': fill_opacity
        },
        name=name_layer, tooltip=folium.GeoJsonTooltip(fields=[column_name], aliases=[f"{column_name}: "])
    ).add_to(initial_map)
    return initial_map

# Parametros
intervalofrecuencia = 15 # minutos
capacidadcoche = 40  # personas promedio sobre el coche
velocidadcaminando = 5 # km/h
longitudturno = 1 # hora
codigoe = "VT_GTRX_1000h9_ExpNormal_Alternativa2"
codigos = "VT_FA_1000h9_ExpNormal_Alternativa2"
directorio = "Venado/Output/Pruebas Nuevos Sistemas/Alternativa2Septiembre/"
sistemanuevo = False
horaminima = 5
coordenadasmapa = [-33.74595770280119, -61.96865467712429]



# Direccion entradas
direccionlineas = "Venado/SHP Files/Lineas/LineasNuevas/LineasNuevas2F.shp"  # direccion de las lineas
direccionrutas = directorio + codigoe + "RutasCompletas.csv"  # direccion de las rutas
direccionparadas = directorio + codigoe + "ParadasPromedio.csv"
direccionfactorecorreccion = "Venado/Output/estadisticosodporlineaMSN.xlsx"
direccionnodos = directorio + codigoe + "nodos.xlsx"
direcciongeolocalizacion = "Venado/Preprocessdata/paradasVT.xlsx"  # direccion de las paradas
direccionmetricas = directorio + codigoe + "Metricas.xlsx"
direccionnodosvelocidad = directorio + "NodosVelocidades.xlsx"
direccionperfilv = "Venado/Preprocessdata/PerfilV3.xlsx"
direccionpmporhora = directorio + codigoe + "CargaPMPorHora.csv"
direccionshpcarrilrapido = ""
# Direccion salidas
directoriosalida = directorio + codigos
direccioncargat = directorio + codigos + "CargaT.csv"

# Carga de datos
start = time.time()
RutasCompletas = pd.read_csv(direccionrutas)
paradasagrupadastpromedio = pd.read_csv(direccionparadas)
factorescorreccion = pd.read_excel(direccionfactorecorreccion)
nodos = pd.read_excel(direccionnodos)
SHP = gpd.read_file(direccionlineas)  # cargo el shp de las lineas si es shp normal si es kml agregar driver='LIBKML'
SHP = SHP.rename(columns={'Linea': 'Lineas'})
SHP = SHP.to_crs(epsg=22185)
Lineas = SHP["Lineas"].unique()
SHP = SHP[['Lineas', 'geometry']]

for i, row in SHP.iterrows():
    length = row['geometry'].length
    SHP.at[i, 'Longitud'] = length

CargaT = pd.merge(RutasCompletas, paradasagrupadastpromedio, on=["Nodo Origen", "Nodo Destino"], how="outer")
CargaT["CargaProporcional"] = CargaT["Cuenta"] * CargaT["p"]
CargaT = CargaT[CargaT["LineaUtilizada"] != "Ninguna"]
CargaT = CargaT.dropna(subset=["Ruta", "CargaProporcional"])
CargaT = CargaT[~CargaT["LineaOriginal"].isin([np.nan, "CONICET-2", "CONICET-9"])]
# CargaT["LineaOriginal"] = CargaT["LineaOriginal"].astype(int) # esta linea depende del sistema evaluado SFE si VT no
CargaT = pd.merge(CargaT, factorescorreccion, left_on="LineaOriginal", right_on="idlinea", how="left")
CargaT["CargaProporcional"] = CargaT["CargaProporcional"] / CargaT["factor"]
CargaT.dropna(subset=["CargaProporcional"], inplace=True)
CargaT["Ruta"] = CargaT["Ruta"].astype(str)
CargaT['Ruta'] = CargaT['Ruta'].apply(ast.literal_eval)
CargaT = CargaT[["LineaUtilizada", "Ruta", "CargaProporcional", "HORA", "MINUTO", "Nodo Origen", "Nodo Destino", "DistanciaCaminando"]]
CargaT = CargaT[CargaT["HORA"] >= horaminima]
CargaSuma = CargaT.CargaProporcional.sum()
CargaPorLinea = CargaT.groupby(["LineaUtilizada"])["CargaProporcional"].sum().reset_index()
CargaT["HexagonoOrigen"] = CargaT["Nodo Origen"].str.split("-", n = 1, expand = True)[0]
CargaPorHexagono = CargaT.groupby(["HexagonoOrigen"])["CargaProporcional"].sum().reset_index()
distancia_caminando = CargaT[["LineaUtilizada", "Ruta", "CargaProporcional", "HORA", "MINUTO", "Nodo Origen", "Nodo Destino", "DistanciaCaminando", "HexagonoOrigen"]]

distancia_caminando["CaminataProporcional"] = distancia_caminando["CargaProporcional"] * distancia_caminando["DistanciaCaminando"]
CaminataPromedio = distancia_caminando["CaminataProporcional"].sum() / (distancia_caminando["CargaProporcional"].sum() * 2)

distancia_caminando_porlinea = distancia_caminando.groupby(["LineaUtilizada"])["CaminataProporcional"].sum().reset_index()
distancia_caminando_porlinea = pd.merge(distancia_caminando_porlinea, CargaPorLinea, on="LineaUtilizada", how="left")
distancia_caminando_porlinea["CaminataProporcional"] = distancia_caminando_porlinea["CaminataProporcional"] / (distancia_caminando_porlinea["CargaProporcional"] * 2)
distancia_caminando_hexagono = distancia_caminando.groupby(["HexagonoOrigen"])["CaminataProporcional"].sum().reset_index()
distancia_caminando_hexagono = pd.merge(distancia_caminando_hexagono, CargaPorHexagono, on="HexagonoOrigen", how="left")
distancia_caminando_hexagono["CaminataProporcional"] = distancia_caminando_hexagono["CaminataProporcional"] / (distancia_caminando_hexagono["CargaProporcional"] * 2)


distancia_caminando_hexagono["geometry"] = distancia_caminando_hexagono["HexagonoOrigen"].apply(lambda x: h3pol(x))
distancia_caminando_hexagono.rename(columns={"HexagonoOrigen": "h3_code"}, inplace=True)
mapa0 = folium.Map(coordenadasmapa, tiles="OpenStreetMap", zoom_start=14)
mapa0hecho = choropleth_map(distancia_caminando_hexagono, "CaminataProporcional", initial_map=mapa0)
mapa0.save(directoriosalida + "MapaCaminataHexagonos.html")



rutasseparadas = CargaT["Ruta"].tolist()
rutasseparadas = [sublist[1:-1] for sublist in rutasseparadas]


nodosvelocidades = pd.read_excel(direccionnodosvelocidad)
pvelocidad = pd.read_excel(direccionperfilv)
nodosvelocidades = pd.merge(nodosvelocidades, pvelocidad, on=["factorx"], how="left")
nodosvelocidades["Velocidad"] = nodosvelocidades["Velocidad"] * nodosvelocidades["correccion"]

dictnodos = nodosvelocidades.set_index(['IDtramo', "HORA"])
dictnodos = dictnodos[["Velocidad", "Distancia"]].to_dict(orient="index")
perfilvelocidades = pd.merge(nodosvelocidades, SHP, left_on="idlinea", right_on="Lineas", how="left")
perfilvelocidades["peso"] = perfilvelocidades["Distancia"] / perfilvelocidades["Longitud"]
perfilvelocidades["Velocidad"] = perfilvelocidades["Velocidad"] * perfilvelocidades["peso"]
perfilvelocidades = perfilvelocidades.groupby(["HORA", "Lineas"])["Velocidad"].sum().reset_index()
perfilvelocidades.to_excel(directoriosalida + "PerfilVelocidades.xlsx")


horasrutas = CargaT["HORA"].tolist()

listatiempos = []
for i, ruta in enumerate(rutasseparadas):
    rutatemporal = []
    a = 0
    hint = horasrutas[i]
    for p in ruta:
        rutatemporal.append(a)
        distancia = dictnodos[p,hint]["Distancia"] / 1000
        velocidad = dictnodos[p,hint]["Velocidad"]
        tiempo = (distancia / velocidad) * 60
        a += tiempo
    listatiempos.append(rutatemporal)
CargaT["RutaT"] = listatiempos


def minutos_a_timedelta(minutos):
    return [pd.Timedelta(f'{m} minutes') for m in minutos]

def ponerfecha(row):
    import datetime
    a = datetime.datetime(2023, 9, 10, row['HORA'], row['MINUTO'])
    return a

CargaT['FECHATRX'] = CargaT.apply(lambda row: ponerfecha(row), axis=1)
CargaT["COTASUPFECHAFIN"] = CargaT["FECHATRX"] + pd.to_timedelta(1, unit='h')

CargaT['fechas'] = CargaT.apply(lambda row: [row['FECHATRX'] + td for td in minutos_a_timedelta(row['RutaT'])], axis=1)
CargaT.to_csv(direccioncargat)
intervalos = pd.interval_range(start=pd.to_datetime('2023-09-10 05:00:00'), end=pd.to_datetime('2023-09-10 23:00:00'),
                               freq=str(intervalofrecuencia) + "T")
listadedf = []
for inte in intervalos:
    inicio = inte.left
    fin = inte.right
    trxenintervalo = CargaT[(CargaT["FECHATRX"] <= fin)]
    trxenintervalo = trxenintervalo[(trxenintervalo["COTASUPFECHAFIN"] >= inicio)]
    trxenintervalo = trxenintervalo.reset_index(drop=True)
    rutasseparadas = trxenintervalo["Ruta"].tolist()
    rutasseparadas = [sublist[1:-1] for sublist in rutasseparadas]
    rutasseparadast = trxenintervalo["fechas"].tolist()
    nodos = list(set([nodo for ruta in rutasseparadas for nodo in ruta]))
    dictinterno = {nodo: 0 for nodo in nodos}
    for i, (r, rt) in enumerate(zip(rutasseparadas, rutasseparadast)):
        for p, (n, t) in enumerate(zip(r, rt)):
            if (t >= inicio) and (t <= fin):
                dictinterno[n] += trxenintervalo["CargaProporcional"].iloc[i]
            else:
                continue
    listadedf.append(dictinterno)

FrecuenciaAnalitica = pd.DataFrame(listadedf).T
FrecuenciaAnalitica = FrecuenciaAnalitica.set_axis(intervalos, axis=1, inplace=False).stack().reset_index()
FrecuenciaAnalitica['inicio'] = FrecuenciaAnalitica.apply(lambda row: row["level_1"].left, axis=1)
FrecuenciaAnalitica['fin'] = FrecuenciaAnalitica.apply(lambda row: row["level_1"].right, axis=1)
FrecuenciaAnalitica["Separar"] = FrecuenciaAnalitica["level_0"]
FrecuenciaAnalitica[["LineaUtilizada", "Orden"]] = FrecuenciaAnalitica["Separar"].str.split("-", n=1, expand=True)
FrecuenciaAnalitica["Orden"] = FrecuenciaAnalitica["Orden"].astype(float)
FrecuenciaAnalitica = FrecuenciaAnalitica[FrecuenciaAnalitica["LineaUtilizada"].isin(Lineas)]
FrecuenciaAnalitica.sort_values(by=["LineaUtilizada", "Orden"], inplace=True)
FrecuenciaAnalitica["inicio"] = pd.to_datetime(FrecuenciaAnalitica["inicio"])
FrecuenciaAnalitica.rename(columns={0: "0", "level_0": "index"}, inplace=True)
FrecuenciaAnalitica = FrecuenciaAnalitica[["inicio", "fin", "LineaUtilizada", "Orden", "0", "index"]]
FrecuenciaAnalitica.to_excel(directoriosalida + "FrecuenciasDesagregadas.xlsx", index=False)

FA_pivot = FrecuenciaAnalitica.pivot(index='index', columns='inicio', values="0")
FA_pivot.reset_index(inplace=True, drop=False)
FA_pivot["Separar"] = FA_pivot["index"]
FA_pivot[["LineaUtilizada", "Orden"]] = FA_pivot["Separar"].str.split("-", n=1, expand=True)
FA_pivot["Orden"] = FA_pivot["Orden"].astype(float)
FA_pivot.sort_values(by=["LineaUtilizada", "Orden"], inplace=True)
FA_pivot.to_excel(directoriosalida + "FrecuenciaPivoteadas.xlsx", index=False)

grupoporhora = FrecuenciaAnalitica.groupby([pd.Grouper(key='inicio', freq=str(longitudturno) + "H", offset="5h"), "LineaUtilizada"])[
    "0"].max().reset_index()
grupoporhora = pd.merge(grupoporhora, SHP, left_on="LineaUtilizada", right_on="Lineas", how="left")
grupoporhora = grupoporhora[['inicio', 'Lineas', '0', 'Longitud']]
grupoporhora["HORA"] = grupoporhora["inicio"].dt.hour
grupoporhora = pd.merge(grupoporhora, perfilvelocidades, left_on=["Lineas", "HORA"], right_on=["Lineas", "HORA"], how="left")
grupoporhora["Tiempoporvuelta"] = (grupoporhora["Longitud"] / (grupoporhora["Velocidad"] * 1000))
grupoporhora["Coches"] = (grupoporhora["0"] * grupoporhora["Tiempoporvuelta"] / capacidadcoche).apply(np.ceil)
grupoporhora["NroVehiculosFicticios"] = grupoporhora["0"] * grupoporhora["Tiempoporvuelta"] / capacidadcoche
grupoporhora["FrecuenciaImplicita"] = (grupoporhora["Tiempoporvuelta"] / grupoporhora["Coches"]) * 60
grupoporhora["kmrecorridos"] = (longitudturno / grupoporhora["Tiempoporvuelta"]) * grupoporhora["Coches"] * grupoporhora[
    "Longitud"] / 1000
cochesporturno = grupoporhora.groupby(["inicio"])["Coches"].sum().reset_index()
cochesporturno.rename(columns={"Coches": "TotalCochesTurno"}, inplace=True)
distanciarecorridoporlinea = grupoporhora.groupby(["Lineas"])["kmrecorridos"].sum().reset_index()
grupoporhora["Utilización"] = grupoporhora["0"] * 100 / (capacidadcoche * grupoporhora["Coches"] / grupoporhora["Tiempoporvuelta"])


Ocio = FrecuenciaAnalitica.groupby(["LineaUtilizada", pd.Grouper(key='inicio', freq=str(longitudturno) + "H")])["0"].max().reset_index()
Ocio["HORA"] = Ocio["inicio"].dt.hour
Ocio.rename(columns={"0": "CargaMaxima"}, inplace=True)
FrecuenciaAnalitica["HORA"] = FrecuenciaAnalitica["inicio"].dt.hour
Ocio = pd.merge(FrecuenciaAnalitica, Ocio, on=["LineaUtilizada", "HORA"], how="left")
Ocio["Ocupacion"] = Ocio["0"] / Ocio["CargaMaxima"]
OcioPromedio = Ocio.groupby(["LineaUtilizada", pd.Grouper(key='inicio_x', freq=str(longitudturno) + "H", offset="5h")])[
    "Ocupacion"].mean().reset_index()
OcioPromedio["Ocupacion"] = OcioPromedio["Ocupacion"] * 100


pmporhora = pd.read_csv(direccionpmporhora)
pmporhora["MINUTO"] = 0
pmporhora['FECHATRX'] = pmporhora.apply(lambda row: ponerfecha(row), axis=1)
pmporhora[['LineaUtilizada', 'Orden']] = pmporhora["Nodo"].str.split("-", n=1, expand=True)
pmax = pmporhora.groupby(["LineaUtilizada", pd.Grouper(key='FECHATRX', freq=str(longitudturno) + "H", offset = "5H")])["Personas"].max().reset_index()
pmax = pmax.rename(columns={"Personas": "MaximoPersonas"})
Homogeneidad = pd.merge(pmporhora, pmax, on=["LineaUtilizada", "FECHATRX"], how="left")
Homogeneidad["Proporcion"] = Homogeneidad["Personas"] / Homogeneidad["MaximoPersonas"]
HomogeneidadTurno = Homogeneidad.groupby(["LineaUtilizada", pd.Grouper(key='FECHATRX', freq=str(longitudturno) + "H", offset = "5H")])["Proporcion"].mean().reset_index()
HomogeneidadTurno.rename(columns={"Proporcion": "Homogeneidad"}, inplace=True)
HomogeneidadTurno["Homogeneidad"] = HomogeneidadTurno["Homogeneidad"]*100
HomogeneidadPromedio = HomogeneidadTurno.Homogeneidad.mean()


Cargaporlinea = CargaT.groupby(["LineaUtilizada", pd.Grouper(key='FECHATRX', freq=str(longitudturno) + "H", offset="5h")])[
    "CargaProporcional"].sum().reset_index()
Cargaporlinea = Cargaporlinea[Cargaporlinea["LineaUtilizada"].isin(Lineas)]

IPK = pd.merge(Cargaporlinea, grupoporhora, left_on=["LineaUtilizada", "FECHATRX"], right_on=["Lineas", "inicio"],
               how="inner")
IPK["IPK"] = IPK["CargaProporcional"] / IPK["kmrecorridos"]
IPK = IPK[["LineaUtilizada", "inicio", "CargaProporcional", "kmrecorridos", "Coches", "NroVehiculosFicticios", "FrecuenciaImplicita", "IPK", "Utilización"]]
IPK = pd.merge(IPK, cochesporturno, left_on="inicio", right_on="inicio", how="left")
IPK["CostoRelativoTurno"] = 0
IPK.replace([np.inf, -np.inf], np.nan, inplace=True)
IPK = IPK.dropna()
IPK.drop(columns=["TotalCochesTurno"], inplace=True)
IPK = pd.merge(IPK, OcioPromedio, left_on=["LineaUtilizada", "inicio"], right_on=["LineaUtilizada", "inicio_x"],
               how="left")

IPK.rename(columns={"inicio": "InicioTurno", "LineaUtilizada": "Linea", "CargaProporcional": "TRXTurno",
                    "kmrecorridos": "DistanciaRecorrida", "Coches": "NroVehiculos",
                    "FrecuenciaImplicita": "Frecuencia"}, inplace=True)
IPK = pd.merge(IPK, SHP, left_on="Linea", right_on="Lineas", how="left")
IPK = pd.merge(IPK, HomogeneidadTurno, left_on=["Linea", "InicioTurno"], right_on=["LineaUtilizada", "FECHATRX"], how="left")
IPK = IPK[["Linea", "InicioTurno", "TRXTurno", "NroVehiculos", "NroVehiculosFicticios", "DistanciaRecorrida", "Longitud", "Frecuencia", "Ocupacion", "Utilización",  "Homogeneidad",
           "CostoRelativoTurno", "IPK"]]

IPKpromedio = IPK.TRXTurno.sum() / IPK.DistanciaRecorrida.sum()
NroCochesMinimos = IPK.groupby(["InicioTurno"])["NroVehiculos"].sum().max()
FrecuenciaPromedio = IPK.groupby(["InicioTurno"])["Frecuencia"].mean().mean()
DistanciaDiaria = IPK["DistanciaRecorrida"].sum()
OcupacionPromedio = IPK["Ocupacion"].mean()
UtilizacionPromedio = IPK["Utilización"].mean()
personasasignadas = IPK["TRXTurno"].sum()
IPK.to_excel(directoriosalida + "ReporteCompletoIPK.xlsx", index=False)

#Reporte por Linea
ReporteLinea = IPK.groupby(["Linea"]).agg({"TRXTurno": "sum", "NroVehiculos": "mean", "DistanciaRecorrida": "sum", "Frecuencia": "mean", "Ocupacion": "mean", "Utilización": "mean", "Homogeneidad": "mean"}).reset_index()
ReporteLinea["IPK"] = ReporteLinea["TRXTurno"] / ReporteLinea["DistanciaRecorrida"]
ReporteLinea.rename(columns={"TRXTurno": "TRX", "NroVehiculos": "NroVehiculosPromedio", "Frecuencia": "FrecuenciaPromedio", "Ocupacion": "OcupacionPromedio", "Utilización": "UtilizacionPromedio", "Homogeneidad": "HomogeneidadPromedio"}, inplace=True)


TablaTuning = IPK[["Linea", "InicioTurno", "Frecuencia"]]
TablaTuning.to_excel(directoriosalida + "TablaTuning.xlsx", index=False)

# Medias

MediaPM = FrecuenciaAnalitica.groupby(["index"])["0"].mean().reset_index()
MediaPM["Separar"] = MediaPM["index"]
MediaPM[["LineaUtilizada", "Orden"]] = MediaPM["Separar"].str.split("-", n=1, expand=True)
MediaPM["Orden"] = MediaPM["Orden"].astype(float)
MediaPM.sort_values(by=["LineaUtilizada", "Orden"], inplace=True)
MediaPM.rename(columns={"0": "Media"}, inplace=True)

DesvioPM = FrecuenciaAnalitica.groupby(["index"])["0"].std().reset_index()
DesvioPM["Separar"] = DesvioPM["index"]
DesvioPM[["LineaUtilizada", "Orden"]] = DesvioPM["Separar"].str.split("-", n=1, expand=True)
DesvioPM["Orden"] = DesvioPM["Orden"].astype(float)
DesvioPM.rename(columns={"0": "Desvio"}, inplace=True)
DesvioPM.sort_values(by=["LineaUtilizada", "Orden"], inplace=True)

MediaPM = pd.merge(MediaPM, DesvioPM, on=["index", "LineaUtilizada", "Orden"], how="inner")
MediaPM.to_excel(directoriosalida + "MediaPM.xlsx", index=False)


# Tiempos de viaje

RutasSimples = RutasCompletas[(RutasCompletas["NroLineas"] == 1) & (RutasCompletas["LineaUtilizada"] != "Ninguna")][["Nodo Origen", "Nodo Destino", "DistanciaCaminando", "LineaUtilizada"]]
RutasCombinadas = RutasCompletas[RutasCompletas["NroLineas"] > 1][["Nodo Origen", "Nodo Destino", "DistanciaCaminando", "LineaUtilizada"]]

RutasSimples["TiempoCaminando"] = (RutasSimples["DistanciaCaminando"] / velocidadcaminando) * (60/1000)
TiempoDeViaje = CargaT[['Nodo Origen', 'Nodo Destino', 'LineaUtilizada', "RutaT"]]
TiempoDeViaje.drop_duplicates(subset=["Nodo Origen", "Nodo Destino", "LineaUtilizada"], inplace=True)
TiempoDeViaje["TiempoDeViaje"] = TiempoDeViaje["RutaT"].apply(lambda x: x[-1])
TiempoDeViaje["TiempoDeViaje"] = TiempoDeViaje["TiempoDeViaje"].astype(float)
TiempoDeViaje = TiempoDeViaje[["Nodo Origen", "Nodo Destino", "LineaUtilizada", "TiempoDeViaje"]]

TiempoDeEspera = grupoporhora[["inicio", "Lineas", "FrecuenciaImplicita"]]
TiempoDeEspera["TiempoDeEspera"] = TiempoDeEspera["FrecuenciaImplicita"]/2

tiemposdecorte = TiempoDeEspera["inicio"].unique()
TiempoDeEspera["inicio"] = TiempoDeEspera["inicio"].astype('datetime64[ns]')
fechaextra = tiemposdecorte[-1] + pd.to_timedelta(longitudturno+0.5, unit='h')
fechaextra = np.array(fechaextra)
tiemposdecorte = np.append(tiemposdecorte, fechaextra).astype('datetime64[ns]')
intervalos = pd.IntervalIndex.from_breaks(breaks=tiemposdecorte, closed='right')
CargaT['intervalo'] = pd.cut(CargaT['FECHATRX'], bins=intervalos, labels=range(len(intervalos)))
CargaT['inicio'] = CargaT["intervalo"].apply(lambda x: x.left)
CargaTparatiempos = CargaT[["Nodo Origen", "Nodo Destino", "LineaUtilizada", "inicio", "CargaProporcional"]]
CargaTparatiempos["inicio"] = CargaTparatiempos["inicio"].astype('datetime64[ns]')
CargaTparatiempos = pd.merge(CargaTparatiempos, TiempoDeEspera, left_on=["inicio", "LineaUtilizada"], right_on=["inicio", "Lineas"], how="left")
CargaTparatiempos = pd.merge(CargaTparatiempos, TiempoDeViaje, left_on=["Nodo Origen", "Nodo Destino", "LineaUtilizada"], right_on=["Nodo Origen", "Nodo Destino", "LineaUtilizada"], how="left")
CargaTparatiempos = pd.merge(CargaTparatiempos, RutasSimples, left_on=["Nodo Origen", "Nodo Destino", "LineaUtilizada"], right_on=["Nodo Origen", "Nodo Destino", "LineaUtilizada"], how="left")
CargaTparatiempos = CargaTparatiempos[['Nodo Origen', 'Nodo Destino', 'LineaUtilizada', 'inicio', 'TiempoDeEspera', 'TiempoDeViaje', 'TiempoCaminando', 'CargaProporcional']]
CargaTparatiempos["TiempoTotal"] = CargaTparatiempos["TiempoDeEspera"] + CargaTparatiempos["TiempoDeViaje"] + CargaTparatiempos["TiempoCaminando"]
CargaSumada = CargaTparatiempos["CargaProporcional"].sum()
CargaTparatiempos["TiempoPonderado"] = CargaTparatiempos["TiempoTotal"] * CargaTparatiempos["CargaProporcional"] / CargaSumada
CargaTparatiempos["TiempoPonderadoC"] = CargaTparatiempos["TiempoCaminando"] * CargaTparatiempos["CargaProporcional"] / CargaSumada
CaminataPromedio = CargaTparatiempos.TiempoPonderadoC.sum()
TiempoPromedio = CargaTparatiempos.TiempoPonderado.sum()
CargaTparatiempos.to_csv(directoriosalida + "TiemposDesagregados.csv", index=False)





end = time.time()
tiempoejecucion = (end - start) / 60
print("Tiempo de ejecución: " + str(tiempoejecucion) + " minutos")

book = load_workbook(direccionmetricas)
writer = pd.ExcelWriter(direccionmetricas, engine = 'openpyxl')
writer.book = book

metricas = pd.DataFrame({
    'IPKp': [IPKpromedio], "PersonasAsignadas" : [personasasignadas],
    'NroCochesMinimo': [NroCochesMinimos], "FrecuenciaPromedio" : [FrecuenciaPromedio],
    'OcupacionPromedio%': [OcupacionPromedio],'UtilizacionPromedio%': [UtilizacionPromedio], "Homogeneidad%": [HomogeneidadPromedio], "KMTotales" : [DistanciaDiaria], "TiempoPromedio" : [TiempoPromedio], "TiempoCaminataPromedio" : [CaminataPromedio],  "TiempoEjecucion": [tiempoejecucion]
}).T

metricas.to_excel(writer, sheet_name = 'MetricasFrecuencia', index=True)
writer.close()
