import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
import pandas as pd
import numpy as np
import folium
import shapely
import warnings
from fiona.drvsupport import supported_drivers
from dividirlinea4 import tramos
import networkx as nx
import math
import time
from shapely.ops import transform
from pyproj import Transformer
from string import digits
import ast
import re
from tqdm import tqdm
import h3
from shapely import Polygon
from geojson import Feature, Point, FeatureCollection
import json
import matplotlib
from branca.colormap import linear

warnings.filterwarnings('ignore')
supported_drivers['LIBKML'] = 'rw'

# Funciones

def generate_color_scale(min_value, max_value):
    colormap = linear.RdYlGn_11.scale(min_value, max_value)
    return colormap

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
            feature = Feature(geometry=geometry_for_row, id=row["h3_code"], properties={column_name: row[column_name], "h3_code": row["h3_code"]})
            list_features.append(feature)
        except:
            print("An exception occurred for hex " + row["h3_code"])

    feat_collection = FeatureCollection(list_features)
    geojson_result = json.dumps(feat_collection)
    return geojson_result

def h3pol(h3_id):
    strings_separados = h3_id.split("-")
    h3code = strings_separados[0]
    return Polygon(h3.h3_to_geo_boundary(h3code, geo_json=True))


def get_color(custom_cm, val, vmin, vmax):
    return matplotlib.colors.to_hex(custom_cm((val - vmin) / (vmax - vmin)))


def choropleth_map(df_aggreg, column_name="value", border_color='black', fill_opacity=0.1, color_map_name="Greens",
                   initial_map=None):
    """
    Creates choropleth maps given the aggregated data. initial_map can be an existing map to draw on top of.
    """
    # colormap
    min_value = df_aggreg[column_name].min()*2
    max_value = df_aggreg[column_name].max()*0.4
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
            'fillColor': get_color(custom_cm, feature['properties'][column_name], vmin=min_value, vmax=max_value),
            'color': border_color,
            'weight': 1,
            'fillOpacity': fill_opacity
        },
        name=name_layer, tooltip=folium.features.GeoJsonTooltip(fields=[column_name, "h3_code"], aliases=[column_name, "h3_code"])
    ).add_to(initial_map)
    return initial_map


# Parametros
distanciatolerancia = 1000  # distancia para conectar las paradas magicas a las geolocalizaciones
distanciatoleranciacombinacion = 200  # distancia para que dos paradas magicas se puedan conectar
costocombinacion = 30000  # big M para el costo de la combinacion
costocaminar = 50000  # big M para el costo de caminar
longitudtramo = 100  # distancia entre paradas magicas
ncaminosalternativos = 3  # nro de caminos alternativos a considerar
penalizacioncaminata = 2  # penalizacion por caminata
degeneracion = 0.25 # degeneracion admitida en la distancia total para un camino alternativo
exponente_suavizacion = 1  # exponente de suavizacion de la funcion de probabilidad
coordenadasmapa = [-33.74595770280119, -61.96865467712429]
lineasexcluidas = ["Cementerio", "Nocturno Sur", "Nocturno Norte", "LA BOCA", "121AB"]  # lineas a excluir
Validacion = False

# Direcciones de entrada
direccionlineas = "SHP Files/Lineas/LineasNuevas/RecorridosF.shp"  # direccion de las lineas
direccionod = "Output/ODexpVT.csv" # direccion de la OD
direcciondfseparaciones = "Rawdata/SeparacionParadas.xlsx"  # direccion del df de separaciones
direccionseparacionkml = "SHP Files/Velocidades/ZonaParadasVT.kml"
direccionfactorecorreccion = "Output/estadisticosodporlineaMSN.xlsx"

# Direcciones de salida
codigo = "VT_GTRX35_LE"
directoriosalida = "Output/Pruebas Nuevos Sistemas/Legado/" + codigo

# cargo la MOD y las paradas
start = time.time()
od = pd.read_csv(direccionod)  # cargo la OD
od["LineaCombinacion"].fillna(0, inplace=True)  # relleno los nulos de la columna LineaCombinacion
# od["LineaCombinacion"] = od["LineaCombinacion"].astype(int)  # convierto la columna LineaCombinacion a int

od['h3_code_ida'] = od.apply(lambda row: h3.geo_to_h3(row['latitude_Origen'], row['longitude_Origen'], 9), axis=1)
od['h3_code_vuelta'] = od.apply(lambda row: h3.geo_to_h3(row['latitude_Destino'], row['longitude_Destino'], 9), axis=1)
od["VueltaAux"] = od["h3_code_vuelta"]

od["h3_code_ida"] = od["h3_code_ida"] + "-" + od["idlinea"].astype(str)
od["h3_code_vuelta"] = od["h3_code_vuelta"] + "-" + od["idlinea"].astype(str)
condicion = [od["LineaCombinacion"] != 0]
valor = [od["VueltaAux"] + "-" + od["LineaCombinacion"].astype(str)]
od["h3_code_vuelta"] = np.select(condicion, valor, default=od["h3_code_vuelta"])
LineasViejas = od["idlinea"].unique().tolist()


clustersagrupados = od.groupby(["h3_code_ida", "h3_code_vuelta"]).size().to_frame("Cuenta").reset_index()  # Agrupo
clustert = od.copy(deep=True)  # me llevo una copia de paradas
clustert["FECHATRX"] = pd.to_datetime(clustert["FECHATRX"])  # transformo la columna en una columna de tipo datetime
clustersagrupadost = clustert.groupby(["h3_code_ida", "h3_code_vuelta", pd.Grouper(key='FECHATRX', freq='15min')]).size().to_frame(
    "Cuenta").reset_index()  # agrupo por paradaida, paradavuelta y por intervalo de 15 minutos



clustersagrupados["geometry_x"] = clustersagrupados["h3_code_ida"].apply(lambda x: h3pol(x))
clustersagrupados["geometry_y"] = clustersagrupados["h3_code_vuelta"].apply(lambda x: h3pol(x))
shpclusters = pd.concat([clustersagrupados["h3_code_ida"], clustersagrupados["h3_code_vuelta"]], axis=0).reset_index()
shpclusters = shpclusters.drop_duplicates().reset_index()
shpclusters.rename(columns={0: "h3_code"}, inplace=True)
shpclusters["geometry"] = shpclusters["h3_code"].apply(lambda x: h3pol(x))
shpclusters = gpd.GeoDataFrame(shpclusters, geometry="geometry")
shpclusters = shpclusters.set_crs(epsg=4326)
clustergjson = shpclusters[['h3_code']]
clustergjson["Valor"] = 10
clustergjson[['h3_code', 'idlinea']] = clustergjson['h3_code'].str.split('-', expand=True)
clustergjson = clustergjson.drop_duplicates().reset_index(drop=True)
mapa0 = folium.Map(coordenadasmapa, tiles="OpenStreetMap", zoom_start=14)
mapa0hecho = choropleth_map(clustergjson, "Valor", initial_map=mapa0)
mapa0.save(directoriosalida + "MapaHexagonos.html")
shpclusters["latitude"] = shpclusters["geometry"].apply(lambda x: x.centroid.y)
shpclusters["longitude"] = shpclusters["geometry"].apply(lambda x: x.centroid.x)
shpclusters = shpclusters[["h3_code", "latitude", "longitude"]]
shpclusters.drop_duplicates(inplace=True)
shpclusters.reset_index(inplace=True, drop=True)



# cargar el shape de la linea y lo preproceso
SHP = gpd.read_file(direccionlineas, crs = 22185)
SHP = SHP.rename(columns={'Linea': 'Lineas'})
SHP = SHP.to_crs(epsg=22185)
Lineas = SHP["Lineas"].unique()
SHP["Linea"] = SHP["Lineas"]
df = SHP.set_geometry(
    SHP.geometry.map(
        lambda linestring: shapely.ops.transform(lambda x, y, *_: (x, y), linestring)
    ))
shpnodos = tramos(df, longitudtramo, direcciondfseparaciones, direccionseparacionkml)  # creo los tramos con la funcion tramos
shpnodos.to_excel(directoriosalida + "nodos.xlsx")

# creo el grafo

G = nx.DiGraph()
# nodos
for l in Lineas:
    # primero agrego el nodo
    dfaux = shpnodos[shpnodos["Lineas"] == l].reset_index()
    nodosaux = dfaux["IDtramo"]
    G.add_nodes_from(nodosaux)
    # segundo tengo que agregar los atributos de cada nodo y en el mismo lazo agrego las aristas
    for i in range(len(nodosaux)):
        G.nodes[nodosaux[i]]["geometry"] = dfaux.at[i, "inicio"]
        G.nodes[nodosaux[i]]["Lineas"] = l
        dist = dfaux.at[i, "Distancia"]
        try:
            G.add_edge(nodosaux[i], nodosaux[i + 1], distanciareal=dist, distanciaCompensada=dist, distanciaCaminata=0,
                       distanciacombinacion=0)
        except KeyError:
            continue

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

# Ahora lo que necesito es calcular para cada parada cuales son los puntos de las lineas nuevas a menos de la tolerancia metros
puntoscluster = gpd.GeoDataFrame(shpclusters,
                                 geometry=gpd.points_from_xy(shpclusters['longitude'], shpclusters['latitude']),
                                 crs="4326").to_crs(epsg=22185)
puntoscluster = puntoscluster.reset_index()

# mapa = folium.Map(coordenadasmapa, tiles="OpenStreetMap", zoom_start=14)
# for l in LineasViejas:
#     fg = folium.FeatureGroup("LINEA" + " " + str(l))
#     puntosparadasinterno = puntosparadas[puntosparadas["Lineas"] == l].reset_index().to_crs(epsg=4326)
#     for j, nmg in enumerate(puntosparadasinterno["geometry"]):
#         fg.add_child(folium.CircleMarker(
#                 location=(nmg.y, nmg.x),
#                 radius=4,
#                 fill=True,
#                 color="#12bc8e", popup=puntosparadasinterno.at[j, "Key_y"]))
#     fg.add_to(mapa)
# folium.LayerControl().add_to(mapa)
# loc = "Paradas Sistema Viejo"
# title_html = '''
#                     <h3 align="center" style="font-size:16px"><b>{}</b></h3>
#                    '''.format(loc)
# mapa.get_root().html.add_child(folium.Element(title_html))
# mapa.save("Venado/Output/ParadasRealesSistemaViejo.html") # Guardo.

puntosnodos = shpnodos.set_geometry("inicio").to_crs(epsg=22185)
ClusterNodos = pd.DataFrame()
ClusterNodos['NodosCercanos'] = None
ClusterNodos['NodosCercanos'] = ClusterNodos['NodosCercanos'].astype('object')
transformer = Transformer.from_crs("epsg:4326",
                                   "epsg:22185")  # para transformar los puntos en long/lat del grafo en una proyeccion que me permita medir metros

for p in range(len(puntoscluster)):
    A = []  # lista de los nodos cercanos a la parada p
    ClusterNodos.at[p, "IDCluster"] = puntoscluster.at[p, "h3_code"]  # agrego el ID de la parada
    for l in Lineas:
        puntosnodosinterno = puntosnodos[puntosnodos["Lineas"] == l].reset_index()  # filtro
        cercanos = {}  # creo un diccionario vacio
        for n in range(len(puntosnodosinterno)):
            distancia = puntoscluster.at[p, "geometry"].distance(puntosnodosinterno.at[n, "inicio"]) * math.sqrt(
                2)  # mido la distancia compensada por raiz de 2
            if distancia <= distanciatolerancia:
                cercanos[puntosnodosinterno.at[
                    n, "IDtramo"]] = distancia  # si es menor a la distancia de tolerancia la agrego al diccionario
            else:
                pass
        cercanos = dict(
            sorted(cercanos.items(), key=lambda x: x[1]))  # ordeno el diccionario por distancia de menor a mayor
        B = list(cercanos.keys())  # creo una lista con los nodos cercanos
        Bbis = []  # lista vacia auxiliar
        try:
            Bbis.append(B[0])  # agrego el mas cercano
            ordeninicial = int(re.findall(r'\d+', B[0])[0])  # extraigo el orden del primer nodo
            for inx, mg in enumerate(B):  # recorro la lista de nodos cercanos
                orden = int(re.findall(r'\d+', mg)[0])  # extraigo el orden del nodo
                if abs(orden - ordeninicial) >= 10:
                    Bbis.append(
                        mg)  # si la diferencia entre los dos ordenes es mayor a 10 agrego el nodo a la lista auxiliar
                    break
                else:
                    pass
            A = A + Bbis
            B = []
            Bbis = []
        except IndexError:
            pass
    A = set(A)
    A = list(A)
    ClusterNodos.at[p, "NodosCercanos"] = A

ParadasExcluidas = ClusterNodos[ClusterNodos["NodosCercanos"].str.len() == 0]
ParadasNodos = ClusterNodos[ClusterNodos[
                                'NodosCercanos'].str.len() > 0]  # borro las paradas que no tienen nodos cercanos idealmente no deberia haber
ParadasNodos = ParadasNodos.reset_index()
G2 = G.copy()  # creo una copia del grafo original
idclusternp = np.array(ParadasNodos["IDCluster"])

for id in idclusternp:
    # tengo que agregar el nodo
    nuevonodo = id
    nodosaconectar = ParadasNodos[ParadasNodos["IDCluster"] == nuevonodo]["NodosCercanos"].tolist()
    geometria = shpclusters[shpclusters["h3_code"] == nuevonodo]["geometry"].tolist()[0]
    G2.add_nodes_from([nuevonodo])
    G2.nodes[nuevonodo]["geometry"] = geometria
    G2.nodes[nuevonodo]["Lineas"] = "Parada"
    for nd in nodosaconectar[0]:
        # para cada posible conexion tengo que agregar una arista
        geometriaparada = puntoscluster[puntoscluster["h3_code"] == nuevonodo]["geometry"].tolist()[0]
        geometrianodo = puntosnodos[puntosnodos["IDtramo"] == nd]["inicio"].tolist()[0]
        distancia = geometriaparada.distance(geometrianodo) * penalizacioncaminata
        distanciatrue = (distancia / penalizacioncaminata) * 1.27
        G2.add_edge(nuevonodo, nd, distanciareal=distanciatrue, distanciaCompensada=costocaminar + distancia,
                    distanciaCaminata=distanciatrue, distanciacombinacion=0)
        G2.add_edge(nd, nuevonodo, distanciareal=distanciatrue, distanciaCompensada=costocaminar + distancia,
                    distanciaCaminata=distanciatrue, distanciacombinacion=0)  # por si es un destino

# ahora agrego la posibilidad de realizar combinaciones
G_nocomb = G2.copy()  # creo una copia del grafo original

for l in Lineas:
    # para cada linea
    naux = np.array(shpnodos[shpnodos["Lineas"] == l]["IDtramo"])  # los nodos que quiero comparar
    naux2 = np.array(shpnodos[shpnodos["Lineas"] != l]["IDtramo"])  # el complemento
    for n1 in naux:
        geometrianodo1 = transform(transformer.transform, G2.nodes[n1]["geometry"])
        # para cada nodo de la linea
        for n2 in naux2:
            # los nodos de las otras lineas
            geometrianodo2 = transform(transformer.transform, G2.nodes[n2]["geometry"])
            distancia = geometrianodo1.distance(geometrianodo2)
            if distancia < distanciatoleranciacombinacion:
                G2.add_edge(n1, n2, distanciareal=distancia, distanciaCompensada=distancia + costocombinacion,
                            distanciaCaminata=distancia, distanciacombinacion=costocombinacion)
            else:
                continue

paradasorigen = np.array(clustersagrupados["h3_code_ida"])
paradasdestino = np.array(clustersagrupados["h3_code_vuelta"])

Aparallenar = []  # lista para luego generar el df
for po, pde in tqdm(zip(paradasorigen, paradasdestino), total=len(paradasorigen)):
    try:
        marcador = False
        s = nx.shortest_path(G_nocomb, source=po, target=pde, weight="distanciaCompensada")
        distanciacompensada = nx.path_weight(G_nocomb, path=s, weight="distanciaCompensada")
        if distanciacompensada >= 3*costocaminar:
            raise nx.NetworkXNoPath
        else:
            l = nx.path_weight(G_nocomb, path=s, weight="distanciareal")
            c = nx.path_weight(G_nocomb, path=s, weight="distanciaCaminata")
            comb = nx.path_weight(G_nocomb, path=s, weight="distanciacombinacion")
            dictinterno = {"Nodo Origen": po,
                           "Nodo Destino": pde,
                           "Ruta": s,
                           "Distancia Total": l,
                           "DistanciaCaminando": c}
            Aparallenar += [dictinterno]
            if po == pde:
                pass
            else:
                lineaoriginal = [G_nocomb.nodes[s[1]]["Lineas"]]
                caminosrestantes = ncaminosalternativos - 1
                while caminosrestantes > 0:
                    marcador = True
                    nodossub = [nodo for nodo, datos in G_nocomb.nodes(data=True) if datos["Lineas"] not in lineaoriginal]
                    subG = G_nocomb.subgraph(nodossub)
                    s = nx.shortest_path(subG, source=po, target=pde, weight="distanciaCompensada")
                    distanciacompensada = nx.path_weight(subG, path=s, weight="distanciaCompensada")
                    if distanciacompensada >= 3*costocaminar:
                        raise nx.NetworkXNoPath
                    else:
                        lalt = nx.path_weight(subG, path=s, weight="distanciareal")
                        calt = nx.path_weight(subG, path=s, weight="distanciaCaminata")
                        combalt = nx.path_weight(subG, path=s, weight="distanciacombinacion")
                        condicion1 = (lalt - l <= l * degeneracion)  # no tiene que empeorar mucho la solucion
                        condicion2 = (comb >= combalt) or (
                                    (comb == 0) and (combalt == 0))  # control combinaciones
                        condicion3 = (
                                    lalt > l)  # si la distancia mejora es que esta haciendo cosas raras con tal de encontrar una solucion factible
                        condicion4 = (calt <= c / degeneracion)
                        if condicion1 and condicion3 and condicion2 and condicion4:
                            dictinterno = {"Nodo Origen": po,
                                           "Nodo Destino": pde,
                                           "Ruta": s,
                                           "Distancia Total": lalt,
                                           "DistanciaCaminando": calt}
                            Aparallenar += [dictinterno]
                            lineaoriginal.append(subG.nodes[s[1]]["Lineas"])
                            caminosrestantes = caminosrestantes - 1
                        else:
                            break
    except nx.NodeNotFound:
        continue
    except nx.NetworkXNoPath:
        try:
            if marcador:
                continue
            else:
                s = nx.shortest_path(G2, source=po, target=pde, weight="distanciaCompensada")
                distanciacompensada = nx.path_weight(G2, path=s, weight="distanciaCompensada")
                if distanciacompensada >= 3*costocaminar:
                    continue
                else:
                    l = nx.path_weight(G2, path=s, weight="distanciareal")
                    c = nx.path_weight(G2, path=s, weight="distanciaCaminata")
                    comb = nx.path_weight(G2, path=s, weight="distanciacombinacion")
                    dictinterno = {"Nodo Origen": po,
                               "Nodo Destino": pde,
                               "Ruta": s,
                               "Distancia Total": l,
                               "DistanciaCaminando": c}
                    Aparallenar += [dictinterno]
        except nx.NetworkXNoPath:
            print("Pase por aca para el nodo " + str(po) + " y " + str(pde) + " y no encontre camino")
            continue

llenar = pd.DataFrame(Aparallenar)
llenar["Orden"] = llenar.groupby(["Nodo Origen", "Nodo Destino"]).cumcount()
llenar2 = llenar.groupby(["Nodo Origen", "Nodo Destino"]).size().to_frame("NroCaminos").reset_index()
llenar3 = llenar.groupby(["Nodo Origen", "Nodo Destino"])["DistanciaCaminando"].sum().reset_index()
llenar3 = llenar3.rename(columns={"DistanciaCaminando": "SumaCaminos"})
llenar3 = pd.merge(llenar2, llenar3, on=["Nodo Origen", "Nodo Destino"], how="inner")
llenar3 = pd.merge(llenar3, llenar, on=["Nodo Origen", "Nodo Destino"], how="inner")
llenar3["Numerador"] = (llenar3["SumaCaminos"] - llenar3["DistanciaCaminando"]) * (
            llenar3["NroCaminos"] - llenar3["Orden"] + 1)
llenar4 = llenar3.groupby(["Nodo Origen", "Nodo Destino"])["Numerador"].sum().reset_index()
llenar4 = llenar4.rename(columns={"Numerador": "SumaPesos"})
llenar4 = pd.merge(llenar3, llenar4, on=["Nodo Origen", "Nodo Destino"], how="left")
llenar4["p"] = llenar4["Numerador"] / llenar4["SumaPesos"]
llenar4["p"].fillna(1, inplace=True)
Rutas = llenar4[["Nodo Origen", "Nodo Destino", "Ruta", "Distancia Total", "DistanciaCaminando", "NroCaminos", "p"]]


def nrolineas(row):
    listinterna = row["Ruta"].copy()  # .copy() para que no modifique la lista original
    if len(listinterna) > 3:
        listinterna.pop(0)  # borro el nodo origen
        listinterna.pop()  # borro el nodo destino
        remove_digits = str.maketrans('', '', digits)  # funcion para eliminar los numeros
        listinterna = [i.translate(remove_digits) for i in listinterna]  # aplicacion de la funcion con lista comprimida
        listf = [s.replace('-', '') for s in listinterna]  # elimino los guiones
        aux = set(listf)  # elimino los duplicados
        nrolineas = len(aux)  # me fijo cuantos elementos quedaron en el set
        return nrolineas
    else:
        return 1


def asignarlinea(nl, s):
    listinterna = s.copy()  # .copy() para que no modifique la lista original
    if nl == 1 and len(s) >= 3:  # si hay una sola linea y hay mas de 3 nodos
        listinterna.pop(0)
        listinterna.pop()
        remove_digits = str.maketrans('', '', digits)
        listinterna = [i.translate(remove_digits) for i in listinterna]
        listf = [s.replace('-', '') for s in listinterna]
        return listf[0]
    elif nl == 2:
        listinterna.pop(0)
        listinterna.pop()
        remove_digits = str.maketrans('', '', digits)
        listinterna = [i.translate(remove_digits) for i in listinterna]
        listf = [s.replace('-', '') for s in listinterna]
        aux = set(listf)
        listf = list(aux)
        listf.sort()
        return listf[0] + "-" + listf[1]
    elif nl == 3:
        listinterna.pop(0)
        listinterna.pop()
        remove_digits = str.maketrans('', '', digits)
        listinterna = [i.translate(remove_digits) for i in listinterna]
        listf = [s.replace('-', '') for s in listinterna]
        aux = set(listf)
        listf = list(aux)
        listf.sort()
        return listf[0] + "-" + listf[1] + "-" + listf[2]
    else:
        return "Ninguna"


Rutas['NroLineas'] = Rutas.apply(lambda row: nrolineas(row), axis=1)  # asigno el nro de lineas utilizadas a cada linea
Rutas['LineaUtilizada'] = Rutas.apply(lambda x: asignarlinea(x['NroLineas'], x['Ruta']),
                                      axis=1)  # asigno las lineas utilizadas a cada linea
Rutas["Tupla"] = Rutas["Nodo Origen"] + "-" + Rutas["Nodo Destino"]
tuplasunicasabarcadas = Rutas["Tupla"].unique()
RutasCombinaciones = Rutas[Rutas["NroLineas"] == 2]  # separo las combinaciones
Rutas = Rutas[Rutas["NroLineas"] == 1]  # separo los viajes simples
combinaciones = RutasCombinaciones.groupby("LineaUtilizada").size().to_frame("Cuenta").reset_index().sort_values(
    by="Cuenta", ascending=False)

clustersagrupados.rename(columns={"h3_code_ida": "Nodo Origen", "h3_code_vuelta": "Nodo Destino"}, inplace=True)
clustersagrupados["Tupla"] = clustersagrupados["Nodo Origen"] + "-" + clustersagrupados["Nodo Destino"]

# Metricas simples de cobertura
tuplasunicasoriginales = clustersagrupados["Tupla"].unique()
clusternocubiertos = clustersagrupados[~clustersagrupados["Tupla"].isin(tuplasunicasabarcadas)]
trxsincubrir = clusternocubiertos.Cuenta.sum()
porcentaje_od_cubierto = ((len(od) - trxsincubrir) / len(od)) * 100
coberturaentuplas = 100 * len(tuplasunicasabarcadas) / len(tuplasunicasoriginales)

# Correccion factores
clusternocubiertos[["Nada", "LineaUtilizada"]] = clusternocubiertos["Nodo Origen"].str.split('-', -1, expand=True)
perdidas = clusternocubiertos.groupby("LineaUtilizada").size().to_frame("CargaLinea")
cargaOD = od.groupby("idlinea").size().to_frame("CargaLinea")
factorperdida = 1 - (perdidas / cargaOD)
factorperdida.rename(columns={"CargaLinea":"Factor"},inplace=True)
factores = pd.read_excel(direccionfactorecorreccion)
factores = pd.merge(factores, factorperdida, left_on="idlinea", right_on="LineaUtilizada", how="left")
factores["Factor"] = factores["Factor"] * factores["factor"]
factores.drop(columns=["factor"], inplace=True)
factores.to_excel(direccionfactorecorreccion, index=False)

# Analisis de carga
Carga = pd.merge(Rutas, clustersagrupados, on=["Nodo Origen", "Nodo Destino"], how="outer")
Carga["CargaProporcional"] = Carga["Cuenta"] * Carga["p"]
cargasimple = Carga.CargaProporcional.sum()

CargaCombinaciones = pd.merge(RutasCombinaciones, clustersagrupados, on=["Nodo Origen", "Nodo Destino"], how="outer")
CargaCombinaciones[['Linea1', 'Linea2']] = CargaCombinaciones["LineaUtilizada"].str.split('-', -1, expand=True)
CargaCombinaciones["CargaProporcional"] = CargaCombinaciones["Cuenta"] * CargaCombinaciones["p"]
CargaCombinaciones["CargaProporcional"].fillna(0, inplace=True)
CargaCombinaciones = CargaCombinaciones[CargaCombinaciones["CargaProporcional"] != 0]
TotalCombinaciones = CargaCombinaciones.CargaProporcional.sum()

cargacombaux1 = CargaCombinaciones.groupby("Linea1")["CargaProporcional"].sum().reset_index()
cargacombaux2 = CargaCombinaciones.groupby("Linea2")["CargaProporcional"].sum().reset_index()
test1 = (cargacombaux1.CargaProporcional.sum() + cargacombaux2.CargaProporcional.sum()) == 2*TotalCombinaciones
print("Test 1:" + str(test1))
cargacombaux2.rename(columns={"Linea2": "Linea1"}, inplace=True)
cargacombtotal = pd.merge(cargacombaux1, cargacombaux2, on="Linea1", how="outer")
cargacombtotal.fillna(0, inplace=True)
cargacombtotal["TotalCombinacion"] = cargacombtotal["CargaProporcional_x"] + cargacombtotal["CargaProporcional_y"]
combinacionestotales = cargacombtotal.TotalCombinacion.sum()

CargaValidacion = Carga.groupby("LineaUtilizada")["CargaProporcional"].sum().reset_index()
CargaValidacion = CargaValidacion[CargaValidacion["LineaUtilizada"].isin(Lineas)]
cargacombtotal.rename(columns={"Linea1": "LineaUtilizada"}, inplace=True)
CargaValidacion = pd.merge(CargaValidacion, cargacombtotal, on="LineaUtilizada", how="left")
CargaValidacion.fillna(0, inplace=True)
CargaValidacion["Total"] = CargaValidacion["TotalCombinacion"] + CargaValidacion["CargaProporcional"]
CargaValidacion = CargaValidacion[["LineaUtilizada", "Total"]]
Total = CargaValidacion.Total.sum()
CargaValidacion["%RelativoGTRX"] = CargaValidacion["Total"] * 100 / Total
if Validacion:
    CargaOriginal = od.groupby("idlinea").size().to_frame("CargaAbsOriginal").reset_index()
    CargaOriginal["%RelativoOriginal"] = CargaOriginal["CargaAbsOriginal"] * 100 / CargaOriginal[
        "CargaAbsOriginal"].sum()
    CargaValidacion = pd.merge(CargaValidacion, CargaOriginal, left_on="LineaUtilizada", right_on="idlinea", how="left")
CargaValidacion.to_excel(directoriosalida + "CargaValidacion.xlsx", index=False)

# Matriz de distribucion de lineas
# por paradas
Rutas[["Nada", "LineaOriginal"]] = Rutas["Nodo Origen"].str.split("-", 1, expand=True)
Matriz = Rutas.groupby(["LineaOriginal", "LineaUtilizada"]).size().to_frame("Cuenta").reset_index()
Matriz = Matriz.pivot(index="LineaOriginal", columns="LineaUtilizada", values="Cuenta")


# por carga
Carga[["Nada", "LineaOriginal"]] = Carga["Nodo Origen"].str.split("-", 1, expand=True)
Carga = Carga[Carga["LineaUtilizada"] != "Ninguna"]
MatrizCarga = Carga.groupby(["LineaOriginal", "LineaUtilizada"])["CargaProporcional"].sum().to_frame(
    "Cuenta").reset_index()
MatrizCarga = MatrizCarga.pivot(index="LineaOriginal", columns="LineaUtilizada", values="Cuenta")
MatrizCarga["Total"] = MatrizCarga.sum(axis=1)
MatrizCarga = MatrizCarga.div(MatrizCarga["Total"], axis=0)
MatrizCarga = MatrizCarga.drop(columns=["Total"])
MatrizCarga = MatrizCarga.mul(100)
if Validacion:
    diagonal = np.diag(MatrizCarga)
    lineasindices = list(MatrizCarga.index)
    for i, d in enumerate(diagonal):
        if d < 75:
            print("Linea: " + lineasindices[i] + " tiene menos de 75% de carga en su propia linea")
        else:
            pass


MatrizCarga.to_excel(directoriosalida + "MatrizCarga.xlsx")

# cargaporhora
puntospuros = od.copy(deep=True)
puntospuros["FECHATRX"] = pd.to_datetime(puntospuros["FECHATRX"])
puntospuros["HORA"] = puntospuros["FECHATRX"].dt.hour
paradasporhora = puntospuros.groupby(["HORA", "h3_code_ida", "h3_code_vuelta"]).size().to_frame("Cuenta").reset_index()
paradasporhora.rename(columns={"h3_code_ida": "Nodo Origen", "h3_code_vuelta": "Nodo Destino"}, inplace=True)


#################################################################################################
# Cargahora = pd.merge(Rutas, paradasporhora, on=["Nodo Origen", "Nodo Destino"], how="outer")
# CargaCombinacionesHora = pd.merge(RutasCombinaciones, paradasporhora, on=["Nodo Origen", "Nodo Destino"], how="outer")
#
# Cargahora["CargaProporcional"] = Cargahora["Cuenta"] * Cargahora["p"]
# CargaCombinacionesHora[['Linea1', 'Linea2']] = CargaCombinacionesHora["LineaUtilizada"].str.split('-', -1, expand=True)
# CargaCombinacionesHora["CargaProporcional"] = CargaCombinacionesHora["Cuenta"] * CargaCombinacionesHora["p"]
# CargaCombinacionesHora["CargaProporcional"].fillna(0, inplace=True)
# CargaCombinacionesHora = CargaCombinacionesHora[CargaCombinacionesHora["CargaProporcional"] != 0]
#
# cargacombaux1hora = CargaCombinacionesHora.groupby(["Linea1", "HORA"])["CargaProporcional"].sum().reset_index()
# cargacombaux2hora = CargaCombinacionesHora.groupby(["Linea2", "HORA"])["CargaProporcional"].sum().reset_index()
# cargacombaux2hora.rename(columns={"Linea2": "Linea1"}, inplace=True)
# cargacombtotalhora = pd.merge(cargacombaux1hora, cargacombaux2hora, on=["Linea1", "HORA"], how="outer")
# cargacombtotalhora.fillna(0, inplace=True)
# cargacombtotalhora["TotalCombinacion"] = cargacombtotalhora["CargaProporcional_x"] + cargacombtotalhora[
#     "CargaProporcional_y"]
#
# CargaValidacionHora = Cargahora.groupby(["HORA", "LineaUtilizada"])["CargaProporcional"].sum().reset_index()
# cargacombtotalhora.rename(columns={"Linea1": "LineaUtilizada"}, inplace=True)
# CargaValidacionHora = pd.merge(CargaValidacionHora, cargacombtotalhora, on=["LineaUtilizada", "HORA"], how="outer")
# CargaValidacionHora["TotalCombinacion"].fillna(0, inplace=True)
# CargaValidacionHora["Total"] = CargaValidacionHora["TotalCombinacion"] + CargaValidacionHora["CargaProporcional"]
# CargaValidacionHora = CargaValidacionHora[["LineaUtilizada", "Total", "HORA"]]
#
# Totalhora = CargaValidacionHora.Total.sum()
#
# CargaValidacionHora["%RelativoGTRX"] = CargaValidacionHora["Total"] * 100 / Totalhora
##############################################################################################
# carga por tramo
filtroaux = Carga[Carga["LineaUtilizada"] != "Ninguna"]
filtroaux.dropna(subset=["Ruta"], inplace=True)
filtroaux = filtroaux[["Nodo Origen", "Nodo Destino", "CargaProporcional", "Ruta"]]
filtroauxcomb = CargaCombinaciones[CargaCombinaciones["LineaUtilizada"] != "Ninguna"]
filtroauxcomb.dropna(subset=["Ruta"], inplace=True)
filtroauxcomb = filtroauxcomb[["Nodo Origen", "Nodo Destino", "CargaProporcional", "Ruta"]]
filtroaux = pd.concat([filtroaux, filtroauxcomb])
filtroaux["Ruta"] = filtroaux["Ruta"].astype(str)
filtroaux['Ruta'] = filtroaux['Ruta'].apply(ast.literal_eval)
rutasseparadas = filtroaux["Ruta"].tolist()
rutasseparadas = [sublist[1:-1] for sublist in rutasseparadas]
cargasparatramos = filtroaux["CargaProporcional"].tolist()
nodos = list(set([nodo for ruta in rutasseparadas for nodo in ruta]))
resultados = {p: 0 for p in nodos}
for i, ruta in enumerate(rutasseparadas):
    for p, r in enumerate(ruta):
        if p + 1 < len(ruta):
            resultados[r] += cargasparatramos[i]
        else:
            continue
# Filtrar los resultados para obtener solo aquellos con un peso distinto de cero
resultados_filtrados = {k: v for k, v in resultados.items() if v != 0}
# Convertir los resultados filtrados en un DataFrame
df_resultados = pd.DataFrame.from_dict(resultados_filtrados, orient='index', columns=['Peso']).reset_index()
df_resultados[['LineaUtilizada', 'Orden']] = df_resultados["index"].str.split('-', -1, expand=True)
df_resultados.to_excel(directoriosalida + "CargaPorPM.xlsx", index=False)

####Salidas
clustersagrupadost.rename(columns={"h3_code_ida": "Nodo Origen", "h3_code_vuelta": "Nodo Destino"}, inplace=True)
clustersagrupadost["DIA"] = clustersagrupadost["FECHATRX"].dt.day
clustersagrupadost["MES"] = clustersagrupadost["FECHATRX"].dt.month
paradasagrupadast_mean = clustersagrupadost.groupby(
    ["Nodo Origen", "Nodo Destino", "DIA", "FECHATRX", "MES"]).mean().reset_index()
paradasagrupadast_mean["HORA"] = paradasagrupadast_mean["FECHATRX"].dt.hour
paradasagrupadast_mean["MINUTO"] = paradasagrupadast_mean["FECHATRX"].dt.minute
ndias = len(paradasagrupadast_mean.groupby(["DIA", "MES"]).size().reset_index())

paradasagrupadastpromedio = paradasagrupadast_mean.groupby(["Nodo Origen", "Nodo Destino", "HORA", "MINUTO"])[
    ["Cuenta"]].sum().reset_index()
paradasagrupadastpromedio["Cuenta"] = paradasagrupadastpromedio["Cuenta"] / ndias

RutasCompletas = pd.concat([Rutas, RutasCombinaciones])
RutasCompletas.to_csv(directoriosalida + "RutasCompletas.csv", index=False)
paradasagrupadastpromedio.to_csv(directoriosalida + "ParadasPromedio.csv", index=False)


#Carga por tramo por hora
filtroaux = Carga[Carga["LineaUtilizada"] != "Ninguna"]
filtroaux.dropna(subset=["Ruta"], inplace=True)
filtroaux = filtroaux[["Nodo Origen", "Nodo Destino", "Ruta", "p", "LineaUtilizada"]]
filtroauxcomb = CargaCombinaciones[CargaCombinaciones["LineaUtilizada"] != "Ninguna"]
filtroauxcomb.dropna(subset=["Ruta"], inplace=True)
filtroauxcomb = filtroauxcomb[["Nodo Origen", "Nodo Destino", "Ruta", "p", "LineaUtilizada"]]
filtroaux = pd.concat([filtroaux, filtroauxcomb])
filtroaux["Ruta"] = filtroaux["Ruta"].astype(str)
filtroaux['Ruta'] = filtroaux['Ruta'].apply(ast.literal_eval)

filtroauxhora = pd.merge(filtroaux, paradasagrupadastpromedio, on = ["Nodo Origen", "Nodo Destino"], how = "outer")
filtroauxhora.dropna(subset=["Ruta"], inplace=True)
filtroauxhora["CargaProporcional"] = filtroauxhora["p"] * filtroauxhora["Cuenta"]
listavacia = []
for h in filtroauxhora.HORA.unique():
    auxint = filtroauxhora[filtroauxhora["HORA"] == h]
    rutasint = auxint["Ruta"].tolist()
    rutasint = [sublist[1:-1] for sublist in rutasint]
    cargasparatramos = auxint["CargaProporcional"].tolist()
    nodosint = list(set([nodo for ruta in rutasint for nodo in ruta]))
    dictint1 = {p: 0 for p in nodosint}
    for i, ruta in enumerate(rutasint):
        for p, r in enumerate(ruta):
            if p + 1 < len(ruta):
                dictint1[r] += cargasparatramos[i]
            else:
                continue
    dictint2 = {k: v for k, v in dictint1.items() if v != 0}
    pmporhora = pd.DataFrame.from_dict(dictint2, orient='index', columns=['Peso']).reset_index()
    pmporhora["HORA"] = h
    pmporhora.set_index(["index", "HORA"], inplace=True)
    pmporhora = pmporhora.T
    dictint2 = pmporhora.to_dict()
    listavacia.append(dictint2)

pmporhora = pd.DataFrame(listavacia).T.stack().reset_index()
pmporhora = pmporhora[["level_0", 0]]
pmporhora["Nodo"] = pmporhora["level_0"].apply(lambda x: x[0])
pmporhora["HORA"] = pmporhora["level_0"].apply(lambda x: x[1])
pmporhora["Personas"] = pmporhora[0].apply(lambda x: x["Peso"])
pmporhora = pmporhora[["Nodo", "HORA", "Personas"]]

pmporhora.to_csv(directoriosalida + "CargaPMPorHora.csv", index=False)




# Mapa de combinaciones

RutasCombinaciones = RutasCompletas[(RutasCompletas["NroLineas"] == 2) & (RutasCompletas["p"] == 1)].reset_index()
RutasCombinaciones = pd.merge(RutasCombinaciones, shpclusters, left_on="Nodo Origen", right_on="h3_code", how="left")
RutasCombinaciones = pd.merge(RutasCombinaciones, shpclusters, left_on="Nodo Destino", right_on="h3_code", how="left")
#
mapa2 = folium.Map(coordenadasmapa, tiles="OpenStreetMap", zoom_start=14)

for i, (no, nd) in enumerate(zip(RutasCombinaciones["Nodo Origen"], RutasCombinaciones["Nodo Destino"])):
    longnd1 = RutasCombinaciones.at[i, "geometry_x"].x
    latnd1 = RutasCombinaciones.at[i, "geometry_x"].y
    longnd2 = RutasCombinaciones.at[i, "geometry_y"].x
    latnd2 = RutasCombinaciones.at[i, "geometry_y"].y
    folium.PolyLine(locations=[[latnd1, longnd1], [latnd2, longnd2]], color="red", weight=0.2, opacity=0.2).add_to(mapa2)
loc = "Combinaciones"
title_html = '''
                    <h3 align="center" style="font-size:16px"><b>{}</b></h3>
                   '''.format(loc)
mapa2.get_root().html.add_child(folium.Element(title_html))
mapa2.save(directoriosalida + "MapaCombinaciones.html")  # Guardo.


# Combinaciones Tipicas
RutasCombinaciones = RutasCompletas[(RutasCompletas["NroLineas"] == 2) & (RutasCompletas["p"] == 1)].reset_index()
RutasCombinaciones = pd.merge(RutasCombinaciones, shpclusters, left_on="Nodo Origen", right_on="h3_code", how="left")
RutasCombinaciones = pd.merge(RutasCombinaciones, shpclusters, left_on="Nodo Destino", right_on="h3_code", how="left")
parescombinaciones = RutasCombinaciones.groupby(["LineaUtilizada"]).size().to_frame("Cuenta").reset_index().sort_values(
    by="Cuenta", ascending=False)
parescombinaciones.to_excel(directoriosalida + "ParesCombinaciones.xlsx", index=False)


mapa3 = folium.Map(coordenadasmapa, tiles="OpenStreetMap", zoom_start=14)

for l in Lineas:
    fg = folium.FeatureGroup("LINEA" + " " + l, show=False)
    dfinterno = df_resultados[df_resultados["LineaUtilizada"] == l]
    dfinterno["Orden"] = dfinterno["Orden"].astype(int)
    dfinterno = dfinterno.sort_values(by="Orden", ascending=True)
    nodosnp = np.array(dfinterno["index"])
    pesosnp = np.array(dfinterno["Peso"])
    maximo = np.amax(pesosnp) * 0.05
    dfinterno["PesoInverso"] = dfinterno['Peso'].rank(ascending=False)
    pesoscolores = np.array(dfinterno["PesoInverso"])
    maximo_colores = dfinterno["PesoInverso"].max()
    minimo_colores = dfinterno["PesoInverso"].min()
    mapa_colores = generate_color_scale(minimo_colores, maximo_colores)
    peso = []
    for i, nd in enumerate(nodosnp):
        try:
            longnd1 = G.nodes[nd]["geometry"].x
            latnd1 = G.nodes[nd]["geometry"].y
            longnd2 = G.nodes[nodosnp[i + 1]]["geometry"].x
            latnd2 = G.nodes[nodosnp[i + 1]]["geometry"].y
            folium.PolyLine([[latnd1, longnd1], [latnd2, longnd2]], color=mapa_colores(pesoscolores[i]), weight=pesosnp[i] / maximo,
                            popup=pesosnp[i]).add_to(fg)
        except IndexError:
            continue
    fg.add_to(mapa3)
folium.LayerControl().add_to(mapa3)
loc = "Mapa de carga por tramos"
title_html = '''
                    <h3 align="center" style="font-size:16px"><b>{}</b></h3>
                   '''.format(loc)
mapa3.get_root().html.add_child(folium.Element(title_html))
mapa3.save(directoriosalida + "MapaCargaPorTramos.html")  # Guardo.


# # Mapa paradas no cubiertas
mapa4 = folium.Map(coordenadasmapa, tiles="OpenStreetMap", zoom_start=14)
clusternocubiertos[['latitude_Origen', 'longitude_Origen']] = clusternocubiertos['geometry_x'].apply(lambda x: pd.Series({'latitude': x.centroid.y, 'longitude': x.centroid.x}))
clusternocubiertos[['latitude_Destino', 'longitude_Destino']] = clusternocubiertos['geometry_y'].apply(lambda x: pd.Series({'latitude': x.centroid.y, 'longitude': x.centroid.x}))
clusternocubiertos['h3_code_ida'] = clusternocubiertos.apply(lambda row: h3.geo_to_h3(row['latitude_Origen'], row['longitude_Origen'], 7), axis=1)
clusternocubiertos['h3_code_vuelta'] = clusternocubiertos.apply(lambda row: h3.geo_to_h3(row['latitude_Destino'], row['longitude_Destino'], 7), axis=1)
clusternocubiertos = clusternocubiertos[['h3_code_ida', 'h3_code_vuelta', "LineaUtilizada", "Cuenta"]]
clusternocubiertos = clusternocubiertos.groupby(["h3_code_ida", "h3_code_vuelta", "LineaUtilizada"])["Cuenta"].sum().reset_index()
clusternocubiertos["geometry_x"] = clusternocubiertos["h3_code_ida"].apply(lambda x: h3pol(x))
clusternocubiertos["geometry_y"] = clusternocubiertos["h3_code_vuelta"].apply(lambda x: h3pol(x))


for l in LineasViejas:
    fg = folium.FeatureGroup("LINEA" + " " + str(l), show=False)
    dfinterno = clusternocubiertos[clusternocubiertos["LineaUtilizada"] == str(l)].reset_index()
    dfinterno = dfinterno[dfinterno["Cuenta"] > len(od)*0.0025]
    dfinterno.reset_index(inplace=True)
    maximo = np.amax(dfinterno["Cuenta"]) *0.25
    for i, (gi, gv) in enumerate(zip(dfinterno["geometry_x"], dfinterno["geometry_y"])):
        longnd1 = gi.centroid.x
        latnd1 = gi.centroid.y
        longnd2 = gv.centroid.x
        latnd2 = gv.centroid.y
        folium.PolyLine([[latnd1, longnd1], [latnd2, longnd2]], color="#12bc8e", weight=dfinterno.at[i, "Cuenta"]/maximo,
                        popup=dfinterno.at[i, "Cuenta"]).add_to(fg)

    fg.add_to(mapa4)
folium.LayerControl().add_to(mapa4)
loc = "Cluster no conectados"
title_html = '''
                    <h3 align="center" style="font-size:16px"><b>{}</b></h3>
                   '''.format(loc)
mapa4.get_root().html.add_child(folium.Element(title_html))
mapa4.save(directoriosalida + "ClustersNoCubiertos.html")

end = time.time()
tiempototal = (end - start) / 60
print("El tiempo transcurrido fue de: " + str(tiempototal))

# metricas

metricas = pd.DataFrame({
    '%OD': [porcentaje_od_cubierto],
    'TRX sin cubrir': [trxsincubrir],
    '%Tuplas Abarcadas': [coberturaentuplas], "Tiempo de ejecucion": [tiempototal],
    "%Combinaciones": [(combinacionestotales / (combinacionestotales + cargasimple)) * 100], "NroCombinaciones": [combinacionestotales/2],
}).T

metricas.to_excel(directoriosalida + "Metricas.xlsx", index=True, sheet_name="ResumenSimulador")
