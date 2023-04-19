import geopandas as gpd
import geojson
import folium
import pandas as pd
from folium.plugins import GroupedLayerControl
from folium.plugins import FeatureGroupSubGroup
lineas=['1A','1B','2A','2B','3A','3B','4A','4B']
A=len(lineas)
#cargo datos de atractores
educacion=gpd.read_file('./rawdata/Educacion_2021/Educacion_2021.shp').to_crs(4326)
administrativos=gpd.read_file('./rawdata/Equipamientos_adminitrativos.gpkg').to_crs(4326)
religiosos= gpd.read_file('./rawdata/Establecimientos_religiosos.gpkg').to_crs(4326)

# cargo datos de infraestructura
censales=gpd.read_file('./rawdata/Radios Censales/Radios censales.shp')
semaforos=gpd.read_file('./rawdata/Semaforos/Semaforos.shp').to_crs(4326)
ciclovia=gpd.read_file('./rawdata/CICLOVÍAS/CICLOVÍAS.dbf').to_crs(4326)

# vamos a leer y mergear nuestros datos de OD
origenes=gpd.read_file('./rawdata/EOD - 2022/Solo Origen.dbf').to_crs(4326)
destinos= gpd.read_file('./rawdata/EOD - 2022/Solo Destino.dbf').to_crs(4326)

# eliminar registros sin geometría
origenes = origenes[origenes.geometry.notnull()]
destinos = destinos[destinos.geometry.notnull()]

# necesito mergearlos en un nuevo archivo
origenes=gpd.GeoDataFrame(origenes)[['ID','LINEA','geometry']]
destinos= gpd.GeoDataFrame(destinos)[['ID','geometry']]
mergeJ= pd.merge(origenes,destinos,how='inner', on='ID')
mergeJ= mergeJ.dropna()
mergeJ=mergeJ.rename(columns={'geometry_x':'O','geometry_y':'D'})
df=mergeJ

# genero el mapa
m = folium.Map(tiles='OpenStreetMap', zoom_start=14,location=[-33.74648, -61.9679])
# leo los recorridos
recorridos= gpd.read_file('./rawdata/TUP/Recorrido_TUP.shp')
# creo lista de recorridos
B=[]
for w in range(A):
    reco = recorridos.iloc[w] # recorrido en tratamiento
    auxiliar = df
    auxiliar=auxiliar[auxiliar['LINEA']==lineas[w]]
        # filtrado por linea

    flg1= folium.FeatureGroup("Recorrido "+str(lineas[w]), show=False, overlay=True, control= True)
    folium.GeoJson(data=reco['geometry']).add_to(flg1)
    flg1.add_to(m)


    flg2=folium.FeatureGroup("Orígenes y Destinos "+str(lineas[w]), show=False, overlay=True, control= True)
    for i,(n1,n2) in enumerate(zip(auxiliar['O'],auxiliar['D'])): # tengo que separar las coordenadas para hacer el polyline
        x1= n1.x
        y1= n1.y
        x2= n2.x
        y2= n2.y
        folium.PolyLine([[y1,x1],[y2,x2]], color = "blue",weight=0.5,opacity=0.6 ).add_to(flg2)
    flg2.add_to(m)
    B.append(flg1)
    B.append(flg2)


#cargo la infraestrucuta en C
C=[]
fg_c= folium.FeatureGroup("Ciclovias",show=False,control=True, overlay=True)
fg_se=  folium.FeatureGroup("Semáforos",show=False,control=True, overlay=True)
fg_cs= folium.FeatureGroup("Radios Censales",show=False,control=True, overlay=True)
                        # cargo los datos de las capas
folium.GeoJson(data=censales).add_to(fg_cs)
folium.GeoJson(data=ciclovia).add_to(fg_c)
folium.GeoJson(data=semaforos).add_to(fg_se)
# modifico iconos de las que me interesen
for i,h in enumerate(semaforos['geometry']):
    folium.Marker(
        location=[h.y,h.x],
        icon= folium.Icon(color='darkgreen',icon="traffic-light",prefix='fa'),
        popup=semaforos['Tipo']
    ).add_to(fg_se)

                        # las agrego al mapa
fg_c.add_to(m)
fg_cs.add_to(m)
fg_se.add_to(m)
                        # las agrego a la lista de 'INFRAESTRUCTURA'
C.append(fg_c)
C.append(fg_cs)
C.append(fg_se)

#cargo los atractores  en D
D=[]

fg_ad= folium.FeatureGroup("Administrativos",show=False,control=True, overlay=True)
fg_e=  folium.FeatureGroup("Establecimientos educativos",show=False,control=True, overlay=True)
fg_r= folium.FeatureGroup("Establecimientos religiosos",show=False,control=True, overlay=True)

folium.GeoJson(data=religiosos).add_to(fg_r)
folium.GeoJson(data=educacion).add_to(fg_e)
folium.GeoJson(data=administrativos).add_to(fg_ad)

for i, h in enumerate(religiosos["geometry"]):
    folium.Marker(
        location=[h.y, h.x],
        icon=folium.Icon(color="green", icon="person-praying", prefix='fa'),
        popup=religiosos["fna"][i]
    ).add_to(fg_r)

for i, h in enumerate(educacion["geometry"]):
    folium.Marker(
        location=[h.y, h.x],
        icon=folium.Icon(color="black", icon="book-open", prefix='fa'),
        popup=educacion["Nombre"][i]
    ).add_to(fg_e)

for i, h in enumerate(administrativos["geometry"]):
    folium.Marker(
        location=[h.y, h.x],
        icon=folium.Icon(color="black", icon="building", prefix='fa'),
        popup=administrativos["Nombre"][i]
    ).add_to(fg_ad)

fg_r.add_to(m)
fg_ad.add_to(m)
fg_e.add_to(m)

D.append(fg_r)
D.append(fg_ad)
D.append(fg_e)

#LIN.add_to(m)
GroupedLayerControl(groups={'Líneas':B,'Infraestructura':C, 'Atractores':D},exclusive_groups=False,collapsed=True).add_to(m)
m.save("./Output/ODMap.html")
print(' mapa generado ')












#for w in range(A):
#    a = gpd.read_file("Preprocessdata/Recorridos/Recorrido_"+ str(lineas[w]) + ".shp")
#    df = pd.read_excel("Output/MODML/linea"+str(lineas[w])+".xlsx")
#    auxiliar = df
#
#    auxiliar['start_lon']=auxiliar['Origen'].str.extract('\((.+?),').astype(float)
#    auxiliar['end_lon']=auxiliar['Destino'].str.extract('\((.+?),').astype(float)
#    auxiliar['start_lat']=auxiliar['Origen'].str.extract(',(.+?)\)').astype(float)
#    auxiliar['end_lat']=auxiliar['Destino'].str.extract(',(.+?)\)').astype(float)
#
#    puntosVecinalIda = gpd.GeoDataFrame(auxiliar, geometry=gpd.points_from_xy(auxiliar['start_lon'], auxiliar['start_lat']), crs="4326")
#    puntosEnVecinalIda= gpd.sjoin_nearest(puntosVecinalIda, SHvecinales, how="inner")
#
#
#    puntosVecinalVuelta = gpd.GeoDataFrame(auxiliar, geometry=gpd.points_from_xy(auxiliar['end_lon'], auxiliar['end_lat']), crs="4326")
#    puntosEnVecinalVuelta= gpd.sjoin_nearest(puntosVecinalVuelta, SHvecinales, how="inner")
#
#    puntosEnIdayVuelta=puntosEnVecinalIda.merge(puntosEnVecinalVuelta, on = "NROTARJETA")
#    puntosEnIdayVuelta = puntosEnIdayVuelta.drop(["Origen_x", "Destino_x", "index_right_x", "Origen_y", "Destino_y", "NUMERO_x", "NUMERO_y", "index_right_y"], axis=1)
#
#    prueba1=puntosEnIdayVuelta.groupby(['VECINAL_x', 'VECINAL_y']).size().to_frame('Cuenta').reset_index()
#    prueba1.to_excel("Output/MODML/Linea" + str(lineas[w])+ "porvecinales.xlsx")
#
#    L = len(prueba1)
#    M = prueba1["Cuenta"].max()
#    puntos = []
#    peso = []
#    peso2 = []
#    for i in range(L):
#        indice1 = prueba1["VECINAL_x"].values[i]
#        indice2 = prueba1["VECINAL_y"].values[i]
#        latitudorigen= SHvecinales.loc[SHvecinales['VECINAL'] == indice1, "lat"].values[0]
#        longitudorigen= SHvecinales.loc[SHvecinales['VECINAL'] == indice1, "long"].values[0]
#        latituddestino = SHvecinales.loc[SHvecinales['VECINAL'] == indice2, "lat"].values[0]
#        longituddestino = SHvecinales.loc[SHvecinales['VECINAL'] == indice2, "long"].values[0]
#        puntos.append([[latitudorigen, longitudorigen],[latituddestino, longituddestino]])
#        peso.append(((prueba1["Cuenta"][i])/(0.35*M)).astype(str))
#        peso2.append(prueba1["Cuenta"][i])
#    fg = folium.FeatureGroup("LINEA" + " " + str(lineas[w]))
#    for m in range(L):
#        folium.PolyLine(puntos[m] , color = "red", weight = peso[m], popup = peso2[m]).add_to(fg)
#    folium.GeoJson(data=a["geometry"]).add_to(fg)
#    fg.add_to(mapa)
#folium.LayerControl().add_to(mapa)
