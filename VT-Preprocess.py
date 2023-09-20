import folium
import pandas as pd
import geopandas as gpd
from fiona.drvsupport import supported_drivers
import numpy as np
import warnings

supported_drivers['kml'] = 'rw' # enable KML support which is disabled by default
supported_drivers['KML'] = 'rw'
warnings.filterwarnings('ignore')
supported_drivers['LIBKML'] = 'rw'



mes = "Mayo"
#Preprocessdata de trx
df= pd.read_csv("Venado/Rawdata/05 Transacciones "+mes+" 2022.csv",sep=";", parse_dates=["Fecha Hora"], dayfirst=True)
df= df[df["Fecha Hora"]>= pd.to_datetime('1-05-2022', format='%d-%m-%Y')]

df["Fecha Hora"] = pd.to_datetime(df["Fecha Hora"])
df["MES"]= df["Fecha Hora"].dt.month
df["DIA"] = df["Fecha Hora"].dt.day
df= df[df["MES"]== 5]

## filtro los duplicados
bool_series= df.duplicated(keep='first')
df= df[~bool_series]

df["Ramal"] = df["Ramal"].str.split(" ",expand=True)[1]

df= df[["Ramal","Tipo Trx","Fecha Hora","Hora","Tarjeta","Latitud","Longitud","DIA","MES"]]
df= df.rename(columns={"Ramal":"idlinea","Fecha Hora":"FECHATRX", "Longitud":"latitude", "Latitud":"longitude", "Tarjeta" : "NROTARJETA", "Tipo Trx" : "CODIGOCONTRATO" })
Lineas = df["idlinea"].unique()
originales = ["1A","1B","2A","2B","3A","3B","4A","4B"]
nuevas=["AA","AB","BA","BB","CA","CB","DA","DB"]
df["idlinea"] = df["idlinea"].replace(originales, nuevas)

TRXperdidas = pd.DataFrame()
dias = df["DIA"].unique()
diasexcluidos = []
PxD = pd.DataFrame()
longitudes = [len(df[df["DIA"]==d]) for d in dias]
longitudprom = np.mean(longitudes)
for d in dias:
    if len(df[df["DIA"]==d]) < 0.05 * longitudprom:
        diasexcluidos.append(d)
        aux = pd.DataFrame({"Tipo perdida":["Día "+str(d)+" excluido"],"Transacciones perdidas":[len(df[df["DIA"]==d])]})
        PxD = pd.concat([PxD,aux])
TRXperdidas = pd.concat([TRXperdidas,PxD])
df = df[~df["DIA"].isin(diasexcluidos)]

dforiginal=df.copy(deep=True )# modificado


df.to_csv("Venado/Preprocessdata/VT_dfbasico_"+mes+".csv")



#Preprocessdata de shp
shp = gpd.read_file("Venado/SHP Files/TUP 2023/TUP - 2023.shp")
shp.rename(columns={"Línea":"Linea"}, inplace=True)
Lineas = shp["Linea"].unique()
#nuevaslineas = [2, 1, 3, 4, 6, 5, 8, 7, 50, 90]
#shp["Linea"] = shp["Linea"].replace(Lineas, nuevaslineas)

for l in Lineas:
    shpinterno = shp[shp["Linea"]== l]
    shpinterno.to_file("Venado/Preprocessdata/Recorridos/Recorrido_"+ str(l) + ".shp")


#Preprocessdata de radios censales
shprc = gpd.read_file("SHP Files/Radios Censales VT/RC.shp")
#shprc = shprc[["Frac", "Radio","Personas", "geometry"]]
shprc.rename(columns={"Frac":"frac", "Radio":"radio"}, inplace=True)
shprc["RC"] = shprc["frac"].astype(str) + "-" + shprc["radio"].astype(str)
shprc.to_file("SHP Files/Radios Censales VT/RC.shp")
# Versión 2 -> Corto los RC
shprc = gpd.read_file("SHP Files/Radios Censales VT/RCpob.shp")
shpres= gpd.read_file("SHP Files/Radios Censales VT/Zona censable con verdes.kml")
shpver= gpd.read_file("SHP Files/Radios Censales VT/Verdes.kml")
# Elimino la parte de la zona residencial que esté en la intersección con las zonas verdes
for index, row in shpver.iterrows():
    interseccion = shpres.intersection(row.geometry) # Busco interseccion
    shpres["geometry"] = shpres.difference(interseccion) # Elimino la interseccion
    # deberia actualizarse en cada corrida la geometria de shpres

shpF = gpd.GeoDataFrame(shprc)
for index,row in shprc.iterrows():
    geometriafinal = shpres.intersection(row.geometry)[0]
    shpF.at[index,"geometry"] = geometriafinal
shpF.to_file("SHP Files/Radios Censales VT/RC_VT_Cortado.shp")
gdfRC= gpd.GeoDataFrame(shpF,geometry="geometry")
destinomapaRC = "SHP Files/Radios Censales VT/mapaRCcortado.html"
mapaRC = folium.Map(location=[-33.7324, -61.9891], zoom_start=12)
gdfRC["geometry"] = gdfRC['geometry'].to_crs(epsg=3857)        # aparto a la geometría
folium.GeoJson(gdfRC).add_to(mapaRC)  # Grafico los polígonos
mapaRC.save(destinomapaRC)
print("Listo los RC cortados")

#Preprocessdata de paradas
paradas = pd.read_excel("Venado/Rawdata/Paradas TUP MTV.xlsx")
import re
import numpy as np

regex = r'Linea\s*\s*\s*'

paradas['Linea'] = paradas['Linea'].apply(lambda x: re.sub(regex, '', x))   # limpio 'Linea '
paradas['Linea'] = paradas['Linea'].apply(lambda x: x.strip())

condicion = [paradas["Linea"] == "1 A", paradas["Linea"] == "1 B", paradas["Linea"] == "2 A", paradas["Linea"] == "2 B", paradas["Linea"] == "3 A", paradas["Linea"] == "3 B", paradas["Linea"] == "4 A" , paradas["Linea"] == "4 B"]
valor = ["unoA" , "unoB", "dosA", "dosB", "tresA", "tresB", "cuatroA", "cuatroB"]
paradas["LineaBien"] = np.select(condicion, valor, paradas["Linea"])
paradas = paradas[["LineaBien", "Parada", "Latitud", "Longitud"]]
ids = range(1, len(paradas)+1)
ids = [str(x).zfill(4) for x in ids]

# Agregamos la columna "ID" al dataframe
paradas['ID'] = ids
paradas.rename(columns={"LineaBien":"Linea"}, inplace=True)
paradas.to_excel("Venado/Preprocessdata/paradasVT.xlsx")


shp = gpd.read_file("Venado/SHP Files/Recorridos/Recorrido_TUP.shp")
shp.rename(columns={"Línea":"Linea"}, inplace=True)

condicion = [shp["Linea"] == "1A", shp["Linea"] == "1B", shp["Linea"] == "2A", shp["Linea"] == "2B", shp["Linea"] == "3A", shp["Linea"] == "3B", shp["Linea"] == "4A" , shp["Linea"] == "4B"]
valor = ["unoA" , "unoB", "dosA", "dosB", "tresA", "tresB", "cuatroA", "cuatroB"]
shp["Linea"] = np.select(condicion, valor, shp["Linea"])

Lineas = shp["Linea"].unique()
nuevaslineas = [2, 1, 3, 4, 6, 5, 8, 7, 50, 90]
shp["Linea"] = shp["Linea"].replace(Lineas, nuevaslineas)

for l in shp.Linea.unique():
    shpinterno = shp[shp["Linea"]== l]
    shpinterno.to_file("Venado/Preprocessdata/Recorridos/Recorrido_"+ str(l) + ".shp")



#preprocess od

od1 = pd.read_csv("Venado/Output/ODfinalmayo.csv")
od2 = pd.read_csv("Venado/Output/ODfinalseptiembre.csv")
od3 = pd.read_csv("Venado/Output/ODfinalnoviembre.csv")

od = pd.concat([od1, od2, od3])
od = od[["idlinea", "NROTARJETA", "FECHATRX", "latitude_Origen", "longitude_Origen", "latitude_Destino", "longitude_Destino", "LineaCombinacion", "Tipificación"]]

conversionlineas = pd.read_excel("Venado/conversionLineas.xlsx")
originales = conversionlineas["original"].tolist()
nuevas = conversionlineas["nueva"].tolist()
od["idlinea"] = od["idlinea"].replace(nuevas, originales)

condicion = [od["idlinea"] == "1A", od["idlinea"] == "1B", od["idlinea"] == "2A", od["idlinea"] == "2B", od["idlinea"] == "3A", od["idlinea"] == "3B", od["idlinea"] == "4A" , od["idlinea"] == "4B"]
valor = ["unoA" , "unoB", "dosA", "dosB", "tresA", "tresB", "cuatroA", "cuatroB"]
od["idlinea"] = np.select(condicion, valor, od["idlinea"])

od.to_csv("Venado/Output/odfinal.csv")



#preprocess factores correccion
c1 = pd.read_csv("Venado/Output/estadisticosodporlineamayo.csv")
c2 = pd.read_csv("Venado/Output/estadisticosodporlineaseptiembre.csv")
c3 = pd.read_csv("Venado/Output/estadisticosodporlineanoviembre.csv")

c = pd.concat([c1, c2, c3])

conversionlineas = pd.read_excel("Venado/conversionLineas.xlsx")
originales = conversionlineas["original"].tolist()
nuevas = conversionlineas["nueva"].tolist()
c["idlinea"] = c["idlinea"].replace(nuevas, originales)

condicion = [c["idlinea"] == "1A", c["idlinea"] == "1B", c["idlinea"] == "2A", c["idlinea"] == "2B", c["idlinea"] == "3A", c["idlinea"] == "3B", c["idlinea"] == "4A" , c["idlinea"] == "4B"]
valor = ["unoA" , "unoB", "dosA", "dosB", "tresA", "tresB", "cuatroA", "cuatroB"]
c["idlinea"] = np.select(condicion, valor, c["idlinea"])

factorescoreccion = c.groupby("idlinea")[["% por linea"]].mean().reset_index()
factorescoreccion.rename(columns={"% por linea":"factor"}, inplace=True)
factorescoreccion.to_excel("Venado/Output/CorreccionCargas.xlsx")
