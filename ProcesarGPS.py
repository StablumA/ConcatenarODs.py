import pandas as pd
import numpy as np
import shapely
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiLineString, MultiPoint, LineString
import folium

# gps = pd.read_csv("Preprocessdata\gpstotal.csv", sep= ",")
# gpsTotal = gpsTotal[["codigoentidad", "idlinea", "interno", "date_time", "longitude", "latitude", "c_control_point", "file_id"]]
# gpsTotal = gpsTotal[gpsTotal['idlinea'].notna()]
# gpsTotal ["idlinea"]=gpsTotal ["idlinea"].astype(int)
# #
# listaoriginal= [763,758,756,759,764,754,751,757,762,750,760,748,755,867,2704]
# lista = [9,18,13,8,15,10,5,16,3,1,14,4,11,2,21]
# gpsTotal["idlinea"] =gpsTotal["idlinea"].replace(listaoriginal, lista).astype(int)
#
# gpsTotal .to_csv("Preprocessdata/gpslimpio.csv")
#
# gps = pd.read_csv("Preprocessdata\gpslimpio.csv", sep= ",",  parse_dates=['date_time'], dayfirst=True)
# gps["MES"] = gps["date_time"].dt.month
# gps["DIA"] = gps["date_time"].dt.day
#
# gpsJunio = gps[gps["MES"] == 6]
#
#
# gpsJunio.to_csv("Preprocessdata/gpsJunio.csv")
# gpsJunio["HORA"] = gpsJunio["date_time"].dt.hour
# ###aux1 con una linea
gpsJunio = pd.read_csv("Preprocessdata\gpsJunio.csv", sep= ",",  parse_dates=['date_time'], dayfirst=True)
gpsJunio.drop(["Unnamed: 0.1", "Unnamed: 0"], axis= 1, inplace=True)
gpsJunio.drop_duplicates(keep= "first", inplace=True)
gpsJunio["HORA"] = gpsJunio["date_time"].dt.hour
lineas = [9]


def calcularvueltas(df, fileid):
    dfint = df[df["file_id"] == fileid].reset_index()  # filtro al chofer
    long = -60.71766136238571
    lat = -31.616246578779528
    d = {"id": ["1"], "geometry": [Point(long, lat)]}
    gdf = gpd.GeoDataFrame(d, crs=4326)
    gdf = gdf.to_crs(epsg=32663)
    vueltas = 0
    M = len(dfint)
    for p in range(M):
        punto = {"id": ["1"], "geometry": [Point(dfint.at[p,"longitude"], dfint.at[p, "latitude"])]}
        gdf2 = gpd.GeoDataFrame(punto, crs=4326)
        gdf2 = gdf2.to_crs(epsg=32663)
        distancia = gdf.distance(gdf2)[0]
        if p != M-1:
            tiempo = (dfint.at[p+1, "date_time"] - dfint.at[p, "date_time"]).seconds // 60
        else:
            tiempo = 0
        if distancia <= 1000 and tiempo > 50:
            vueltas = vueltas + 1
        else:
            continue
    return vueltas


for l in lineas:
    gpsJuniointerno= gpsJunio[gpsJunio["idlinea"] == l]
    #si filtro en un dia y busco los valores unicos de file id puedo sacar cuantas vueltas de esa linea hubo en el dia

    dias = gpsJuniointerno["DIA"].unique()
    dias = np.sort(dias)
    L = len(dias)
    dfinterno = pd.DataFrame(columns= ["Dia","NroVueltas", "NroVehiculospromedioporhora", "DistanciaRecorrida"])

    B= gpd.read_file("SHP Files/Recorridos/Recorrido_"+ str(l) + "ilineas.shp")
    B = B.to_crs(epsg=32663)
    C = gpd.read_file("SHP Files/Recorridos/Recorrido_"+ str(l) + "vlineas.shp")
    C = C.to_crs(epsg=32663)
    dist1 = B.at[0, "geometry"].distance(B.at[len(B)-1, "geometry"])
    dist2 = C.at[0, "geometry"].distance(C.at[len(C)-1, "geometry"])
    distotal = dist1+ dist2
    for i in range(L):
        #por cada dia del mes
        auxiliar = gpsJuniointerno[gpsJuniointerno["DIA"] == dias[i]] #filtro el dia
        Choferes = auxiliar["file_id"].unique()
        contador = 0
        for w in Choferes:
            contador = contador + math.ceil((calcularvueltas(auxiliar,w) / 2))
        NroVueltas = contador
        aux2 = auxiliar.groupby(["HORA", "file_id"]).size().to_frame('Vehiculosporhora').reset_index()
        aux3 = aux2.groupby(["HORA"]).size().to_frame('Vehiculosporhora').reset_index() #Me da la cantidad de vehiculos que circulan en cada hora
        NroVehiculospromedioporhora = aux3["Vehiculosporhora"].mean()
        DistanciaRecorrida = distotal * NroVueltas
        dict = {"Dia": dias[i], "NroVueltas" : NroVueltas,"NroVehiculospromedioporhora": NroVehiculospromedioporhora, "DistanciaRecorrida" : DistanciaRecorrida}
        dfinterno = dfinterno.append(dict, ignore_index=True)

    dfinterno["NroVueltas"] =dfinterno["NroVueltas"].astype(int)
    dfinterno["DistanciaRecorrida"] =dfinterno["DistanciaRecorrida"].astype(int)

# dfinterno.to_excel("Output/AnalisisGPS/linea" + str(l)+ ".xlsx")


#Prueba para ver donde terminan los coles
mapit = folium.Map( location= [-31.633150571408255, -60.714748612677475])
folium.GeoJson(data=SHrecorridos["geometry"]).add_to(mapit)
mapit.save("Output/AnalisisGPS/Pruebas.html")