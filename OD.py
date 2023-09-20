import pandas as pd
import geopandas as gpd
import warnings
import numpy as np
import shapely
from shapely.ops import transform
from pyproj import Transformer
from shapely.geometry import LineString
import random
import time
from dividirlinea import tramos_basico
warnings.filterwarnings('ignore')




inicio= time.time()
#Parametros
cotainftiempo = 5 #tiempo maximo para que dos trx se consideren de dos personas distintas
cotasuptiempo = 40 #tiempo min para que dos trx se consideren independientes (no son una combinacion ni dos personas)
distanciatol = 1000 #distancia de tolerancia para caminata
lineasincluidas = ["AA","AB","BA","BB","CA","CB","DA","DB"]
lineasexcluidas = []
diasferiados = [25]
distanciaclustering = 100
clustering = True
perdidatolerancia = [] # Perdidas por tolerancia
direcciondf = "Venado/Preprocessdata/VT_dfbasico_Mayo.csv"
direccionshp = "Venado/Preprocessdata/Recorridos/Recorrido_" #recorridos
direccionzonificacion = "Venado/SHP Files/Radios Censales VT/RC_VT_Cortado.shp" #zonificacion predecir n=1
direcciondestinood = "Venado/Output/ODfinalMayo.csv"
direcciondestinot1 = "Venado/Output/estadisticosodporlineaMayo.xlsx"
direcciondestinoreporte = "Venado/Output/reporteODMayo.xlsx"
zona = "RC"
tipodedia = "Habil"
diaparapromediar = "Monday"
#Parte general a cualquier od nx
#Preprocesamiento
df = pd.read_csv(direcciondf, sep=",", parse_dates=['FECHATRX'], dayfirst=True) #Leo el archivo con los datos

# Corroboramos nombres de las columnas df
columnas = df.columns
Caux = [col for col in columnas if col not in ["NROTARJETA","FECHATRX","MES","longitude","latitude","idlinea","CODIGOCONTRATO"]]
if len(Caux) > 0:
    print('Error: DataFrame con columnas erróneas ' + str(Caux))
# Corroboramos que los idlinea no tengan números
for l in df["idlinea"].unique():
    for caracter in l:
        if caracter.isdigit():
            print('Error: El nombre de la línea ' + str(l) +' contiene el número ' + str(caracter))


# Corroboramos ingreso de MES como número y no texto
if df["MES"].dtype != "int64":
    try:
        df["MES"] = df["MES"].astype(int)
    except ValueError:
        print("Error: Columna MES no ingresada como número entero")

try:
    df['FECHATRX'] = pd.to_datetime(df['FECHATRX'])
except ValueError:
    print("Error: Columna 'FECHATRX' en formato incorrecto")

try:
    df = df.drop(["Unnamed: 0"], axis=1)
except KeyError:
     df = df

df["MINUTO"]=df['FECHATRX'].dt.minute
df['HORA']=df['FECHATRX'].dt.hour
df['DIA']=df['FECHATRX'].dt.day
df['nDIA']=df['FECHATRX'].dt.day_name() #le agrego los dias por nombre
df = df[~df["idlinea"].isin(lineasexcluidas)]
df = df[df["idlinea"].isin(lineasincluidas)]
df['xy']="("+df['longitude'].astype(str)+", "+df['latitude'].astype(str)+")" # Genera una columna (long, lat)
df['latitude'] = (df['latitude']).astype(str).str.replace(',', '.').astype(float) # modificado
df['longitude'] = df['longitude'].astype(str).str.replace(',', '.').astype(float) # modificado

df = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']), crs='EPSG:4326')
df = df.to_crs(epsg=22185)
#df = df.astype({'idlinea':'int'}) #Esto porque estaban como float  # no irÃ­a ahora que sean string
Lineas = df['idlinea'].unique()


if tipodedia == "Habil":
    dffiltrado = df[df['nDIA'].isin(['Wednesday', 'Thursday', 'Friday', 'Monday',
           'Tuesday'])]
    dffiltrado = dffiltrado[~dffiltrado["DIA"].isin(diasferiados)]
    dffiltrado = gpd.GeoDataFrame(dffiltrado, geometry=gpd.points_from_xy(dffiltrado['longitude'], dffiltrado['latitude']), crs="4326").to_crs(epsg=22185) #Armo un geodataframe con mi informacion
    dias = dffiltrado["DIA"].unique()
elif tipodedia == "Sabado":
    dffiltrado = df[df['nDIA'].isin(['Saturday'])]
    dffiltrado = dffiltrado[~dffiltrado["DIA"].isin(diasferiados)]
    dffiltrado = gpd.GeoDataFrame(dffiltrado,
                                  geometry=gpd.points_from_xy(dffiltrado['longitude'], dffiltrado['latitude']),
                                  crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
    dias = dffiltrado["DIA"].unique()
else:
    dffiltrado = df[df['nDIA'].isin(['Sunday'])]
    dffiltrado = dffiltrado[~dffiltrado["DIA"].isin(diasferiados)]
    dffiltradoaux = df[df['nDIA'].isin(['Wednesday', 'Thursday', 'Friday', 'Monday',
                                     'Tuesday'])]
    dffiltradoaux = dffiltradoaux[dffiltradoaux["DIA"].isin(diasferiados)]
    dffiltrado = pd.concat([dffiltrado, dffiltradoaux])
    dffiltrado = gpd.GeoDataFrame(dffiltrado,
                                  geometry=gpd.points_from_xy(dffiltrado['longitude'], dffiltrado['latitude']),
                                  crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
    dias = dffiltrado["DIA"].unique()
dffiltradoOriginal = dffiltrado.copy(deep=True) # modif
dffiltrado = dffiltrado[dffiltrado["latitude"] != 0] # modif
dffiltrado = dffiltrado[dffiltrado["longitude"] != 0] # modif

perdidagps = len(dffiltrado) - len(dffiltradoOriginal) # modif

df = []
#funciones generales
def casosn(dataframe, n):
    """Funcion para filtrar por numero de transacciones por dia"""
    dfaux1 = dataframe['NROTARJETA'].value_counts().reset_index()
    dfaux1.columns = ['ID', "Cantidad"]  # Renombra las columnas de dfc, creada en la linea anterior
    dfaux2=dfaux1['Cantidad'].value_counts().reset_index() # Crea una nueva auxiliar, para ver la distribucion de cantidad de viajes.
    dfaux2.columns=['Viajes en dia', "Cantidad"] #Renombra las columnas
    dfaux1=dfaux1[dfaux1['Cantidad']==n] # Filtra viajes igual a 2.
    salida=dataframe[dataframe['NROTARJETA'].isin(dfaux1['ID'])] #Crea una nueva dataframe, filtrando la dataframe llamada auxiliar para aquellos NROTARJETA contenidos en dfc.
    salida=salida.sort_values(by=['NROTARJETA','FECHATRX']) #Ordena los valores
    salida = salida.reset_index()
    return salida

def minutototal(registro):
    "Calcular el minuto total de una transaccion"
    auxiliar = 60*int(registro["HORA"]) + int(registro["MINUTO"])
    return auxiliar



def ExtraerMODdl(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-1):
        prop1 = dataframe.at[index+1,"idlinea"] != dataframe.at[index,"idlinea"] #Distinta linea
        prop2 = (dataframe.at[index+1,"minutototal"] - dataframe.at[index,"minutototal"]) >= cotasuptiempo #Al menos 40 minutos de diferencia entre las transacciones
        prop3 = dataframe.at[index+1,"NROTARJETA"] == dataframe.at[index,"NROTARJETA"] #Misma tarjeta
        if prop1 and prop2 and prop3:
            dictaux = {"NROTARJETA": dataframe.at[index,"NROTARJETA"], "FECHATRX": dataframe.at[index,"FECHATRX"], "longitude": dataframe.at[index,"longitude"], "latitude": dataframe.at[index,"latitude"], "idlinea": dataframe.at[index,"idlinea"], "CODIGOCONTRATO": dataframe.at[index,"CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index+1,"NROTARJETA"], "FECHATRX": dataframe.at[index+1,"FECHATRX"], "longitude": dataframe.at[index+1,"longitude"], "latitude": dataframe.at[index+1,"latitude"], "idlinea": dataframe.at[index+1,"idlinea"], "CODIGOCONTRATO": dataframe.at[index+1,"CODIGOCONTRATO"]}
            copiainterna.append(dictaux)
            copiainterna.append(dictaux2)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]

    a1 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="first")
    a2 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="last")
    a3 = a1.merge(a2, on='NROTARJETA', how='inner', suffixes=("_Inicial", "_Final"))
    return a3

def ExportarMODdl(dataframe):                                                                               # modificado
    B = ExtraerMODdl(dataframe)
    #Distancia entre linea final y punto inicial
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_Final"] ==l] #filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp") #cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        a = a.set_geometry(
            a.geometry.map(
        lambda linestring: shapely.ops.transform(lambda x, y, *_: (x, y), linestring)
))
        shpinterno = gpd.GeoDataFrame(dfinterno, geometry=gpd.points_from_xy(dfinterno['longitude_Inicial'], dfinterno['latitude_Inicial']), crs="4326").to_crs(epsg=22185) #Armo un geodataframe con mi informacion
        axt = tramos_basico(a,100) # divido la linea en tramos
        axt = axt.set_geometry("inicio")
        axt= axt[["inicio"]]
        axt=axt.to_crs(epsg=22185)
        cruce= gpd.sjoin_nearest(shpinterno,axt,how='left', distance_col= "DistanciaconLineaFinal")
        cruce['geometry'] = [axt['inicio'][i] for i in cruce['index_right']]
        cruce["geometry"]=cruce["geometry"].to_crs(epsg=4326)
        cruce["longitude_Inicial"],cruce["latitude_Inicial"]=cruce["geometry"].x,cruce["geometry"].y  # reemplazo por los nuevos puntos
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_Final", "FECHATRX_Inicial"],keep="first", inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["geometry","index_right"], axis=1, inplace=True)


    #Distancia entre linea inicial y punto Final
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Inicial"] ==l] #filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        a = a.set_geometry(
            a.geometry.map(
                lambda linestring: shapely.ops.transform(lambda x, y, *_: (x, y), linestring)
))
        shpinterno = gpd.GeoDataFrame(dfinterno, geometry=gpd.points_from_xy(dfinterno['longitude_Final'], dfinterno['latitude_Final']), crs="4326").to_crs(epsg=22185) #Armo un geodataframe con mi informacion
        axt = tramos_basico(a,100) # divido la linea en tramos
        axt =  axt.set_geometry("inicio")
        axt= axt[["inicio"]]
        axt = axt.to_crs(epsg=22185)
        cruce= gpd.sjoin_nearest(shpinterno,axt,how='left', distance_col= "DistanciaconLineaInicial")
        cruce['geometry'] = [axt['inicio'][i] for i in cruce['index_right']]
        cruce["geometry"]=cruce["geometry"].to_crs(epsg=4326)
        cruce["longitude_Final"],cruce["latitude_Final"]=cruce["geometry"].x,cruce["geometry"].y  # reemplazo por los nuevos puntos
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_Final", "FECHATRX_Inicial"],keep="first", inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["geometry","index_right"], axis=1, inplace=True)

    Bini = len(F)
    F = F[F["DistanciaconLineaInicial"] <= distanciatol]
    F = F[F["DistanciaconLineaFinal"] <= distanciatol]
    Bfin = len(F) # perdida x tolerancia
    Bfuera=Bini-Bfin
    perdidatolerancia.append(Bfuera)


    dfx = F[F.filter(regex = ("._Inicial|latitude_Final|longitude_Final|NROTARJETA")).columns]
    dfx.rename(columns = {'FECHATRX_Inicial':'FECHATRX','idlinea_Inicial':'idlinea', 'latitude_Inicial':'latitude_Origen', 'longitude_Inicial':'longitude_Origen', 'latitude_Final':'latitude_Destino', 'longitude_Final':'longitude_Destino', 'CODIGOCONTRATO_Inicial':'CODIGOCONTRATO'}, inplace = True) #cambio el nombre por uno mas feliz

    dfy = F[F.filter(regex = ("._Final|latitude_Inicial|longitude_Inicial|NROTARJETA")).columns]
    dfy.rename(columns = {'FECHATRX_Final':'FECHATRX','idlinea_Final':'idlinea', 'latitude_Final':'latitude_Origen', 'longitude_Final':'longitude_Origen', 'latitude_Inicial':'latitude_Destino', 'longitude_Inicial':'longitude_Destino', 'CODIGOCONTRATO_Final':'CODIGOCONTRATO'}, inplace = True) #cambio el nombre por uno mas feliz

    A1 = [dfy,dfx]
    A2 = pd.concat(A1)
    A2["Tipificación"] = "N2Dl"
    return A2

def ExtraerMODml(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-1):
        prop1 = dataframe.at[index+1,"idlinea"] == dataframe.at[index,"idlinea"] #misma linea
        prop2 = (dataframe.at[index+1,"minutototal"] - dataframe.at[index,"minutototal"]) >= cotasuptiempo #Al menos 40 minutos de diferencia entre las transacciones
        prop3 = dataframe.at[index+1,"NROTARJETA"] == dataframe.at[index,"NROTARJETA"] #Misma tarjeta
        if prop1 and prop2 and prop3:
            dictaux = {"NROTARJETA": dataframe.at[index,"NROTARJETA"], "FECHATRX": dataframe.at[index,"FECHATRX"], "longitude": dataframe.at[index,"longitude"], "latitude": dataframe.at[index,"latitude"], "idlinea": dataframe.at[index,"idlinea"], "CODIGOCONTRATO": dataframe.at[index,"CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index+1,"NROTARJETA"], "FECHATRX": dataframe.at[index+1,"FECHATRX"], "longitude": dataframe.at[index+1,"longitude"], "latitude": dataframe.at[index+1,"latitude"], "idlinea": dataframe.at[index+1,"idlinea"], "CODIGOCONTRATO": dataframe.at[index+1,"CODIGOCONTRATO"]}
            copiainterna.append(dictaux)
            copiainterna.append(dictaux2)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    a1 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="first")
    a2 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="last")
    a3 = a1.merge(a2, on='NROTARJETA', how='inner', suffixes=("_Inicial", "_Final"))
    return a3


def ExportarMODml(dataframe):
    B = ExtraerMODml(dataframe)
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_Final"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno, geometry=gpd.points_from_xy(dfinterno['longitude_Inicial'],
                                                                             dfinterno['latitude_Inicial']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaconLineaFinal")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_Final", "FECHATRX_Inicial"], keep="first", inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea inicial y punto Final
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Inicial"] == l]  # filtro la linea inicial
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno, geometry=gpd.points_from_xy(dfinterno['longitude_Final'],
                                                                             dfinterno['latitude_Final']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaconLineaInicial")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_Final", "FECHATRX_Inicial"], keep="first", inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    transformer = Transformer.from_crs("epsg:4326", "epsg:22185")


    def get_distance(k):
        lstr = LineString([(k.longitude_Inicial, k.latitude_Inicial), (k.longitude_Final, k.latitude_Final)])
        line2 = transform(transformer.transform, lstr)
        return line2.length

    F["Distancia"] = F.apply(get_distance, axis=1)
    Bini = len(F)
    G = F[F["Distancia"] >=500]
    dfx = G[G.filter(regex=("._Inicial|latitude_Final|longitude_Final|NROTARJETA")).columns]
    dfx.rename(
        columns={'FECHATRX_Inicial': 'FECHATRX', 'idlinea_Inicial': 'idlinea', 'latitude_Inicial': 'latitude_Origen',
                 'longitude_Inicial': 'longitude_Origen', 'latitude_Final': 'latitude_Destino',
                 'longitude_Final': 'longitude_Destino', 'CODIGOCONTRATO_Inicial': 'CODIGOCONTRATO'},
        inplace=True)  # cambio el nombre por uno mas feliz

    dfy = G[G.filter(regex=("._Final|latitude_Inicial|longitude_Inicial|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Final': 'FECHATRX', 'idlinea_Final': 'idlinea', 'latitude_Final': 'latitude_Origen',
                        'longitude_Final': 'longitude_Origen', 'latitude_Inicial': 'latitude_Destino',
                        'longitude_Inicial': 'longitude_Destino', 'CODIGOCONTRATO_Final': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz
    Bfin = len(G)
    #filtro = ~ Bini.isin(Bfin)
    #Bfuera = Bini[filtro].dropna()
    #Bfuera["Tipificación"] = "Ptol" # perdida x tolerancia
    Bfuera= Bini-Bfin
    perdidatolerancia.append(Bfuera)

    A1 = [dfx, dfy]
    A2 = pd.concat(A1)
    A2["Tipificación"] = "N2ML"
    return A2

def ExtraerMODn3(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-2):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) >= cotasuptiempo)  # Las primeras 2 trx tienen que darse en mas de 40 minutos
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index+1, "minutototal"]) >= cotasuptiempo)
        prop3 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[index + 1, "NROTARJETA"])  # Tiene que tratarse de la misma tarjeta
        prop4 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[index + 2, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        if prop4 and prop1 and prop3 and prop2:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX" : dataframe.at[index, "FECHATRX"], "longitude" : dataframe.at[index, "longitude"], "latitude" : dataframe.at[index, "latitude"], "idlinea" : dataframe.at[index, "idlinea"], "CODIGOCONTRATO" : dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"], "FECHATRX" : dataframe.at[index + 1, "FECHATRX"], "longitude" : dataframe.at[index + 1, "longitude"], "latitude" : dataframe.at[index + 1, "latitude"], "idlinea" : dataframe.at[index + 1, "idlinea"], "CODIGOCONTRATO" : dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"], "FECHATRX" : dataframe.at[index + 2, "FECHATRX"], "longitude" : dataframe.at[index + 2, "longitude"], "latitude" : dataframe.at[index + 2, "latitude"], "idlinea" : dataframe.at[index + 2, "idlinea"], "CODIGOCONTRATO" : dataframe.at[index + 2, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)

        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X','variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X", 'variable6': 'CODIGOCONTRATO_Y', 'variable7': 'FECHATRX_Y',
                                "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y", 'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z"})
    return df1

def ExportarMODn3(dataframe):
    B = ExtraerMODn3(dataframe).reset_index()
    # Distancia entre linea inicial y punto medio. Coloquialmente: quiero saber si el punto medio corresponde a un destino inicial valido
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_X"] == l]  # filtro la linea media
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaconLineaInicial")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea media y punto final
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Y"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaconLineaMedia")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    G = []
    # me dirijo al origen ?
    for l in Lineas:
        dfinterno = F[F["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaconLineaFinal")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    I = []
    # me dirijo al punto medio
    for l in Lineas:
        dfinterno = H[H["idlinea_Z"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaconLineaFinalconPuntoMedio")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        I.append(cruce)

    J = pd.concat(I)
    J.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    conditions = [
        (J["DistanciaconLineaInicial"] <= distanciatol) & (J["DistanciaconLineaMedia"] <= distanciatol) & (
                    J["DistanciaconLineaFinal"] <= distanciatol),
        (J["DistanciaconLineaInicial"] <= distanciatol) & (J["DistanciaconLineaMedia"] <= distanciatol) & (
                    J["DistanciaconLineaFinal"] > distanciatol) & (J["DistanciaconLineaFinalconPuntoMedio"] <= distanciatol),
        (J["DistanciaconLineaInicial"] <= distanciatol) & (J["DistanciaconLineaMedia"] <= distanciatol) & (
                    J["DistanciaconLineaFinal"] > distanciatol) & (J["DistanciaconLineaFinalconPuntoMedio"] > distanciatol),
        (J["DistanciaconLineaInicial"] > distanciatol) & (J["DistanciaconLineaMedia"] <= distanciatol) & (
                    J["DistanciaconLineaFinal"] <= distanciatol),
        (J["DistanciaconLineaInicial"] <= distanciatol) & (J["DistanciaconLineaMedia"] > distanciatol) & (
                    J["DistanciaconLineaFinal"] <= distanciatol),
        (J["DistanciaconLineaInicial"] > distanciatol) & (J["DistanciaconLineaMedia"] <= distanciatol) & (
                    J["DistanciaconLineaFinal"] > distanciatol),
        (J["DistanciaconLineaInicial"] > distanciatol) & (J["DistanciaconLineaMedia"] > distanciatol) & (
                    J["DistanciaconLineaFinal"] > distanciatol),
        (J["DistanciaconLineaInicial"] > distanciatol) & (J["DistanciaconLineaMedia"] > distanciatol) & (
                    J["DistanciaconLineaFinal"] <= distanciatol),
        (J["DistanciaconLineaInicial"] <= distanciatol) & (J["DistanciaconLineaMedia"] > distanciatol) & (
                    J["DistanciaconLineaFinal"] > distanciatol),
    ]
    valores = ["Triangulo", "Cuchara", "Vibora", "Ene", "Ce", "Casa", "Nada", "Escudo", "Recta"]
    J["Tipificación"] = np.select(conditions,
                                  valores)  # calculo la distancia promedio que se va a recorrer usando esa linea

    Triangulos = J[J["Tipificación"] == "Triangulo"]

    dfx = Triangulos[Triangulos.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Y': 'latitude_Destino',
                        'longitude_Y': 'longitude_Destino', "CODIGOCONTRATO_X": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfy = Triangulos[Triangulos.filter(regex=("._Y|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', "CODIGOCONTRATO_Y": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfz = Triangulos[Triangulos.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', "CODIGOCONTRATO_Z": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz

    A1 = [dfx, dfy, dfz]
    A2triangulos = pd.concat(A1)
    A2triangulos["Tipificación"] = "Triangulo"

    Cucharas = J[J["Tipificación"] == "Cuchara"]

    dfx = Cucharas[Cucharas.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Y': 'latitude_Destino',
                        'longitude_Y': 'longitude_Destino', "CODIGOCONTRATO_X": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfy = Cucharas[Cucharas.filter(regex=("._Y|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', "CODIGOCONTRATO_Y": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfz = Cucharas[Cucharas.filter(regex=("._Z|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_Y': 'latitude_Destino',
                        'longitude_Y': 'longitude_Destino', "CODIGOCONTRATO_Z": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz

    A1 = [dfx, dfy, dfz]
    A2cucharas = pd.concat(A1)
    A2cucharas["Tipificación"] = "Cuchara"

    Enes = J[J["Tipificación"] == "Ene"]

    dfy = Enes[Enes.filter(regex=("._Y|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', "CODIGOCONTRATO_Y": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfz = Enes[Enes.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', "CODIGOCONTRATO_Z": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz

    A1 = [dfy, dfz]
    A2enes = pd.concat(A1)
    A2enes["Tipificación"] = "Ene"

    Ces = J[J["Tipificación"] == "Ce"]

    dfx = Ces[Ces.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Y': 'latitude_Destino',
                        'longitude_Y': 'longitude_Destino', "CODIGOCONTRATO_X": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfz = Ces[Ces.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', "CODIGOCONTRATO_Z": "CODIGOCONTRATO"},
               inplace=True)  # cambio el nombre por uno mas feliz

    A1 = [dfx, dfz]
    A2ces = pd.concat(A1)
    A2ces["Tipificación"] = "Ce"

    A2total = [A2ces,A2enes, A2cucharas, A2triangulos]
    return pd.concat(A2total)

def ExtraerMODcf(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-2):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) >= cotasuptiempo)
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) <= cotasuptiempo)
        prop3 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[index + 1, "NROTARJETA"])
        prop4 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[index + 2, "NROTARJETA"]
        prop5 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) >= cotainftiempo)
        prop6 = (dataframe.at[index + 2, "idlinea"]) != (dataframe.at[index + 1, "idlinea"])
        if prop4 and prop1 and prop3 and prop2 and prop5 and prop6:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"], "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"], "idlinea": dataframe.at[index, "idlinea"], "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"], "FECHATRX": dataframe.at[index + 1, "FECHATRX"], "longitude": dataframe.at[index + 1, "longitude"], "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"], "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"], "FECHATRX": dataframe.at[index + 2, "FECHATRX"], "longitude": dataframe.at[index + 2, "longitude"], "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"], "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X",'variable6': 'CODIGOCONTRATO_Y', 'variable7': 'FECHATRX_Y',
                                "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y",'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z"})
    return df1


def ExportarMODn3cf(dataframe):
    B = ExtraerMODcf(dataframe).reset_index()

    # distancia entre linea media y punto final
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_Y"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea final y punto inicial
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea inicial y punto medio
    G = []
    for l in Lineas:
        dfinterno = F[F["idlinea_X"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    Bini = len(H)
    H = H[H["DistanciaLXconY"] <= distanciatol]
    H = H[H["DistanciaLZconX"] <= distanciatol]    ## para que pueda volver con la Z?
    H = H[H["DistanciaLYconZ"] <= distanciatol]
    Bfin = len(H)
    Bfuera= Bini-Bfin
    perdidatolerancia.append(Bfuera)



    # si es una combinacion al final el origen es el punto y y el destino el punto z

    dfx = H[H.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen',
                        'latitude_Y': 'latitude_Destino', 'longitude_Y': 'longitude_Destino',
                        'CODIGOCONTRATO_X': 'CODIGOCONTRATO'}, inplace=True)
    dfx["Tipificación"] = "n3cfps"
    dfz = H[H.filter(regex=("._Y|latitude_X|longitude_X|NROTARJETA|idlinea_Z")).columns]
    dfz.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen',
                        'latitude_X': 'latitude_Destino', 'longitude_X': 'longitude_Destino',
                        'CODIGOCONTRATO_Y': 'CODIGOCONTRATO', "idlinea_Z" : "LineaCombinacion"}, inplace=True)
    dfz["Tipificación"] = "n3cf"
    out = [dfx, dfz]
    salida = pd.concat(out)
    return salida

def ExtraerMODci(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-2):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) <= cotasuptiempo) # que la 2da sea combinacion
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) >= cotasuptiempo) # la 3ra no sea combinacion
        prop3 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[index + 1, "NROTARJETA"]) # misma tarjeta
        prop4 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[index + 2, "NROTARJETA"] #misma tarjeta
        prop5 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) >= cotainftiempo) # no sea prestada
        prop6 = (dataframe.at[index + 1, "idlinea"]) != (dataframe.at[index, "idlinea"]) # distinta linea entre x e y
        if prop4 and prop1 and prop3 and prop2 and prop5 and prop6:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X",'variable6': 'CODIGOCONTRATO_Y', 'variable7': 'FECHATRX_Y',
                                "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y",'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z"})
    return df1


def ExportarMODn3ci(dataframe):
    B = ExtraerMODci(dataframe).reset_index()
    # distancia entre linea media y punto final
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_Y"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea final y punto inicial
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea inicial y punto medio
    G = []
    for l in Lineas:
        dfinterno = F[F["idlinea_X"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    Bini = len(H)
    H = H[H["DistanciaLXconY"] <= distanciatol]
    H = H[H["DistanciaLZconX"] <= distanciatol]
    H = H[H["DistanciaLYconZ"] <= distanciatol]
    Bfin = len( H)
    Bfuera= Bini-Bfin
    perdidatolerancia.append(Bfuera)


    # si es una combinacion al principio el origen es el punto inicial y el destino el punto z

    dfx = H[H.filter(regex=("._X|latitude_Z|longitude_Z|NROTARJETA|idlinea_Y")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen',
                        'latitude_Z': 'latitude_Destino', 'longitude_Z': 'longitude_Destino',
                        'CODIGOCONTRATO_X': 'CODIGOCONTRATO', "idlinea_Y" : "LineaCombinacion"}, inplace=True)
    dfx["Tipificación"] = "n3ci"
    dfz = H[H.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen',
                        'latitude_X': 'latitude_Destino', 'longitude_X': 'longitude_Destino',
                        'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'}, inplace=True)
    dfz["Tipificación"] = "n3cips"
    out = [dfx, dfz]
    salida = pd.concat(out)
    return salida


def ExtraerMODdpf(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-2):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) >= cotasuptiempo) # Y no es comb
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) <= cotainftiempo) # que vayan 2 personas juntas
        prop3 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[index + 1, "NROTARJETA"])  # con la misma tarjeta
        prop4 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[index + 2, "NROTARJETA"]
        if prop4 and prop1 and prop3 and prop2:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X",'variable6': 'CODIGOCONTRATO_Y', 'variable7': 'FECHATRX_Y',
                                "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y",'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z"})
    return df1


def ExportarMODdpf(dataframe):
    B = ExtraerMODdpf(dataframe).reset_index()
    # distancia entre linea media y punto final
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_Y"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea final y punto inicial
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea inicial y punto medio
    G = []
    for l in Lineas:
        dfinterno = F[F["idlinea_X"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    Bini = len(H)
    H = H[H["DistanciaLXconY"] <= distanciatol]
    H = H[H["DistanciaLZconX"] <= distanciatol]
    H = H[H["DistanciaLYconZ"] <= distanciatol]
    Bfin = len(H)
    Bfuera= Bini-Bfin
    perdidatolerancia.append(Bfuera)




    # si es una combinacion al principio el origen es el punto inicial y el destino el punto z

    dfx = H[H.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen',
                        'latitude_Y': 'latitude_Destino', 'longitude_Y': 'longitude_Destino',
                        'CODIGOCONTRATO_X': 'CODIGOCONTRATO'}, inplace=True)

    dfz = H[H.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen',
                        'latitude_X': 'latitude_Destino', 'longitude_X': 'longitude_Destino',
                        'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'}, inplace=True)
    out = [dfx, dfz]
    salida = pd.concat(out)
    salida["Tipificación"] = "n3dpf"
    return salida


def ExtraerMODdpi(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-2):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) <= cotainftiempo)  #  comparten tarjeta
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) >= cotasuptiempo) # no es combinacion
        prop3 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[index + 1, "NROTARJETA"]) # misma tarjeta
        prop4 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[index + 2, "NROTARJETA"] #misma tarjeta
        if prop4 and prop1 and prop3 and prop2:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X",'variable6': 'CODIGOCONTRATO_Y', 'variable7': 'FECHATRX_Y',
                                "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y",'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z"})
    return df1


def ExportarMODn3pi(dataframe):
    B = ExtraerMODdpi(dataframe).reset_index()
    # distancia entre linea media y punto final
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_Y"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea final y punto inicial
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Distancia entre linea inicial y punto medio
    G = []
    for l in Lineas:
        dfinterno = F[F["idlinea_X"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    Bini= len(H)
    H = H[H["DistanciaLXconY"] <= distanciatol]
    H = H[H["DistanciaLZconX"] <= distanciatol]
    H = H[H["DistanciaLYconZ"] <= distanciatol]
    Bfin = len(H)
    Bfuera= Bini-Bfin
    perdidatolerancia.append(Bfuera)



    # si es una combinacion al principio el origen es el punto inicial y el destino el punto z

    dfx = H[H.filter(regex=("._X|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen',
                        'latitude_Z': 'latitude_Destino', 'longitude_Z': 'longitude_Destino',
                        'CODIGOCONTRATO_X': 'CODIGOCONTRATO'}, inplace=True)

    dfz = H[H.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen',
                        'latitude_X': 'latitude_Destino', 'longitude_X': 'longitude_Destino',
                        'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'}, inplace=True)
    out = [dfx, dfz]
    salida = pd.concat(out)
    salida["Tipificación"] = "n3dpi"
    return salida

def ExtraerMODn4s(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-3):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) >= cotasuptiempo) # no es combinacion
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) >= cotasuptiempo)  # no es combinacion
        prop3 = ((dataframe.at[index + 3, "minutototal"] - dataframe.at[index + 2, "minutototal"]) >= cotasuptiempo) # no es combinacion
        prop4 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[
            index + 1, "NROTARJETA"])  # Tiene que tratarse de la misma tarjeta
        prop5 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[
            index + 2, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        prop6 = dataframe.at[index + 2, "NROTARJETA"] == dataframe.at[
            index + 3, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        if prop4 and prop1 and prop3 and prop2 and prop6 and prop5:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            dictaux4= {"NROTARJETA": dataframe.at[index + 3, "NROTARJETA"], "FECHATRX": dataframe.at[index + 3, "FECHATRX"], "longitude": dataframe.at[index + 3, "longitude"], "latitude": dataframe.at[index + 3, "latitude"], "idlinea": dataframe.at[index + 3, "idlinea"], "CODIGOCONTRATO": dataframe.at[index + 3, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
            copiainterna.append(dictaux4)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X",'variable6': 'CODIGOCONTRATO_Y', 'variable7': 'FECHATRX_Y',
                                "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y",'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z", 'variable16': 'CODIGOCONTRATO_W', 'variable17': 'FECHATRX_W', "variable18": "idlinea_W",
                              "variable19": "latitude_W", "variable20": "longitude_W"})
    return df1


def ExportarMODn4s(dataframe):
    B = ExtraerMODn4s(dataframe).reset_index()
    # Distancia entre linea inicial y segundo punto. Coloquialmente: quiero saber si el segundo punto corresponde a un destino inicial valido
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_X"] == l]  # filtro la linea media
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z", "FECHATRX_W"],
                              keep="first", inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Me volvi al inicio?
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Y"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    G = []
    # me dirijo al cuarto punto?
    for l in Lineas:
        dfinterno = F[F["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_W'], dfinterno['latitude_W']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconW")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    I = []
    # me vuelvo?
    for l in Lineas:
        dfinterno = H[H["idlinea_W"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLWconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        I.append(cruce)

    J = pd.concat(I)
    J.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    K = []
    # y a z
    for l in Lineas:
        dfinterno = J[J["idlinea_Y"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        K.append(cruce)

    M = pd.concat(K)
    M.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    N = []
    # w a x
    for l in Lineas:
        dfinterno = M[M["idlinea_W"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLWconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        N.append(cruce)

    O = pd.concat(N)
    O.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    transformer = Transformer.from_crs("epsg:4326", "epsg:22185")

    def get_distancexz(k):
        lstr = LineString([(k.longitude_X, k.latitude_X,), (k.longitude_Z, k.latitude_Z)])
        line2 = transform(transformer.transform, lstr)
        return line2.length

    O["DistanciaXZ"] = O.apply(get_distancexz, axis=1)

    def get_distanceyw(k):
        lstr = LineString([(k.longitude_Y, k.latitude_Y,), (k.longitude_W, k.latitude_W)])
        line2 = transform(transformer.transform, lstr)
        return line2.length

    O["DistanciaYW"] = O.apply(get_distanceyw, axis=1)

    conditions = [
        (O["DistanciaXZ"] <= 500) & (O["DistanciaYW"] <= 500) & (O["idlinea_X"] == O["idlinea_Z"]) & (
                    O["idlinea_Y"] == O["idlinea_W"]),
        ((O["idlinea_X"] == O["idlinea_Y"]) | (O["DistanciaLYconX"] <= distanciatol)) & (O["DistanciaLWconX"] <= distanciatol) & (
                    O["DistanciaLZconW"] <= distanciatol)
    ]
    valores = ["DosVecesElMismoViaje", "Dos viajes distintos"]
    O["Tipificación"] = np.select(conditions, valores)

    DosVeces = O[O["Tipificación"] == "DosVecesElMismoViaje"]

    dfx = DosVeces[DosVeces.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Y': 'latitude_Destino',
                        'longitude_Y': 'longitude_Destino', 'CODIGOCONTRATO_X': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfy = DosVeces[DosVeces.filter(regex=("._Y|latitude_X|longitude_X|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_Y': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfz = DosVeces[DosVeces.filter(regex=("._Z|latitude_W|longitude_W|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_W': 'latitude_Destino',
                        'longitude_W': 'longitude_Destino', 'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfw = DosVeces[DosVeces.filter(regex=("._W|latitude_X|longitude_X|NROTARJETA")).columns]
    dfw.rename(columns={'FECHATRX_W': 'FECHATRX', 'idlinea_W': 'idlinea', 'latitude_W': 'latitude_Origen',
                        'longitude_W': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_W': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz

    A1 = [dfx, dfy, dfz, dfw]
    A2dosveces = pd.concat(A1)
    A2dosveces["Tipificación"] = "DosVeces"


    DosDistintos = O[O["Tipificación"] == "Dos viajes distintos"]

    dfx = DosDistintos[DosDistintos.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Y': 'latitude_Destino',
                        'longitude_Y': 'longitude_Destino', 'CODIGOCONTRATO_X': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfy = DosDistintos[DosDistintos.filter(regex=("._Y|latitude_X|longitude_X|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_Y': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfz = DosDistintos[DosDistintos.filter(regex=("._Z|latitude_W|longitude_W|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_W': 'latitude_Destino',
                        'longitude_W': 'longitude_Destino', 'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz
    dfw = DosDistintos[DosDistintos.filter(regex=("._W|latitude_X|longitude_X|NROTARJETA")).columns]
    dfw.rename(columns={'FECHATRX_W': 'FECHATRX', 'idlinea_W': 'idlinea', 'latitude_W': 'latitude_Origen',
                        'longitude_W': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_W': 'CODIGOCONTRATO'},
               inplace=True)  # cambio el nombre por uno mas feliz

    A1 = [dfx, dfy, dfz, dfw]
    A2dosdistintos = pd.concat(A1)   # me arma un df con todos los origenes - destinos, por linea. especificando codigo contrato, nro tarjeta, tipificacion
    A2dosdistintos["Tipificación"] = "DosDistintos"

    Otros = O[O["Tipificación"] == "0"]
    return pd.concat([A2dosdistintos,A2dosveces])

def ExtraerMOD4c(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-3):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) <= cotasuptiempo)   # ser combinacion Y
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) >= cotasuptiempo) # no ser combinacion
        prop3 = ((dataframe.at[index + 3, "minutototal"] - dataframe.at[index + 2, "minutototal"]) <= cotasuptiempo) # ser combinacion tambien
        prop4 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[
            index + 1, "NROTARJETA"])  # Tiene que tratarse de la misma tarjeta
        prop5 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[
            index + 2, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        prop6 = dataframe.at[index + 2, "NROTARJETA"] == dataframe.at[
            index + 3, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        prop7 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) > 5)  # no sean dos personas
        prop8 = ((dataframe.at[index + 3, "minutototal"] - dataframe.at[index + 2, "minutototal"]) > 5)  # no sean dos personas
        if prop4 and prop1 and prop3 and prop2 and prop6 and prop5 and prop8 and prop7:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            dictaux4 = {"NROTARJETA": dataframe.at[index + 3, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 3, "FECHATRX"],
                        "longitude": dataframe.at[index + 3, "longitude"],
                        "latitude": dataframe.at[index + 3, "latitude"], "idlinea": dataframe.at[index + 3, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 3, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
            copiainterna.append(dictaux4)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X",'variable6': 'CODIGOCONTRATO_Y', 'variable7': 'FECHATRX_Y',
                                "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y",'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z", 'variable16': 'CODIGOCONTRATO_W', 'variable17': 'FECHATRX_W', "variable18": "idlinea_W",
                              "variable19": "latitude_W", "variable20": "longitude_W"})
    return df1


def ExportarMODn4c(dataframe):
    B = ExtraerMOD4c(dataframe).reset_index()

    # Distancia entre linea inicial y segundo punto. Coloquialmente: quiero saber si el segundo punto corresponde a un destino inicial valido
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_X"] == l]  # filtro la linea media
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z", "FECHATRX_W"],
                              keep="first", inplace=True)
        C.append(cruce)            # coloquialmente: quiero ver si donde me puedo bajar de la X estÃ¡ cerca de donde me tomo la Y

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Y"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    G = []
    # me dirijo al cuarto punto?
    for l in Lineas:
        dfinterno = F[F["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_W'], dfinterno['latitude_W']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconW")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    I = []
    for l in Lineas:
        dfinterno = H[H["idlinea_W"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLWconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        I.append(cruce)
    J = pd.concat(I)
    J.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    conditions = [
        (J["DistanciaLXconY"] <= distanciatol) & (J["DistanciaLYconZ"] <= distanciatol) & (J["DistanciaLZconW"] <= distanciatol) & (
                    J["DistanciaLWconX"] <= distanciatol)]
    valores = ["n4combinacion"]
    J["Tipificación"] = np.select(conditions, valores)

    n4combinacion = J[J["Tipificación"] == "n4combinacion"]

    dfx = n4combinacion[n4combinacion.filter(regex=("._X|latitude_Z|longitude_Z|NROTARJETA|idlinea_Y")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', 'CODIGOCONTRATO_X': 'CODIGOCONTRATO', "idlinea_Y" : "LineaCombinacion"}, inplace=True)

    dfz = n4combinacion[n4combinacion.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA|idlinea_W")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_Z': 'CODIGOCONTRATO', "idlinea_W" : "LineaCombinacion"}, inplace=True)

    out = [dfx, dfz]
    A = pd.concat(out)
    A["Tipificación"] = "n4combinacion"
    return A

def ExtraerMOD4dp(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-3):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) <= cotainftiempo)   # se presta tarjeta
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index+1, "minutototal"]) >= cotasuptiempo) # no es comb ni prestada
        prop3 = ((dataframe.at[index + 3, "minutototal"] - dataframe.at[index + 2, "minutototal"]) <= cotainftiempo) # vuelvo a prestar tarjeta
        prop4 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[index + 1, "NROTARJETA"])  # Tiene que tratarse de la misma tarjeta
        prop5 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[index + 2, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        prop6 = dataframe.at[index + 2, "NROTARJETA"] == dataframe.at[index + 3, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        if prop4 and prop1 and prop3 and prop2 and prop6 and prop5:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            dictaux4 = {"NROTARJETA": dataframe.at[index + 3, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 3, "FECHATRX"],
                        "longitude": dataframe.at[index + 3, "longitude"],
                        "latitude": dataframe.at[index + 3, "latitude"], "idlinea": dataframe.at[index + 3, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 3, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
            copiainterna.append(dictaux4)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X", 'variable6': 'CODIGOCONTRATO_Y',
                              'variable7': 'FECHATRX_Y',
                              "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y", 'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z", 'variable16': 'CODIGOCONTRATO_W',
                              'variable17': 'FECHATRX_W', "variable18": "idlinea_W",
                              "variable19": "latitude_W", "variable20": "longitude_W"})
    return df1


def ExportarMOD4dp(dataframe):
    B = ExtraerMOD4dp(dataframe).reset_index()

    # Distancia entre linea inicial y segundo punto. Coloquialmente: quiero saber si el segundo punto corresponde a un destino inicial valido
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_X"] == l]  # filtro la linea media
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z", "FECHATRX_W"],
                              keep="first", inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Y"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    G = []
    # me dirijo al cuarto punto?
    for l in Lineas:
        dfinterno = F[F["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_W'], dfinterno['latitude_W']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconW")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    I = []
    for l in Lineas:
        dfinterno = H[H["idlinea_W"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLWconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        I.append(cruce)

    J = pd.concat(I)
    J.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    conditions = [
        (J["DistanciaLYconZ"] <= distanciatol) & (J["DistanciaLWconX"] <= distanciatol)]  # si Z es el destino de Y y X el destino de W
    valores = ["n4dp"]
    J["Tipificación"] = np.select(conditions, valores)

    n4dp = J[J["Tipificación"] == "n4dp"]

    dfx = n4dp[n4dp.filter(regex=("._X|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', 'CODIGOCONTRATO_X': 'CODIGOCONTRATO'}, inplace=True)

    dfy = n4dp[n4dp.filter(regex=("._Y|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', 'CODIGOCONTRATO_Y': 'CODIGOCONTRATO'}, inplace=True)

    dfz = n4dp[n4dp.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'}, inplace=True)

    dfw = n4dp[n4dp.filter(regex=("._W|latitude_X|longitude_X|NROTARJETA")).columns]
    dfw.rename(columns={'FECHATRX_W': 'FECHATRX', 'idlinea_W': 'idlinea', 'latitude_W': 'latitude_Origen',
                        'longitude_W': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_W': 'CODIGOCONTRATO'}, inplace=True)

    out = [dfx, dfz, dfw, dfy]
    A = pd.concat(out)
    A["Tipificación"] = "n4dp"
    return A

def ExtraerMOD2pf(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-3):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) >= cotasuptiempo)  # no son prestada ni combin
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) >= cotasuptiempo) # no son prestada ni combi
        prop3 = ((dataframe.at[index + 3, "minutototal"] - dataframe.at[index + 2, "minutototal"]) <= cotainftiempo) # Z y W son la misma
        prop4 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[
            index + 1, "NROTARJETA"])  # Tiene que tratarse de la misma tarjeta
        prop5 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[
            index + 2, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        prop6 = dataframe.at[index + 2, "NROTARJETA"] == dataframe.at[
            index + 3, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        if prop4 and prop1 and prop3 and prop2 and prop6 and prop5:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            dictaux4 = {"NROTARJETA": dataframe.at[index + 3, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 3, "FECHATRX"],
                        "longitude": dataframe.at[index + 3, "longitude"],
                        "latitude": dataframe.at[index + 3, "latitude"], "idlinea": dataframe.at[index + 3, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 3, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
            copiainterna.append(dictaux4)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X", 'variable6': 'CODIGOCONTRATO_Y',
                              'variable7': 'FECHATRX_Y',
                              "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y", 'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z", 'variable16': 'CODIGOCONTRATO_W',
                              'variable17': 'FECHATRX_W', "variable18": "idlinea_W",
                              "variable19": "latitude_W", "variable20": "longitude_W"})
    return df1


def ExportarMOD42pf(dataframe):
    B = ExtraerMOD2pf(dataframe).reset_index()

    # Distancia entre linea inicial y segundo punto. Coloquialmente: quiero saber si el segundo punto corresponde a un destino inicial valido
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_X"] == l]  # filtro la linea media
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z", "FECHATRX_W"],
                              keep="first", inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Me volvi al inicio?
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Y"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    G = []
    # me dirijo al cuarto punto?
    for l in Lineas:
        dfinterno = F[F["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_W'], dfinterno['latitude_W']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconW")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    I = []
    for l in Lineas:
        dfinterno = H[H["idlinea_W"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLWconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        I.append(cruce)

    J = pd.concat(I)
    J.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    conditions = [
        (J["DistanciaLXconY"] <= distanciatol) & (J["DistanciaLYconZ"] <= distanciatol) & (J["DistanciaLZconW"] <= distanciatol) & (
                    J["DistanciaLWconX"] <= distanciatol)]
    valores = ["n4dpf"]
    J["Tipificación"] = np.select(conditions, valores)

    n4dpi = J[J["Tipificación"] == "n4dpf"]

    dfx = n4dpi[n4dpi.filter(regex=("._X|latitude_Y|longitude_Y|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Y': 'latitude_Destino',
                        'longitude_Y': 'longitude_Destino', 'CODIGOCONTRATO_X': 'CODIGOCONTRATO'}, inplace=True)

    dfy = n4dpi[n4dpi.filter(regex=("._Y|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', 'CODIGOCONTRATO_Y': 'CODIGOCONTRATO'}, inplace=True)

    dfz = n4dpi[n4dpi.filter(regex=("._Z|latitude_X|longitude_X|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'}, inplace=True)

    dfw = n4dpi[n4dpi.filter(regex=("._W|latitude_X|longitude_X|NROTARJETA")).columns]
    dfw.rename(columns={'FECHATRX_W': 'FECHATRX', 'idlinea_W': 'idlinea', 'latitude_W': 'latitude_Origen',
                        'longitude_W': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_W': 'CODIGOCONTRATO'}, inplace=True)

    out = [dfx, dfz, dfw, dfy]
    A = pd.concat(out)
    A["Tipificación"] = "n4dpf"
    return A

def ExtraerMOD2pi(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-3):
        prop1 = ((dataframe.at[index + 1, "minutototal"] - dataframe.at[index, "minutototal"]) <= cotainftiempo) # comparto tarjeta
        prop2 = ((dataframe.at[index + 2, "minutototal"] - dataframe.at[index + 1, "minutototal"]) >= cotasuptiempo) # ni comb ni compartida
        prop3 = ((dataframe.at[index + 3, "minutototal"] - dataframe.at[index + 2, "minutototal"]) >= cotasuptiempo) # ni comb ni compartida
        prop4 = (dataframe.at[index, "NROTARJETA"] == dataframe.at[
            index + 1, "NROTARJETA"])  # Tiene que tratarse de la misma tarjeta
        prop5 = dataframe.at[index + 1, "NROTARJETA"] == dataframe.at[
            index + 2, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        prop6 = dataframe.at[index + 2, "NROTARJETA"] == dataframe.at[
            index + 3, "NROTARJETA"]  # Tiene que tratarse de la misma tarjeta
        if prop4 and prop1 and prop3 and prop2 and prop6 and prop5:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            dictaux3 = {"NROTARJETA": dataframe.at[index + 2, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 2, "FECHATRX"],
                        "longitude": dataframe.at[index + 2, "longitude"],
                        "latitude": dataframe.at[index + 2, "latitude"], "idlinea": dataframe.at[index + 2, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 2, "CODIGOCONTRATO"]}
            dictaux4 = {"NROTARJETA": dataframe.at[index + 3, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 3, "FECHATRX"],
                        "longitude": dataframe.at[index + 3, "longitude"],
                        "latitude": dataframe.at[index + 3, "latitude"], "idlinea": dataframe.at[index + 3, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 3, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
            copiainterna.append(dictaux3)
            copiainterna.append(dictaux4)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    df1 = (copiainterna.set_index(['NROTARJETA', copiainterna.groupby('NROTARJETA').cumcount()])
           .unstack()
           .sort_index(axis=1, level=1))
    df1.columns = [f'variable{x}' for x in range(1, len(df1.columns) + 1)]
    df1 = df1.rename(columns={'variable1': 'CODIGOCONTRATO_X', 'variable2': 'FECHATRX_X', "variable3": "idlinea_X",
                              "variable4": "latitude_X", "variable5": "longitude_X", 'variable6': 'CODIGOCONTRATO_Y',
                              'variable7': 'FECHATRX_Y',
                              "variable8": "idlinea_Y", "variable9": "latitude_Y",
                              "variable10": "longitude_Y", 'variable11': 'CODIGOCONTRATO_Z',
                              'variable12': 'FECHATRX_Z', "variable13": "idlinea_Z",
                              "variable14": "latitude_Z", "variable15": "longitude_Z", 'variable16': 'CODIGOCONTRATO_W',
                              'variable17': 'FECHATRX_W', "variable18": "idlinea_W",
                              "variable19": "latitude_W", "variable20": "longitude_W"})
    return df1


def ExportarMOD42pi(dataframe):
    B = ExtraerMOD2pi(dataframe).reset_index()

    # Distancia entre linea inicial y segundo punto. Coloquialmente: quiero saber si el segundo punto corresponde a un destino inicial valido
    C = []
    for l in Lineas:
        dfinterno = B[B["idlinea_X"] == l]  # filtro la linea media
        a = gpd.read_file(direccionshp + str(l) + ".shp",
                          crs=4326)  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Y'], dfinterno['latitude_Y']),
                                      crs="EPSG:4326").to_crs(epsg=22185)
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLXconY")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z", "FECHATRX_W"],
                              keep="first", inplace=True)
        C.append(cruce)

    D = pd.concat(C)
    D.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    # Me volvi al inicio?
    E = []
    for l in Lineas:
        dfinterno = D[D["idlinea_Y"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_Z'], dfinterno['latitude_Z']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLYconZ")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        E.append(cruce)

    F = pd.concat(E)
    F.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    G = []
    # me dirijo al cuarto punto?
    for l in Lineas:
        dfinterno = F[F["idlinea_Z"] == l]  # filtro la linea final
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_W'], dfinterno['latitude_W']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLZconW")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        G.append(cruce)

    H = pd.concat(G)
    H.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)

    I = []
    for l in Lineas:
        dfinterno = H[H["idlinea_W"] == l]
        a = gpd.read_file(direccionshp + str(l) + ".shp")  # cargo el shp de la linea
        a = a[a["geometry"] != None]
        a = a.to_crs(epsg=22185)
        a = a.reset_index()
        a = a[["Linea", "geometry"]]
        shpinterno = gpd.GeoDataFrame(dfinterno,
                                      geometry=gpd.points_from_xy(dfinterno['longitude_X'], dfinterno['latitude_X']),
                                      crs="4326").to_crs(epsg=22185)  # Armo un geodataframe con mi informacion
        cruce = gpd.sjoin_nearest(shpinterno, a, how="left", distance_col="DistanciaLWconX")
        cruce.drop_duplicates(subset=["NROTARJETA", "FECHATRX_X", "FECHATRX_Y", "FECHATRX_Z"], keep="first",
                              inplace=True)
        I.append(cruce)

    J = pd.concat(I)
    J.drop(["Linea", "index_right", "geometry"], axis=1, inplace=True)
    conditions = [
        (J["DistanciaLYconZ"] <= distanciatol) & (J["DistanciaLZconW"] <= distanciatol) & (J["DistanciaLWconX"] <= distanciatol)]
    valores = ["n4dpi"]
    J["Tipificación"] = np.select(conditions, valores)

    n4dpi = J[J["Tipificación"] == "n4dpi"]

    dfx = n4dpi[n4dpi.filter(regex=("._X|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfx.rename(columns={'FECHATRX_X': 'FECHATRX', 'idlinea_X': 'idlinea', 'latitude_X': 'latitude_Origen',
                        'longitude_X': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', 'CODIGOCONTRATO_X': 'CODIGOCONTRATO'}, inplace=True)

    dfy = n4dpi[n4dpi.filter(regex=("._Y|latitude_Z|longitude_Z|NROTARJETA")).columns]
    dfy.rename(columns={'FECHATRX_Y': 'FECHATRX', 'idlinea_Y': 'idlinea', 'latitude_Y': 'latitude_Origen',
                        'longitude_Y': 'longitude_Origen', 'latitude_Z': 'latitude_Destino',
                        'longitude_Z': 'longitude_Destino', 'CODIGOCONTRATO_Y': 'CODIGOCONTRATO'}, inplace=True)

    dfz = n4dpi[n4dpi.filter(regex=("._Z|latitude_W|longitude_W|NROTARJETA")).columns]
    dfz.rename(columns={'FECHATRX_Z': 'FECHATRX', 'idlinea_Z': 'idlinea', 'latitude_Z': 'latitude_Origen',
                        'longitude_Z': 'longitude_Origen', 'latitude_W': 'latitude_Destino',
                        'longitude_W': 'longitude_Destino', 'CODIGOCONTRATO_Z': 'CODIGOCONTRATO'}, inplace=True)

    dfw = n4dpi[n4dpi.filter(regex=("._W|latitude_X|longitude_X|NROTARJETA")).columns]
    dfw.rename(columns={'FECHATRX_W': 'FECHATRX', 'idlinea_W': 'idlinea', 'latitude_W': 'latitude_Origen',
                        'longitude_W': 'longitude_Origen', 'latitude_X': 'latitude_Destino',
                        'longitude_X': 'longitude_Destino', 'CODIGOCONTRATO_W': 'CODIGOCONTRATO'}, inplace=True)

    out = [dfx, dfz, dfw, dfy]
    A = pd.concat(out)
    A["Tipificación"] = "n4dpi"
    return A

def ExtraerMODn2c(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-1):
        prop1 = dataframe.at[index+1,"idlinea"] != dataframe.at[index,"idlinea"] #Distinta linea
        prop2 = (dataframe.at[index+1,"minutototal"] - dataframe.at[index,"minutototal"]) > cotainftiempo  # no es tarjeta compartida
        prop3 = dataframe.at[index+1,"NROTARJETA"] == dataframe.at[index,"NROTARJETA"] #Misma tarjeta
        prop4 = (dataframe.at[index+1,"minutototal"] - dataframe.at[index,"minutototal"]) < cotasuptiempo # si o si combinacion
        if prop1 and prop2 and prop3 and prop4:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}

            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    a1 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="first")
    a2 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="last")
    a3 = a1.merge(a2, on='NROTARJETA', how='inner', suffixes=("_Inicial", "_Final"))
    return a3


def ExportarMODn2c(dataframe):
    B = ExtraerMODn2c(dataframe)
    C = B[['NROTARJETA', 'FECHATRX_Inicial', 'longitude_Inicial',
           'latitude_Inicial', 'idlinea_Inicial', 'CODIGOCONTRATO_Inicial', 'longitude_Final', 'latitude_Final', 'idlinea_Final']]
    C.rename(
        columns={'FECHATRX_Inicial': 'FECHATRX', 'idlinea_Inicial': 'idlinea', 'latitude_Inicial': 'latitude_Origen',
                 'longitude_Inicial': 'longitude_Origen', 'CODIGOCONTRATO_Inicial': 'CODIGOCONTRATO', 'latitude_Final' : "latitude_Destino",  'longitude_Final' : "longitude_Destino", "idlinea_Final" : "LineaCombinacion" }, inplace=True)
    C["Tipificación"] = "n2c"
    return C


def obtenern2():
    A = []
    for d in dias:
        print("Analizando dia " + str(d))
        #casosn2
        aux = dffiltrado[dffiltrado["DIA"] == d]
        #N2
        n2 = casosn(aux,2)
        aml = ExportarMODml(n2)
        adl = ExportarMODdl(n2)
        A.append(aml)
        A.append(adl)               # los meto en una lista
        try:
            acomb = ExportarMODn2c(n2)
        except KeyError:
            print("No hay casos de n2c en el dia " + str(d))
            acomb = pd.DataFrame()
            A.append(acomb)
    A = pd.concat(A)
    return A                 # lista con todos los casos de n2 ya tipificados

start = time.time()
n2 = obtenern2()
conteon2 = len(n2)
def obtenern3n4():
    A = []
    for d in dias:
        print("Analizando dia " + str(d))
        aux = dffiltrado[dffiltrado["DIA"] == d]
        #N3
        n3 = casosn(aux,3)
        an3 = ExportarMODn3(n3)
        A.append(an3)
        try:
            an3ci = ExportarMODn3ci(n3)
            A.append(an3ci)
        except KeyError:
            print("No hay casos de n3ci en el dia " + str(d))
            an3ci = pd.DataFrame()
            A.append(an3ci)
        try:
            an3cf = ExportarMODn3cf(n3)
            A.append(an3cf)
        except KeyError:
            print("No hay casos de n3cf en el dia " + str(d))
            an3cf = pd.DataFrame()
            A.append(an3cf)
        try:
            an3dpf = ExportarMODdpf(n3)
            A.append(an3dpf)
        except KeyError:
            print("No hay casos de n3dpf en el dia " + str(d))
            an3dpf = pd.DataFrame()
            A.append(an3dpf)
        try:
            an3dpi = ExportarMODn3pi(n3)
            A.append(an3dpi)
        except KeyError:
            print("No hay casos de n3dpi en el dia " + str(d))
            an3dpi = pd.DataFrame()
            A.append(an3dpi)
        #N4
        n4 = casosn(aux,4)
        try:
            an4s = ExportarMODn4s(n4)
            A.append(an4s)
        except KeyError:
            print("No hay casos de n4s en el dia " + str(d))
            an4s = pd.DataFrame()
            A.append(an4s)
        an4dp = ExportarMOD4dp(n4)
        A.append(an4dp)
        try:
            an4dpf = ExportarMOD42pi(n4)
            A.append(an4dpf)
        except KeyError:
            print("No hay casos de n4dpf en el dia " + str(d))
            an4dpf = pd.DataFrame()
            A.append(an4dpf)
        try:
            an4ci = ExportarMODn4c(n4)
            A.append(an4ci)
        except KeyError:
            print("No hay casos de n4c en el dia " + str(d))
            an4ci = pd.DataFrame()
            A.append(an4ci)
    return pd.concat(A)

n3n4 = obtenern3n4()
conteon3n4 = len(n3n4)
B = pd.concat([n2, n3n4])
tiempon234 = (time.time() - inicio)/60
print("n2, n3, n4 obtenidos,", "Tiempo de procesamiento:",str(tiempon234),"minutos")


def ExtraerMODn1den2(dataframe):
    """dataframe tiene que estar ordenado por nro de tarjeta y fecha"""
    dataframe['minutototal'] = dataframe.apply(lambda row: minutototal(row), axis = 1)
    copiainterna = []
    for index in range(len(dataframe)-1):
        prop1 = dataframe.at[index+1,"idlinea"] == dataframe.at[index,"idlinea"] #Misma linea
        prop2 = (dataframe.at[index+1,"minutototal"] - dataframe.at[index,"minutototal"]) <= cotainftiempo
        prop3 = dataframe.at[index+1,"NROTARJETA"] == dataframe.at[index,"NROTARJETA"] #Misma tarjeta
        if prop1 and prop2 and prop3:
            dictaux1 = {"NROTARJETA": dataframe.at[index, "NROTARJETA"], "FECHATRX": dataframe.at[index, "FECHATRX"],
                        "longitude": dataframe.at[index, "longitude"], "latitude": dataframe.at[index, "latitude"],
                        "idlinea": dataframe.at[index, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index, "CODIGOCONTRATO"]}
            dictaux2 = {"NROTARJETA": dataframe.at[index + 1, "NROTARJETA"],
                        "FECHATRX": dataframe.at[index + 1, "FECHATRX"],
                        "longitude": dataframe.at[index + 1, "longitude"],
                        "latitude": dataframe.at[index + 1, "latitude"], "idlinea": dataframe.at[index + 1, "idlinea"],
                        "CODIGOCONTRATO": dataframe.at[index + 1, "CODIGOCONTRATO"]}
            copiainterna.append(dictaux1)
            copiainterna.append(dictaux2)
        else:
            continue
    copiainterna = pd.DataFrame(copiainterna)
    copiainterna = copiainterna[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    a1 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="first")
    a2 = copiainterna.drop_duplicates(subset=["NROTARJETA"], keep="last")
    a3 = a1.merge(a2, on='NROTARJETA', how='inner', suffixes=("_Inicial", "_Final"))
    return a3


def ExportarMODn1den2(dataframe):
    B = dataframe
    C = B[['NROTARJETA', 'FECHATRX_Inicial', 'longitude_Inicial',
           'latitude_Inicial', 'idlinea_Inicial', 'CODIGOCONTRATO_Inicial']]
    D = B[["NROTARJETA", 'FECHATRX_Final', 'longitude_Final', 'latitude_Final', 'idlinea_Final',
           'CODIGOCONTRATO_Final']]

    C.rename(
        columns={'FECHATRX_Inicial': 'FECHATRX', 'idlinea_Inicial': 'idlinea', 'latitude_Inicial': 'latitude_Origen',
                 'longitude_Inicial': 'longitude_Origen', 'CODIGOCONTRATO_Inicial': 'CODIGOCONTRATO'}, inplace=True)
    D.rename(columns={'FECHATRX_Final': 'FECHATRX', 'idlinea_Final': 'idlinea', 'latitude_Final': 'latitude_Origen',
                      'longitude_Final': 'longitude_Origen', 'CODIGOCONTRATO_Final': 'CODIGOCONTRATO'}, inplace=True)
        # ambos van con origen y final
    E = pd.concat([C, D])
    return E


def obtenern1():
    n1aux= pd.DataFrame()
    n1n2aux = pd.DataFrame()
    for d in dias:
        auxiliar = dffiltrado[dffiltrado["DIA"] == d]
        auxiliarn2n1 = casosn(auxiliar,2)
        auxiliarn1 = casosn(auxiliar,1)
        n1aux = pd.concat([n1aux, auxiliarn1])
        try:
            auxn2n1 = ExtraerMODn1den2(auxiliarn2n1)
            auxn2n1 = ExportarMODn1den2(auxn2n1)
            n1n2aux = pd.concat([n1n2aux, auxn2n1])       ### los que van de a 2 pero son 1
        except ValueError:
            continue

    n1 = n1aux
    n1 = n1[["NROTARJETA", "FECHATRX", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
    n1.rename(columns={'latitude': 'latitude_Origen',
                        'longitude': 'longitude_Origen'}, inplace=True)
    n1n2 = n1n2aux
    n1 = pd.concat([n1,n1n2])
    return n1

n1 = obtenern1()
conteon1 = len(n1)


def predecirn1(dfn1, dfnx):
    #aca el shape es particular y se hacen procesamiento particulares
    SHzona = gpd.read_file(direccionzonificacion)
    try:
        SHzona[zona]= SHzona.apply(lambda x:  str(x['frac']) +'-'+ str(x['radio']) ,axis=1)#### agregado para santa fe
    except KeyError:
        print("Error: Columnas erróneas en zonificación")
    SHzona=SHzona.drop_duplicates(subset = zona)
    #SHzona.rename(columns={'RC': "Vecinales"}, inplace=True) # no harÃ­a mÃ¡s falta
    SHzona= SHzona.to_crs(crs="4326")
    SHzona["centroid"]=SHzona.representative_point()
    SHzona["long"] = SHzona.centroid.map(lambda p: p.x)
    SHzona["lat"]  = SHzona.centroid.map(lambda p: p.y)
    dfnx = dfnx.reset_index()
    A = pd.DataFrame()
    A2 = pd.DataFrame()
    dfnx["HORA"] = dfnx["FECHATRX"].dt.hour
    for l in Lineas:
        auxiliar = dfnx[dfnx["idlinea"] == l]
        puntosZonaIda = gpd.GeoDataFrame(auxiliar, geometry=gpd.points_from_xy(auxiliar['longitude_Origen'],
                                                                                  auxiliar['latitude_Origen']),
                                            crs="4326")
        puntosEnZonaIda = gpd.sjoin_nearest(puntosZonaIda, SHzona, how="inner")
        puntosZonaVuelta = gpd.GeoDataFrame(auxiliar, geometry=gpd.points_from_xy(auxiliar['longitude_Destino'],
                                                                                     auxiliar['latitude_Destino']),
                                               crs="4326")
        puntosEnZonaVuelta = gpd.sjoin_nearest(puntosZonaVuelta, SHzona, how="inner")
        puntosEnIdayVuelta = puntosEnZonaIda.merge(puntosEnZonaVuelta, on="index")
        puntosEnIdayVuelta = puntosEnIdayVuelta[
            [zona+'_x', zona+'_y', "idlinea_x", "CODIGOCONTRATO_x", "HORA_x", "NROTARJETA_x", "FECHATRX_x"]]
        puntosEnIdayVuelta = puntosEnIdayVuelta.drop_duplicates(subset=["idlinea_x", "NROTARJETA_x", "HORA_x", "FECHATRX_x"],keep='first')
        prueba1 = puntosEnIdayVuelta.groupby([zona+'_x', zona+'_y', "idlinea_x", "HORA_x"]).size().to_frame(
            'Cuenta').reset_index()
        A = pd.concat([A, prueba1])
        A2 = pd.concat([A2, puntosEnIdayVuelta])
    zonasnx = A                 ## cuenta de vecinales origenes y destinos para vecinales nx
    B2 = A2                         ## puntos de ida y vuelta de transacciones tipo nx

    idn1 = dfn1["NROTARJETA"].unique()  # me quedo con las tarjetas de las n1
    idnx = B2["NROTARJETA_x"].unique()  # y las tarjetas nx

    n1aux = dfn1[dfn1["NROTARJETA"].isin(idnx)] # si tengo una n1 que tambien ha hecho transacciones nx
    puntosn1 = gpd.GeoDataFrame(n1aux, geometry=gpd.points_from_xy(n1aux['longitude_Origen'],
                                                                             n1aux['latitude_Origen']), crs="4326")
    n1aux = gpd.sjoin_nearest(puntosn1, SHzona, how="inner")  # me agrega a cada punto una coolumna vecinal con la vecinal mas cercana
    n1aux = n1aux[["NROTARJETA", "FECHATRX", "idlinea", zona, "CODIGOCONTRATO"]]
    n1aux.drop_duplicates(subset=["NROTARJETA", "FECHATRX"], inplace=True)
    n1aux["CODIGO"] = n1aux[zona]
    B2.rename(columns = {"NROTARJETA_x" : "NROTARJETA"}, inplace = True) #cambio el nombre por uno mas feliz

    nxutiles = B2[B2["NROTARJETA"].isin(idn1)] # me quedo con las tarjetas de nx que estÃ©n en n1

    AA = pd.DataFrame()
    for i in n1aux["NROTARJETA"].unique():
        aux = nxutiles[nxutiles["NROTARJETA"] == i] #me centro en un nro de tarjeta
        aux2 = aux.copy(deep=True)
        aux2[[zona+'_x', zona+'_y']] =  aux2[[zona+'_y', zona+'_x']]
        aux3 = pd.concat([aux,aux2]) # idas y vueltas identicas
        aux4 = aux3.groupby([zona+'_x', zona+'_y']).size().to_frame("Cuenta").reset_index() ## cuenta cuantas veces voy de x a y
        aux5 = aux4.groupby([zona+'_x']).sum().reset_index()  # cuantas veces partÃ­ de x
        aux5 = aux5[[zona+'_x','Cuenta']]  # Quitado la columna RC_y ###
        aux5 = aux5.rename(columns= {"Cuenta" : "Total"})
        aux6 = pd.merge(aux4, aux5, how="left", on=[zona+'_x']) # de todas la veces que parti de x cuantas fui a y ### MODIFICADO on=
        aux6["p"] = aux6["Cuenta"] / aux6["Total"]  # proporcion de veces que voy de x a y
        aux6["CODIGO"] = aux6[zona+'_x']          # cambio nombre a codigo
        filtron1 = n1aux[n1aux["NROTARJETA"] == i]  # busco esta tarjeta en n1aux
        codigospart = filtron1["CODIGO"].unique()   #desde que vecinal salgo
        A = pd.DataFrame()
        t = 0
        for c in codigospart:
            t += 1
            aux7 = filtron1[filtron1["CODIGO"] == c]  # desde donde parto
            prob = aux6[aux6["CODIGO"] == c].sort_values(by='Cuenta', ascending=True).reset_index(drop=True) # la probabilidad
            if len(prob) > 0:
                totals = prob.Cuenta.cumsum().tolist()  # acumulada de cantidad de veces
                lista = prob[zona+'_y'].tolist()   # cuales son esos destinos
                a = random.choices(lista, cum_weights=totals, k=len(aux7)) # elijo con el criterio acumulado de cantidad de vces
                aux7[zona+'_y'] = a   #  le asigno esta elegida como destino de mi n1
                A = pd.concat([A, aux7])  #
            else:
                continue
        try:
            P = A
        except ValueError:
            continue
        AA = pd.concat([AA, P])
    BB = AA ## BB son todas las transacciones n1 transformadas a n2 con sus destinos probablisticamente elegidos
    zonasn1 = BB[["NROTARJETA", "FECHATRX", zona+'_y',zona, "idlinea", "CODIGOCONTRATO" ]]
    zonasn1["FECHATRX"] = pd.to_datetime(zonasn1["FECHATRX"])  # las zonas n1  obtenidas de BB
    ## ahora tengo en cuenta la hora de la transaccion
    C = zonasnx.groupby([zona+'_x', "idlinea_x", "HORA_x"]).sum().reset_index()
    C = C.rename(columns={'Cuenta': 'Total'})
    C = C[[zona+'_x','idlinea_x','HORA_x','Total']]
    D = pd.merge(zonasnx, C, on=[zona+'_x', "idlinea_x", "HORA_x"], how="left")

    D["p"] = D["Cuenta"] / D["Total"]

    auxiliar = zonasn1[["NROTARJETA", "FECHATRX"]]     ## vuevlo con estas transacciones ya asignadas a zonas
    auxiliar = dfn1.merge(auxiliar, on=["NROTARJETA", "FECHATRX"], how='outer', indicator=True)  # #
    n1paraprobas = auxiliar[auxiliar["_merge"] == "left_only"]  # ME QUEDO CON LAS Q NO ANALICE YA
    puntosn1 = gpd.GeoDataFrame(n1paraprobas, geometry=gpd.points_from_xy(n1paraprobas['longitude_Origen'],
                                                                   n1paraprobas['latitude_Origen']), crs="4326")
    n1aux = gpd.sjoin_nearest(puntosn1, SHzona, how="inner")
    n1aux = n1aux[["NROTARJETA", "FECHATRX", "idlinea", zona, "CODIGOCONTRATO" ]]
    n1aux.drop_duplicates(subset=["NROTARJETA", "FECHATRX"], inplace=True)

    n1aux["HORA"] =  n1aux["FECHATRX"].dt.hour
    n1aux["CODIGO"] = n1aux[zona] + "-" + n1aux["idlinea"].astype(str) + "-" + "h" + n1aux["HORA"].astype(str)
    D["CODIGO"] = D[zona+'_x'] + "-" + D["idlinea_x"].astype(str) + "-" + "h" + D["HORA_x"].astype(str)

    codigos = n1aux["CODIGO"].unique()
    A = pd.DataFrame()
    t = 0
    # codigo para distribuir por probabilidad particular
    for c in codigos:
        t += 1
        aux = n1aux[n1aux["CODIGO"] == c]
        prob = D[D["CODIGO"] == c].sort_values(by='Cuenta', ascending=True).reset_index(drop=True)
        if len(prob) > 0:
            totals = prob.Cuenta.cumsum().tolist()
            lista = prob[zona+'_y'].tolist()
            a = random.choices(lista, cum_weights=totals, k=len(aux))
            aux[zona+'_y'] = a
            A = pd.concat([A, aux])
        else:
            continue
    P = A
    P = P.rename(columns={zona: zona+'_x'})
    P = P.drop("CODIGO", axis=1)
    zonasn1["HORA"] = zonasn1["FECHATRX"].dt.hour
    zonasn1 = zonasn1.rename(columns={zona: zona+'_x'})
    zonasnx = zonasnx.rename(columns={'idlinea_x': 'idlinea', "HORA_x": "HORA"})
    zonasn1["Tipificación"] = "n1log"
    P["Tipificación"] = "n1proba"
    n1tot = pd.concat([zonasn1, P])
    SHzona = SHzona.rename(columns={zona: zona+'_x'})
    n1tot = pd.merge(n1tot, SHzona, on=zona+'_x', how="left")
    SHzona = SHzona.rename(columns={zona+'_x': zona+'_y'}) # para que coincida con las n1tot
    n1tot = pd.merge(n1tot, SHzona, on=zona+'_y', how="left") # para poder ubicar las transacciones
    n1tot = n1tot.drop(["geometry_x", "centroid_x", "geometry_y", "centroid_y", zona+'_x', zona+'_y', "HORA"], axis=1)
    n1tot = n1tot.rename(columns={'long_x': 'longitude_Origen','long_y': 'longitude_Destino', 'lat_x': 'latitude_Origen', 'lat_y': 'latitude_Destino' }) # asigno origen y destino a c/u
    return (n1tot,zonasnx)







n1d, nx = predecirn1(n1,B)
odtotal = pd.concat([n1d,B])
odtotal.to_csv(direcciondestinood)

##########################################################################################################################
#ESTADISTICOS
#trx perdidas por combinacion
listacomb = ["n3ci", "n3cf", "n4combinacion", "n2c"]
condicion = [odtotal["Tipificación"].isin(listacomb)]
valor  = [1]
odtotal["Combinacion"] = np.select(condicion, valor, default=0)



enes= [5,6,7,8,9,10]
PC = []
for d in dias:
    aux = dffiltrado[dffiltrado["DIA"] == d]
    for n in enes:
        K = casosn(aux,n) # casos n en el dÃ­a d
        PC.append(len(K)) # cantidad de transacciones de caso n por dÃ­a
sumaPC = sum(PC) # transacciones perdidas por no tomar casos n= 5,6,7,8,9,10,11,12

perdidasporcombinacion = odtotal.Combinacion.sum()
perdidatolerancia = sum(perdidatolerancia) #trx perdidas por tolerancia superada
perdidaconocida= perdidatolerancia+perdidasporcombinacion + perdidagps # modif
perdidatotal= len(dffiltradoOriginal) - len(odtotal)


cobertura = ((len(odtotal) + perdidasporcombinacion) / len(dffiltrado)) * 100
print("La OD generada explica al menos el " + str(round(cobertura)) + "% de las transacciones")


tiempo= (time.time() - inicio)/60
print('Tiempo de procesamiento: ' + str(tiempo) + ' minutos')

#estadisticas por linea
ODporlinea = odtotal.groupby("idlinea").size().to_frame("CuentaOD").reset_index()
ODporlinea = ODporlinea[["idlinea", "CuentaOD"]]
dfporlinea = dffiltradoOriginal.groupby("idlinea").size().to_frame("CuentaDF").reset_index() #modif # Estoy contando las trx incluso de las que tenian lat,long 0
dfporlinea = dfporlinea[["idlinea", "CuentaDF"]]

# Reporte
NDatos =['Transacciones Procesadas', 'N1','N2','N3-N4','Transacciones Perdidas','%Cobertura','Nº Combinaciones','% Combinaciones' , 'Pérdidas por tolerancia', '%Pérdidas por tolerancia','Pérdidas casos N5 a N10', 'Tiempo de Procesamiento']
Datos = [len(dffiltrado),conteon1,conteon2,conteon3n4,perdidatotal,cobertura,perdidasporcombinacion*2, round((perdidasporcombinacion*2/len(dffiltrado) *100),2),perdidatolerancia,round(perdidatolerancia/len(dffiltrado)*100),sumaPC,tiempo]
dfReporte = pd.DataFrame(index=NDatos,data=Datos,columns=['Datos'])
dfReporte.to_excel(direcciondestinoreporte, index=True)



t1 = pd.merge(ODporlinea, dfporlinea, on="idlinea", how="inner")
t1["% por linea"] = (t1["CuentaOD"] / t1["CuentaDF"])*100
t1["Factor"] = (t1["CuentaOD"] / t1["CuentaDF"])

for i in range(len(t1)):
    if t1.at[i,"Factor"] < 0.7:
        print("Atención: Línea "+str(t1.at[i,"idlinea"])+ " Tiene un factor de " +str(t1.at[i,"Factor"]))


t1.to_excel(direcciondestinot1)
