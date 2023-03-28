import pandas as pd
#Cargo los datos de gps
# gps1 = pd.read_csv("Rawdata/gpsEnero.csv", sep=";", parse_dates=['date_time'], dayfirst=True)
# gps1 = gps1[gps1['date_time']>=pd.to_datetime('01-01-2022',format="%d-%m-%Y")]
#
#
# gps2 = pd.read_csv("Rawdata/gpsMarzo.csv", sep=";", parse_dates=['date_time'], dayfirst=True)
# gps2 = gps2[gps2['date_time']>=pd.to_datetime('01-01-2022',format="%d-%m-%Y")]
#
# gps3 = pd.read_csv("Rawdata/gpsAbril.csv", sep=";", parse_dates=['date_time'], dayfirst=True)
# gps3 = gps3[gps3['date_time']>=pd.to_datetime('01-01-2022',format="%d-%m-%Y")]
#
# gps4 = pd.read_csv("Rawdata/gpsMayo.csv", sep=";", parse_dates=['date_time'], dayfirst=True)
# gps4 = gps4[gps4['date_time']>=pd.to_datetime('01-01-2022',format="%d-%m-%Y")]
#
# gps5 = pd.read_csv("Rawdata/gpsJunio.csv", sep=";", parse_dates=['date_time'], dayfirst=True)
# gps5 = gps5[gps5['date_time']>=pd.to_datetime('01-01-2022',format="%d-%m-%Y")]
#
# gps6 = pd.read_csv("Rawdata/gpsJulio.csv", sep=";", parse_dates=['date_time'], dayfirst=True)
# gps6 = gps6[gps6['date_time']>=pd.to_datetime('01-01-2022',format="%d-%m-%Y")]
#
# lista = [gps1,gps2,gps3,gps4,gps5,gps6]
# gpstotal = pd.concat(lista)
# bool_series = gpstotal.duplicated(keep='first')
# gpstotal=gpstotal[~bool_series]
# # gpstotal.to_csv("gpstotal.csv")


#Leo el archivo con todos los gps y un mes en particular y algo el merge para ese mes
######################################################################################
ts = pd.read_csv("Rawdata/transJunio.csv", sep=";", parse_dates=['FECHATRX'], dayfirst=True)
gps = pd.read_csv("Preprocessdata/gpstotal.csv", sep=",", parse_dates=['date_time'], dayfirst=True)

def crearkeyts(registro): #creo funciones para generar las key que se usan para matchear las tablas
    auxiliar = str(registro["IDARCHIVOINTERCAMBIO"]) + str(registro["ID_POSICIONAMIENTO"])
    return auxiliar

def crearkeygps(registro):
    auxiliar = str(registro["file_id"]) + str(registro["c_control_point"])
    return auxiliar

ts['keyunica'] = ts.apply(lambda row: crearkeyts(row), axis=1)
gps['keyunica'] = gps.apply(lambda row: crearkeygps(row), axis=1)

mergeJ= pd.merge(ts, gps, how="left", on=["keyunica"])
mergeJ = mergeJ.drop(['codigoentidad', "IDLINEA", "interno", "device", "dtsn"], axis=1) #estos atributos se repiten

df = mergeJ

df["FECHATRX"] = pd.to_datetime(df["FECHATRX"])
df['MES']=df['FECHATRX'].dt.month
df = df[df["MES"] == 11]
listaoriginal= [763,758,756,759,764,754,751,757,762,750,760,748,755,867,2704]
lista = [9,18,13,8,15,10,5,16,3,1,14,4,11,2,21]
df = df.rename(columns={"IDLINEA" : "idlinea"})
df["idlinea"] =df["idlinea"].replace(listaoriginal, lista).astype(int)

listatarifas= [602, 621, 681, 690, 691, 692, 693, 695, 696, 716, 717, 718, 719, 722, 749]
listatarifasnueva = ["Tarifa Plana","Atributo Nacional","Secundario Prioridad 1","Escolares","Terciario y Universitario","Seguro","Jubilados - Pensiones ley 5110","Empleado Municipal","Empleado Municipal Despacho","BEG Docente Provincial","BEG No Docente Provincial","BEG Primario Provincial","BEG Secundario Provincial","BEG Terciario Provincial","BEG Universitario Provincial"]
#lista = ["Tarifa Plana", "Atributo Nacional", "Est Sec", "Est Prim", "Est Ter/UNI", "Seguro", "Jubilado", "Emp Muni", "Emp Muni", "Docente", "No Docente", "Est Prim", "Est Sec", "Est Ter/UNI", "Est Ter/UNI" ] #Esta es una clasificacion auxiliar un poco mas compacta.
df["CODIGOCONTRATO"] =df["CODIGOCONTRATO"].replace(listatarifas, listatarifasnueva).astype(str)

da = df #DF va a ser la version resumida de los datos para procesamiento basico. #DA mantiene todas las columnas para procesamientos mas especificos
df = df[["NROTARJETA", "FECHATRX", "MES", "longitude", "latitude", "idlinea", "CODIGOCONTRATO"]]
# df.to_csv('Preprocessdata/DataTRXjuniodfbasico.csv')
# da.to_csv('Preprocessdata/DataTRXjuniodfcompleto.csv')


