import pandas as pd
import numpy as np
import openpyxl
import datetime
#### importo los datos
Ts_oct= pd.read_csv("10 Transacciones Octubre 2022.csv", sep= ';', parse_dates=["Fecha Hora"], dayfirst=True)
Ts_oct= Ts_oct[Ts_oct["Fecha Hora"]>= pd.to_datetime('1-10-2022', format='%d-%m-%Y')]
df=Ts_oct

# renombro y corrijo
df= df.rename(columns={"Fecha Hora":"FECHATRX"})
df= df.rename(columns={"Hora":"HORA"})
df= df.rename(columns={"Minutos":"MINUTO"})
df= df.rename(columns={"Segundos":"SEGUNDO"})
df["MES"]= df["FECHATRX"].dt.month
df["AÑO"]= df["FECHATRX"].dt.year
df["DIA"]=df["FECHATRX"].dt.day
df['FECHATRX'] = df.apply(lambda row: datetime.datetime(row["AÑO"], row["MES"], row["DIA"], row['HORA'], row['MINUTO'], row['SEGUNDO']), axis=1)
#df['FECHATRX'] = pd.to_datetime(df['DIA'].astype(str) + '/' + df['MES'].astype(str) + '/' + df['AÑO'].astype(str) + ' ' + df['HORA'].astype(str) + ':' + df['MINUTO'].astype(str) + ':' + df['SEGUNDO'].astype(str),format='%d/%m/%Y %H:%M:%S')


df= df[df["MES"]== 10]
## filtro los duplicados
print(df)
print(len(df))
df= df.drop_duplicates()
df= df
print(len(df))

da = df

### saco el ! a las líneas
def limp(numer):
    auxiliar=str(numer["Linea"])[1:]
    return auxiliar
df['LINEA']=df.apply(lambda row: limp(row), axis=1)

da = df
# me llevo lo que no se repite
df= df.rename(columns={"Ramal":"idlinea","Tarjeta":"NROTARJETA","Longitud":"longitude","Latitud":"latitude","Tipo Trx":"CODIGOCONTRATO"})
df= df[["NROTARJETA","FECHATRX","MES","latitude","longitude","idlinea","CODIGOCONTRATO"]]
da.to_csv("VT_OCT.csv")
df.to_csv("./Outputs/10-TX-SData.csv")
print('GENERADO VT_OCT.csv EN CARPETA OUTPUTS')
