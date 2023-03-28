###Script para generar un reporte mensual
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import geopandas as gpd
import contextily as cx
import folium
from folium import plugins

import numpy
from shapely.geometry import Point
import pygeos
import branca.colormap as cm
from collections import defaultdict
from jinja2 import Template
from folium.map import Layer

#CargarDFbasico
df = pd.read_csv('Preprocessdata/DataTRXnoviembre2021dfbasico.csv', sep=",", parse_dates=['FECHATRX'], dayfirst=True) #Leo el archivo con los datos
df = df.drop(["Unnamed: 0"], axis=1)
df["FECHATRX"] = pd.to_datetime(df["FECHATRX"])
df["MINUTO"]=df['FECHATRX'].dt.minute
df['HORA']=df['FECHATRX'].dt.hour
df['DIA']=df['FECHATRX'].dt.day
df['nDIA']=df['FECHATRX'].dt.day_name() #le agrego los dias por nombre
df['nDIA']=df['FECHATRX'].dt.day_name() #le agrego los dias por nombre
df = df[df["idlinea"] != 21]
df = df[df["idlinea"] != 761]

#Filtrar por tipo de dia. Esto no se puede automatizar debido a los feriados
dfhabil = df[df['nDIA'].isin(['Wednesday', 'Thursday', 'Friday', 'Monday',
       'Tuesday'])] #filtro los dias habiles a priori
dfhabil = dfhabil[dfhabil["DIA"] != 20] #quito el 20 por ser feriado
dfsabados = df[df['nDIA'].isin(["Saturday"])] #filtro los sabados
dfdomfer = df[df['nDIA'].isin(["Sunday", "Monday"])]
dfdomfer = dfdomfer[~dfdomfer['DIA'].isin([6,13,27])]

#EstadisticasBasicas

#LINEAS
pptPar=dfhabil.groupby(['idlinea', "DIA"]).size().to_frame('NumerodeTrans').reset_index()
thabiles = len(pptPar["DIA"].unique())
PromedioHabiles=pptPar.groupby(["idlinea"]).apply(lambda NumerodeTrans: (NumerodeTrans.sum()/thabiles)).sort_values(by = "NumerodeTrans", ascending=False)
PromedioHabiles=PromedioHabiles.drop(["DIA"], axis=1)
PromedioHabiles = PromedioHabiles.astype({'NumerodeTrans':'int'}) #Esto porque estaban como float
PromedioHabiles = PromedioHabiles.astype({'idlinea':'int'}) #Esto porque estaban como float

#si se cambia dfhabil por otro de los df generados se obtiene la distribución en el uso de la linea para otro tipo de dia.

#grafico asociado
n = PromedioHabiles["NumerodeTrans"].sum()

factorinterno = 100 / n
fig, ax = plt.subplots(figsize=(16,9))
ax.barh(PromedioHabiles["idlinea"].astype(str), PromedioHabiles["NumerodeTrans"], color=["#12bc8e"])
for c in ax.containers:
   labels = [f'{int(w*factorinterno)}%' if (w := v.get_width()) > 0 else f"{int(-v.get_width()*factorinterno)}" for v in c]
   ax.bar_label(c, labels=labels, label_type='edge', padding=0.3 )
ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=5)
ax.grid(b=True, color='green',
            linestyle='-.', linewidth=0.5,
            alpha=0.2)
ax.invert_yaxis()
ax.set_title("Barplot de las transacciones promedio por linea para un dia habil del mes de Junio del año 2022",   loc='left', fontname="Gotham Rounded", fontsize=20, pad=10, color="dimgrey")
plt.suptitle(
        "Transacciones promedio por Linea",
        fontsize=35,
        fontweight="bold", fontname="Gotham Rounded",
        x=0.123,
        y=0.975,
        ha="left")
for tick in ax.get_xticklabels():
   tick.set_fontname("Lato")
   tick.set_fontsize("15")
for tick in ax.get_yticklabels():
   tick.set_fontname("Lato")
   tick.set_fontsize("15")

ax.text(0.5, 0.5, 'NonLinear', transform=ax.transAxes,
            fontsize=90, color='grey', alpha=0.3,
            ha='center', va='center', rotation=30)

plt.annotate('Fuente: Datos SUBE-Santa Fe, Junio 2022', (0, 0), (0, -30), fontsize=9,
                 xycoords='axes fraction', textcoords='offset points', va='top')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
logo = plt.imread('Rawdata/LogoMCSF.jpg')
fig.figimage(logo, 1800, 800, zorder=3)
plt.show()

#POR ATRIBUTO
pptPar=dfhabil.groupby(['CODIGOCONTRATO', "DIA"]).size().to_frame('NumerodeTrans').reset_index()
thabiles = len(pptPar["DIA"].unique())
PromedioHabiles=pptPar.groupby(["CODIGOCONTRATO"]).apply(lambda NumerodeTrans: (NumerodeTrans.sum()/thabiles)).reset_index().sort_values(by = "NumerodeTrans", ascending=False)
PromedioHabiles=PromedioHabiles.drop(["DIA"], axis=1)
PromedioHabiles = PromedioHabiles.astype({'NumerodeTrans':'int'}) #Esto porque estaban como float


natributo = PromedioHabiles["NumerodeTrans"].sum()

factorinterno = 100 / natributo
fig, ax = plt.subplots(figsize=(16,9 ))
ax.barh(PromedioHabiles["CODIGOCONTRATO"].astype(str), PromedioHabiles["NumerodeTrans"], color=["#12bc8e"])
for c in ax.containers:
   labels = [f'{int(w*factorinterno)}%' if (w := v.get_width()) > 0 else f"{int(-v.get_width()*factorinterno)}" for v in c]
   ax.bar_label(c, labels=labels, label_type='edge', padding=0.3 )
ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=5)
ax.grid(b=True, color='green',
            linestyle='-.', linewidth=0.5,
            alpha=0.2)
ax.invert_yaxis()
ax.set_title("Barplot de las transacciones por atributo para un dia promedio del mes de Junio del año 2022",   loc='left', fontname="Gotham Rounded", fontsize=20, pad=10, color="dimgrey")
plt.suptitle(
        "Transacciones por Atributo",
        fontsize=35,
        fontweight="bold", fontname="Gotham Rounded",
        x=0.19,
        y=0.975,
        ha="left")
for tick in ax.get_xticklabels():
   tick.set_fontname("Lato")
   tick.set_fontsize("15")
for tick in ax.get_yticklabels():
   tick.set_fontname("Lato")
   tick.set_fontsize("15")

ax.text(0.5, 0.5, 'NonLinear', transform=ax.transAxes,
            fontsize=90, color='grey', alpha=0.3,
            ha='center', va='center', rotation=30)

plt.annotate('Fuente: Datos SUBE-Santa Fe, Junio 2022', (0, 0), (0, -30), fontsize=9,
                 xycoords='axes fraction', textcoords='offset points', va='top')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
logo = plt.imread('Rawdata/LogoMCSF.jpg')
fig.figimage(logo, 1800, 800, zorder=3)
plt.tight_layout()
plt.show()


#MATRIZ DE CALOR POR LINEA Y HORA
pptPar= dfhabil.groupby(["idlinea", "HORA"]).size().to_frame("NumerodeTrans").reset_index()
pptPar["NumerodeTrans"] = pptPar["NumerodeTrans"] / len(dfhabil["DIA"].unique())
ax = sn.heatmap(pptPar.pivot("idlinea", "HORA", "NumerodeTrans"), linewidth=0.5, vmin=pptPar["NumerodeTrans"].min(), vmax=pptPar["NumerodeTrans"].max(), cmap="Greens")
ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=5)
ax.grid(b=True, color='green',
            linestyle='-.', linewidth=0.5,
            alpha=0.2)
ax.set_title("Matriz de calor de transacciones por linea y hora",   loc='left', fontname="Gotham Rounded", fontsize=20, pad=10, color="dimgrey")
plt.suptitle(
        "Transacciones por Linea y Hora",
        fontsize=35,
        fontweight="bold", fontname="Gotham Rounded",
        x=0.122,
        y=0.975,
        ha="left")
for tick in ax.get_xticklabels():
   tick.set_fontname("Lato")
   tick.set_fontsize("15")
for tick in ax.get_yticklabels():
   tick.set_fontname("Lato")
   tick.set_fontsize("15")

ax.text(0.5, 0.5, 'NonLinear', transform=ax.transAxes,
            fontsize=90, color='grey', alpha=0.3,
            ha='center', va='center', rotation=30)

plt.annotate('Fuente: Datos SUBE-Santa Fe, Junio 2022', (0, 0), (0, -30), fontsize=9,
                 xycoords='axes fraction', textcoords='offset points', va='top')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
logo = plt.imread('Rawdata/LogoMCSF.jpg')
fig.figimage(logo, 1800, 800, zorder=3)
plt.show()


#MAPA DE CALOR
df=df[df['DIA']==15]
lineas= df["idlinea"].unique()
lineas.astype(numpy.int64)
lineas = [4,5]
listadelista = []
style1 = {'fillColor': '#5b5b5f', 'color': '#5b5b5f'}

class HeatMapWithTimeAdditional(Layer):
   _template = Template("""
       {% macro script(this, kwargs) %}
           var {{this.get_name()}} = new TDHeatmap({{ this.data }},
               {heatmapOptions: {
                   radius: {{this.radius}},
                   minOpacity: {{this.min_opacity}},
                   maxOpacity: {{this.max_opacity}},
                   scaleRadius: {{this.scale_radius}},
                   useLocalExtrema: {{this.use_local_extrema}},
                   defaultWeight: 1,
                   {% if this.gradient %}gradient: {{ this.gradient }}{% endif %}
               }
           }).addTo({{ this._parent.get_name() }});
       {% endmacro %}
   """)

   def __init__(self, data, name=None, radius=15,
                min_opacity=0, max_opacity=0.6,
                scale_radius=False, gradient=None, use_local_extrema=False,
                overlay=True, control=True, show=True):
      super(HeatMapWithTimeAdditional, self).__init__(
         name=name, overlay=overlay, control=control, show=show
      )
      self._name = 'HeatMap'
      self.data = data

      # Heatmap settings.
      self.radius = radius
      self.min_opacity = min_opacity
      self.max_opacity = max_opacity
      self.scale_radius = 'true' if scale_radius else 'false'
      self.use_local_extrema = 'true' if use_local_extrema else 'false'
      self.gradient = gradient
#Esta clase hay que agregarla porque la libreria folium por defecto no permite hacer heatmapwithtime con layers.

mapa = folium.Map([-31.63238910, -60.69945910], tiles="OpenStreetMap", zoom_start=14 )

for indice in lineas: #para cada linea
   h = int(indice)
   auxiliar= df[df["idlinea"]==indice] #Filtrar por linea

   dfHeatTime=auxiliar.groupby(['latitude','longitude','HORA']).size().to_frame('NumerodeTrans').reset_index().sort_values(by=['HORA'],ascending=True) # Agrupo por lat long y horas.
   dfHeatTime['NumerodeTrans']=dfHeatTime['NumerodeTrans']/max(dfHeatTime['NumerodeTrans'])*10
   dfHeatTime['HORA'] = dfHeatTime['HORA'].sort_values(ascending=True) # Ordeno de manera ascedente, lo pide la funcion HeatMap.

   for i in range(24):#Esto es para incorporar un dato para cada hora, sino luego faltan horas. El NTRX es infinitesimal.
      dfAux=pd.DataFrame({'latitude':[-31.63238910],'longitude':[-60.69945910],'HORA':i,'NumerodeTrans':[0.000001]})
      dfHeatTime=pd.concat([dfHeatTime,dfAux])

   data = [] #Aca genero una lista vacia, que sera la lista de listas, de una dimension = 24 (porque hay 24 hs)
   for _, d in dfHeatTime.groupby('HORA'):
      data.append([[row['latitude'], row['longitude'], row['NumerodeTrans']] for _, row in d.iterrows()])
   listadelista += [data]


for i in range(2):
    a = gpd.read_file("Preprocessdata/Recorridos/Recorridos2021/Recorrido_" + str(int(lineas[i])) + ".shp")  # Leer el archivo con el recorrido
    fg = folium.FeatureGroup("LINEA" + " " + str(int(lineas[i])))
    if i == 0:
        plugins.HeatMapWithTime(data=listadelista[i], index=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], auto_play=True, display_index=True, name= str(int(lineas[i]))).add_to(mapa)
        folium.GeoJson(data=a["geometry"], style_function=lambda x:style1).add_to(fg)
    else:
        HeatMapWithTimeAdditional(data=listadelista[i], name= str(int(lineas[i]))).add_to(fg)
        folium.GeoJson(data=a["geometry"], style_function=lambda x: style1).add_to(fg)
    fg.add_to(mapa)
folium.LayerControl().add_to(mapa)
loc = "Mapa de calor lineas 4 y 5 noviembre 2021"
title_html = '''
                <h3 align="center" style="font-size:16px"><b>{}</b></h3>
               '''.format(loc)
mapa.get_root().html.add_child(folium.Element(title_html))
mapa.save("Output/MapadeCalorNov2021"+".html") # Guardo.
print ("Mapa Generado")


#TablaCruzada casos n=x
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

TablaCruzada= pd.DataFrame()
dias = dfhabil["DIA"].unique()
cases = [1,2,3,4]
TablaCruzada["NroDeTrx"] = None

for c in range(len(cases)):
    a = 0
    for d in dias:
        aux = dfhabil[dfhabil["DIA"] == d]
        a += len(casosn(aux,cases[c]))
    TablaCruzada.at[c,"NroDeTrx"] = a / len(dfhabil["DIA"].unique())
    TablaCruzada.at[c,"Caso"] = cases[c]

TablaCruzada["Porcentaje relativo"] = TablaCruzada["NroDeTrx"] / TablaCruzada["NroDeTrx"].sum()


#Tabla Cruzada Casos nx y atributo

TablaN = pd.DataFrame()
TablaN["Tipo de tarifa"] = df["CODIGOCONTRATO"].unique()
TablaN["N = 1"] = None
TablaN["N = 2"] = None
TablaN["N = 3"] = None
TablaN["N = 4"] = None

cases = [1,2,3,4]
#
for i in range(len(cases)):
    for w in range(len(df["CODIGOCONTRATO"].unique())):
        a = 0
        for d in dias:
            aux = dfhabil[dfhabil["DIA"] == d]
            auxiliar1 = casosn(aux, cases[i])
            auxiliar2 = auxiliar1[auxiliar1["CODIGOCONTRATO"] == TablaN.iloc[w, 0]]
            a += len(auxiliar2)
        TablaN.iloc[w, i+1] = a


