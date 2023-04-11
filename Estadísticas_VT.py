import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import geopandas as gpd
import folium
import openpyxl
from folium import plugins
from folium.map import Layer
from jinja2 import Template
df = pd.read_csv("10-TX-SData.csv",sep=',',parse_dates=['FECHATRX'],dayfirst=True)
df['FECHATRX']=pd.to_datetime(df['FECHATRX'])

# elimino primer columna
df=df.drop(df.columns[[0]],axis=1)

df['Hora']=df['FECHATRX'].dt.hour
# filtro por tipo de día
df['Dia']=df['FECHATRX'].dt.day
df['nDia']=df['FECHATRX'].dt.day_name()
#saco los sabados y domingos
dfhabil= df[df['nDia'].isin(['Monday','Tuesday','Wednesday','Thursday','Friday'])]
# saco los feriados
feriados = [7,10]
dfhabil= dfhabil[dfhabil['Dia']!= (i for i in feriados)]

## Veo num. de transacciones por  'Linea' en el total de días habiles
PromedioHabiles= dfhabil.groupby(['Ramal']).size().to_frame('NumerodeTrans').reset_index()
Thabiles= len(dfhabil['Dia'].unique())
# obtengo el promedio por día x linea

PromedioHabiles['NTransPromedio']= PromedioHabiles['NumerodeTrans'].apply(lambda NTrans: (NTrans/Thabiles))
#PromedioHabiles= PromedioHabiles.astype({"NumerodeTrans:"int"})
print(PromedioHabiles)

total = PromedioHabiles["NTransPromedio"].sum()
factor = 100/total

PromedioHabiles['FracUso']= PromedioHabiles["NTransPromedio"].apply(lambda cant: (cant*factor))
PromedioHabiles= PromedioHabiles.astype({'FracUso':'int'})
PromedioHabiles=PromedioHabiles.sort_values(by ="NTransPromedio", ascending=False).reset_index()

## realizo un gráfico
fig, ax = plt.subplots(figsize=(8,4))
ax.barh(PromedioHabiles['Ramal'].astype('str'), PromedioHabiles["NTransPromedio"], color=["#12bc8e"])

rects = ax.patches
# For each bar: Place a label
for i, rect in enumerate(rects):
    x_value = rect.get_width() / 2
    y_value = rect.get_y() + rect.get_height() / 2
    space = 5
    ha = 'left'
    # If value of bar is negative: Place label left of bar
    if x_value < 0:
        space *= -1
    label = "{:.0f}%".format(PromedioHabiles.at[i,'FracUso'])
    plt.annotate(label,                      # Use `label` as label
        (x_value, y_value),         # Place label at end of the bar
        xytext=(space, 0),          # Horizontally shift label by `space`
        textcoords="offset points", # Interpret `xytext` as offset in points
        va='center',                # Vertically center label
        ha=ha, fontsize=12)


ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=5) # dist de las etiquetas
ax.grid(color='green',linestyle='-.',linewidth=0.5,alpha=0.15)
ax.invert_yaxis()
ax.set_title("Barplot de las transacciones promedio en un día hábil, por línea.",
             loc='center', fontname='Arial',fontsize='12', pad=5,color='black')
plt.suptitle("Transacciones promedio por línea",
              fontname='Arial',fontsize='35',color='black', fontweight='bold',y=0.990)
ax.text(0.5, 0.5, 'NonLinear', transform=ax.transAxes,
            fontsize=80, color='gray', alpha=0.25,
            ha='center', va='center', rotation=25)
plt.annotate('Fuente: Datos Venado Tuerto, Octubre 2022', (0, 0), (0, -25), fontsize=9,
                 xycoords='axes fraction', textcoords='offset points', va='top')


plt.show()
#####

#### MAPA DE CALOR
pptPar= dfhabil.groupby(['Ramal','Hora']).size().to_frame('NumerodeTrans').reset_index()
pptPar['NumerodeTrans'] = pptPar['NumerodeTrans'] / len(df['Dia'].unique())
ax= sn.heatmap(pptPar.pivot(index="Ramal",columns="Hora",values="NumerodeTrans")
,linewidth=1,vmin=(pptPar['NumerodeTrans']).min(),vmax=(pptPar['NumerodeTrans']).max()
,cmap='Greens')
ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=5)
ax.grid( color='green',
            linestyle='-.', linewidth=0.5,
            alpha=0.2)
ax.set_title("Matriz de calor de transacciones por línea y hora", loc='center',fontsize=20,pad=10,color='darkblue')
plt.show()
####








###
# Cuantas personas toman 1 2 3 4 colectivos por día
def casosn(dataframe,n):
    """ Filtra numero de transacciones por día, por usuario"""
    dfaux1= dataframe['Tarjeta'].value_counts().reset_index() # cant viajes por mes
    dfaux1.columns=['ID','Cantidad'] # las renombro
    dfaux2= dfaux1['Cantidad'].value_counts().reset_index()
    dfaux2.columns=['Viajes por día','Cantidad']
    dfaux1=dfaux1[dfaux1['Cantidad']==n]
    salida= dataframe[dataframe['Tarjeta'].isin(dfaux1['ID'])]
    salida= salida.sort_values(by=['Tarjeta','FECHATRX'])
    salida= salida.reset_index()
    return salida

TablaCruzada= pd.DataFrame()
dias= dfhabil['Dia'].unique()
cases = [1,2,3,4]
TablaCruzada["NumerodeTrans"]=None

## note que hay algunas repetidas, debo quitarlas
for c in range(len(cases)):
    a = 0
    for d in dias:
        aux=dfhabil[dfhabil['Dia']==d]
        a += len(casosn(aux,cases[c]))
    TablaCruzada.at[c,'NumerodeTrans'] = a    # trans promedios por día de cda caso
    TablaCruzada.at[c,"Caso"] = cases[c]

TablaCruzada['PorcentajeRelativo'] = TablaCruzada['NumerodeTrans'] / TablaCruzada['NumerodeTrans'].sum()
TablaCruzada.to_excel('TablaCruzada.xlsx')













#########################
# Habrá que filtrar las que no tienen datos de Lat y Long
# pasemos las coordenadas de texto a números


for _,row in df.iterrows(): # para cada renglon
    if row['latitude']=='0' or row['longitude']=='0': # que tenga latitud o long
        df=df.drop(index=row.name) # quito el dato

df['latitude']=df['latitude'].str.replace(',','.').astype(float)
df['longitude']=df['longitude'].str.replace(',','.').astype(float)
df.reset_index().drop('index',axis=1) # reseteo indices y borro la columna que se me generó

#### Mapa de calor para un día hábil en particular
lineas=df['Ramal'].dropna().unique()
style1= style1 = {'fillColor' : '#5b5b5f', 'color': '#5b5b5f'}
listadelista = []
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
                 scale_radius=False, gradient=None,
                 use_local_extrema=False,overlay=True,
                 control=True, show=True):
        super(HeatMapWithTimeAdditional, self).__init__(
            name=name, overlay=overlay, control=control, show=show
        )
        self._name = 'HeatMap'
        self.data = data
        # heatmap settings
        self.radius = radius
        self.min_opacity = min_opacity
        self.max_opacity = max_opacity
        self.scale_radius = 'true' if scale_radius else 'false'
        self.use_local_extrema = 'true' if use_local_extrema else 'false'
        self.gradient = gradient

mapa = folium.Map([-33.74648, -61.9679], tiles="OpenStreetMap", zoom_start=14 )

for indice in lineas:
    auxiliar = df[df["Ramal"]==indice]
    dfHeatTime= auxiliar.groupby(["latitude","longitude","Hora"]).size().to_frame('NumerodeTrans').reset_index().sort_values(by=["Hora"], ascending=True)
    dfHeatTime['NumerodeTrans']=(dfHeatTime['NumerodeTrans']*10)/dfHeatTime['NumerodeTrans'].max()
    dfHeatTime['Hora']= dfHeatTime['Hora'].sort_values(ascending=True)

    for i in range(24): # incorporo un dato a cada hora para evitar error
        dfAux= pd.DataFrame({'latitude':[-33.74648], 'longitude':[-61.9679], 'Hora':i, 'NumerodeTrans':[0.00001]})
        dfHeatTime=pd.concat([dfHeatTime,dfAux])

    data = [] #Aca genero una lista vacia, que sera la lista de listas, de una dimension = 24 (porque hay 24 hs)
    for _, d in dfHeatTime.groupby('Hora'):
       data.append([[row['latitude'], row['longitude'], row['NumerodeTrans']] for _, row in d.iterrows()])
    listadelista += [data]

recorridos= gpd.read_file("Recorrido_TUP.shp")
recorridos= recorridos.drop(index=8)
recorridos= recorridos.drop(index=9)
#quité los recorridos de cementerio ida y vuelta

for _,reco in recorridos.iterrows():
    fg = folium.FeatureGroup('Línea'+' '+str(reco['Línea']))
    if _==0:
        plugins.HeatMapWithTime(data=listadelista[_], index=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23],
        auto_play=True, display_index=True, name= str(lineas[_])).add_to(mapa)
        folium.GeoJson(data=reco["geometry"], style_function=lambda x:style1).add_to(fg)
    else:
        HeatMapWithTimeAdditional(data=listadelista[_], name= str(lineas[_])).add_to(fg),
        folium.GeoJson(data=reco["geometry"], style_function=lambda x: style1).add_to(fg)
    fg.add_to(mapa)
folium.LayerControl().add_to(mapa)
loc = "Mapa de calor lineas 1 2 3 4 VT, Octubre 2022"
title_html = '''
                <h3 align="center" style="font-size:16px"><b>{}</b></h3>
               '''.format(loc)
mapa.get_root().html.add_child(folium.Element(title_html))
mapa.save("Mapacalor"+".html")
print ("Mapa Generado")
