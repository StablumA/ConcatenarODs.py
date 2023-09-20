import pandas as pd
import geopandas as gpd



def carrilrapido(lineas,velrapida,direccionshp,nv):
    """Ingresar Lista de líneas, velocidad promedio del carril rápido, shp del carril, nodosvelocidades"""
    shpaux = gpd.read_file(direccionshp) # Cargo polígono de carril rápido
    nvaux = nv[nv["idlinea"].isin(lineas)] # Filtro las líneas que son afectadas por el carril rápido

    nvaux =  gpd.GeoDataFrame(nvaux,geometry=gpd.points_from_xy(nvaux["longitude"],nvaux["latitude"]),crs="EPSG:4326")
    nvrapido = gps.sjoin(nvaux,shpaux,how='inner') # Los nodos dentro del carril rápido
    vpromedioaux = nvrapido["Velocidad"].mean() # Calculo su velocidad promedio actual
    nvrapido["FV"] = nvrapido["Velocidad"] / vpromedioaux # Calculo factor
    nvrapido["Velocidad"] = nvrapido["FV"] * velrapida # Recalculo con la nueva velocidad promedio
    coloriginales = nvaux.columns.to_list()
    nvrapido = nvrapido[coloriginales] # Quito las columnas del shp
    nvcomplemento = nvaux[~nvaux.index.isin(nvrapido.index)] # Me quedo con los nodos que no estuvieron dentro del shp
    nvsalida = pd.concat([nvrapido,nvcomplemento]) # Uno ambos nodos
    nvsalida = nvsalida.drop_duplicates(["IDtramo","HORA"],keep="last").sort_index()
    return nvsalida
