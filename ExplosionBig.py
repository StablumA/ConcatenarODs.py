import pandas as pd
import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
import shapely
import numpy as np
import warnings
from shapely.geometry import Point
import pointpats
import time
from multiprocessing import Pool

warnings.filterwarnings('ignore')
def puntosenpoligono(args):
    poligono, n = args
    coor = pointpats.random.poisson(poligono, size=int(n))
    x = coor[0]
    y = coor[1]
    return x, y

def explotar(args):
    x, y,  zonificacion = args
    circulo = Point(x, y).buffer(0.002)
    interseccion = zonificacion[zonificacion.geometry.intersects(circulo)]
    if len(interseccion) == 0:
        print("ERROR")
        return x, y
    else:
        print("OK")
        interseccion["centrocirculo"] = Point(x, y)
        interseccion["distancia"] = interseccion.apply(lambda x: Point(x["long"], x["lat"]).distance(x["centrocirculo"]),
                                                       axis=1)
        interseccion["Importancia"] = interseccion["Personas"] / interseccion["distancia"]
        interseccion["p"] = interseccion["Importancia"] / interseccion["Importancia"].sum()
        poligonoelegido = np.random.choice(interseccion["geometry"], 1, p=interseccion["p"])
        nuevopunto = puntosenpoligono((poligonoelegido[0], 1))
        return nuevopunto

def parallel_explotar(args):
    return explotar(args)

if __name__ == '__main__':
    start = time.time()
    od = pd.read_csv("Venado/Output/OD_MaySepNovFINAL.csv")
    zona = gpd.read_file("Venado/SHP Files/Radios Censales/RC_VT_Cortado.shp", crs = "22185")
    zona = zona.to_crs("EPSG:4326")
    zona["long"] = zona.geometry.centroid.x
    zona["lat"] = zona.geometry.centroid.y
    zona["long"] = zona["long"] + 0.00001
    odaux = od.copy(deep=True)

    num_processes = 6  # Number of processes to run in parallel
    pool = Pool(num_processes)
    results = pool.map(parallel_explotar, zip(odaux["longitude_Origen"], odaux["latitude_Origen"], [zona]*len(odaux)))
    pool.close()
    pool.join()

    odaux[['longitude_Origen', 'latitude_Origen']] = pd.DataFrame(results, columns=['long', 'lat'])

    num_processes = 6  # Number of processes to run in parallel
    pool = Pool(num_processes)
    results = pool.map(parallel_explotar, zip(odaux["longitude_Destino"], odaux["latitude_Destino"], [zona] * len(odaux)))
    pool.close()
    pool.join()

    odaux[['longitude_Destino', 'latitude_Destino']] = pd.DataFrame(results, columns=['long', 'lat'])

    odaux.to_csv("Venado/Output/OD_MNS_EBIG.csv", index=True)
    end = time.time()
    print("Execution time:", ((end - start)/60), "minutos")
