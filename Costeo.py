import pandas as pd
import numpy as np
from Turnero import calcularturnos


def costos_km(dfc):
    aux1 = dfc[dfc['Unidad'] == '$/N']  # Costos por unidaad
    aux2 = dfc[dfc['Unidad'] == 'N/km']  # Consumo por km
    aux3 = dfc[dfc['Unidad'] == 'km/N']  # Vida útil
    aux4 = pd.merge(aux1, aux2, how='inner', on=['Concepto'])
    aux5 = pd.merge(aux1, aux3, how='inner', on=['Concepto'])
    aux4['Valor'] = aux4.apply(lambda x: x['Valor_x'] * x['Valor_y'], axis=1)
    aux5['Valor'] = aux5.apply(lambda x: x['Valor_x'] / x['Valor_y'], axis=1)
    aux6 = pd.concat([aux4, aux5])
    aux6['Unidad'] = '$/km'
    aux6['Tipo'] = 'Costo'
    aux6 = aux6[['Tipo', 'Concepto', 'Valor', 'Unidad']]
    return aux6


def costos_vmes(dfc):
    aux1 = dfc[dfc['Unidad'] == '$/Vmes']  # Costo por vehiculo por mes
    aux2 = dfc[dfc['Unidad'] == '$/V']  # Costo por vehículo
    aux3 = dfc[dfc['Unidad'] == 'mes']  # Vida útil en meses
    aux01 = dfc[dfc['Unidad'] == '$/N']  # No hay casos pero podría
    aux02 = dfc[dfc['Unidad'] == 'N/Vmes']  # No hay casos pero podría
    aux4 = pd.merge(aux2, aux3, how='inner', on=['Concepto'])
    aux5 = pd.merge(aux01, aux02, how='inner', on=['Concepto'])
    aux4['Valor'] = aux4.apply(lambda x: x['Valor_x'] / x['Valor_y'], axis=1)
    aux5['Valor'] = aux5.apply(lambda x: x['Valor_x'] * x['Valor_y'], axis=1)
    aux6 = pd.concat([aux1, aux4, aux5])
    aux6['Unidad'] = '$/Vmes'
    aux6['Tipo'] = 'Costo'
    aux6 = aux6[['Tipo', 'Concepto', 'Valor', 'Unidad']]
    return aux6


def costossalariales_vmes(dfC, dfparam):
    """Calculo de costos salariales excluyendo a los choferes, con Contribuciones Sociales"""
    aux1 = dfC[(dfC['Tipo'] == 'Salario') & (dfC['Unidad'] == '$/Pmes')]
    unidad = ['P/V', 'P/V', 'P/V']
    Tipo = ['Promedio de personal por vehículo', 'Promedio de personal por vehículo',
            'Promedio de personal por vehículo']
    Concepto = ['Mantenimiento', 'Tráfico', 'Administrativo']
    N_Mantenimiento_V = dfparam['Personal de mantenimiento por coche'].iloc[
        0]  # Personal promedio de mantenimiento por coche (obtenido de matriz costos 2022)
    N_Trafico_V = dfparam['Personal de tráfico por coche'].iloc[
        0]  # Personal promedio de tráfico por coche (obtenido de matriz costos 2022)
    N_Administrativo_V = dfparam['Personal de administración por coche'].iloc[
        0]  # Personal promedio de administración por coche (obtenido de matriz costos 2022)
    CS = dfparam['Contribuciones Sociales'].iloc[0]
    Valores = [N_Mantenimiento_V, N_Trafico_V, N_Administrativo_V]
    aux2 = pd.DataFrame(data={'Tipo': Tipo, 'Concepto': Concepto, 'Valor': Valores, 'Unidad': unidad})
    aux3 = pd.merge(aux1, aux2, how='inner', on=['Concepto'])
    aux3['Valor'] = aux3.apply(lambda x: x['Valor_x'] * x['Valor_y'] * (1 + CS), axis=1)
    aux3['Unidad'] = '$/Vmes'
    aux3['Tipo'] = 'Salario (con CS)'
    aux3 = aux3[['Tipo', 'Concepto', 'Valor', 'Unidad']]
    return aux3


def CostosEstructura_Vmes(dfC):  # estos van a sumar al costo final total
    """Leo los inputs de costos por estructura empresarial si los hubiera"""
    # Filtrar los elementos de la lista que contengan la palabra clave (sin distinguir mayúsculas y minúsculas)
    Inmuebles = dfC[(dfC['Concepto'].str.contains('inmuebles', case=False))]
    Vigilancia = dfC[dfC['Concepto'].str.contains('vigilancia', case=False)]
    aux3 = pd.concat([Inmuebles, Vigilancia])
    aux3['Unidad'] = '$/Vmes'
    aux3['Tipo'] = 'Costo Amortización'
    aux3 = aux3[['Tipo', 'Concepto', 'Valor', 'Unidad']]
    return aux3


def costeototal(direccioncostos, direccionparametros, dfIPK, direcciondestinodf, tipo, TARIFA):
    """ Ingresar direccion df con Costos , direccion df con parametros, df con datos por turno, direccion destino,
      'mes' o 'dia" , tarifa """
    data = pd.read_excel(direccioncostos) # Leo costos y consumos
    dfIPK = dfIPK  # cargo df con km recorridos, num vehículos e IPK
    tiposturnos = [8]
    dfChofer = calcularturnos(dfIPK, tiposturnos)
    dfparam = pd.read_excel(direccionparametros,header=None).T
    dfparam.columns = dfparam.iloc[0]
    dfparam.drop(0,axis=0,inplace=True)
    lineas = dfIPK['Linea'].unique()
    lineasexcluidas = dfparam['Líneas Excluídas'].to_list()
    lineas = [linea for linea in lineas if (linea not in lineasexcluidas)]
    nhab = dfparam['Hábiles'].iloc[0]
    sab = dfparam['Sábados'].iloc[0]
    domfer = dfparam['Domingos y Feriados'].iloc[0]
    # Calculo número de días regulares
    ndias_mes = nhab + 0.6 * domfer + 0.8 * sab
    CS = dfparam['Contribuciones Sociales'].iloc[0]  # Contribuciones sociales
    VREPUESTO = dfparam['Vehículos Repuesto'].iloc[
        0]  # Cantidad de vehículos de 'repuesto' x línea además de los necesarios para realizar los recorridos
    FEGasoil = dfparam["Factor emisión Gasoil"].iloc[0]
    N_Gasoil_km = data[(data["Tipo"]=="Consumo") & (data["Concepto"]=="Gasoil")]["Valor"].iloc[0] #Consumo de Gasoil por km
    E_km = FEGasoil * N_Gasoil_km   # Emisión de kg eqCO2 / km
    E_km = {"Tipo":"Emisión","Concepto":"GEI","Valor":E_km,"Unidad":"kg.eqCO2/km"}

    N_Mantenimiento_V = dfparam['Personal de mantenimiento por coche'].iloc[
        0]  # Personal promedio de mantenimiento por coche (obtenido de matriz costos 2022)
    N_Trafico_V = dfparam['Personal de tráfico por coche'].iloc[
        0]  # Personal promedio de tráfico por coche (obtenido de matriz costos 2022)
    N_Administrativo_V = dfparam['Personal de administración por coche'].iloc[
        0]  # Personal promedio de administración por coche (obtenido de matriz costos 2022)
    capacidad = 40 # Capacidad máxima de cada colectivo.
    # Creo hoja a cargar datos de salida
    writer = pd.ExcelWriter(direcciondestinodf, engine='xlsxwriter')


    # Comienzo a preprocesar y agrupar los costos
    costosconsiderados = pd.DataFrame()
    dfcostos_km = costos_km(data)  # recordar que lo multiplicado por km ya incluye multiplicación por Num Vehiculos
    costosconsiderados = pd.concat([costosconsiderados,dfcostos_km])
    dfcostos_Vmes = costos_vmes(data)
    costosconsiderados = pd.concat([costosconsiderados,dfcostos_Vmes])
    dfcostossalariales_Vmes = costossalariales_vmes(data, dfparam)
    costosconsiderados = pd.concat([costosconsiderados,dfcostossalariales_Vmes])

    def Costeo(linea):
        dfaux = dfIPK.loc[dfIPK['Linea'] == linea]  # Filtro la línea a trabajar
        dfaux = dfaux[['Linea', 'TRXTurno', 'NroVehiculos', 'DistanciaRecorrida','IPK','TiempoDeViaje']]  # filtro columnas de interés
        dfaux2 = dfChofer[dfChofer['Linea'] == linea]  # Leo los choferes para esta línea
        V = dfaux['NroVehiculos'].max() + VREPUESTO  # Número vehiculos máximo utilizado más de repuesto
        km_dia = dfaux['DistanciaRecorrida'].sum()  # Distancia recorrida por día
        km_mes = km_dia * ndias_mes
        N_Choferes = dfaux2['NroVehiculos'].iloc[0]


        ocupacionprom_turno= (( dfaux["TRXTurno"] / dfaux["NroVehiculos"] ) * dfaux["TiempoDeViaje"]/60)
        paxs_dia = dfaux['TRXTurno'].sum()
        paxs_mes = paxs_dia * ndias_mes
        IPK = dfaux["IPK"].mean()
        IPKs.append(IPK)

        E_paxTurno = {"Concepto":"Huella de CO2 x pasajero","Valor":[(E_km["Valor"]*km_dia/(V*ocupacionprom_turno))*1000],"Unidad":"gr.eqCO2/pasajero","Linea":linea}
        E_pax = {"Concepto":"Huella de CO2 x pasajero","Valor":[((E_km["Valor"]*km_dia/(V*ocupacionprom_turno))*1000).mean()],"Unidad":"gr.eqCO2/pasajero","Linea":linea}
        E_kmpaxTurno = {"Concepto":"Huella de CO2 x pasajero x km","Valor":(E_km["Valor"]/ocupacionprom_turno)*1000,"Unidad":"gr.eqCO2/pasajero.km","Linea":linea}
       ## Aún las por turno no las usé para volverlas a incluir en el IPK -> Coming Soon?
        E_kmpax = {"Concepto":"Huella de CO2 x pasajero x km","Valor":((E_km["Valor"]/ocupacionprom_turno).mean())*1000,"Unidad":"gr.eqCO2/pasajero.km","Linea":linea}
        dfE_kmpax = pd.DataFrame({"Concepto":"Huella de CO2 x pasajero x km","Valor":[((E_km["Valor"]/ocupacionprom_turno).mean())*1000],"Unidad":"gr.eqCO2/pasajero.km","Linea":linea})
        dfE_pax = pd.DataFrame(E_pax)
        E_dia = {"Concepto":"Huella de CO2","Valor":[E_km["Valor"]*km_dia],"Unidad":"kg.eqCO2/dia"}
        E_mes = {"Concepto":"Huella de CO2","Valor":[E_km["Valor"]*km_mes],"Unidad":"kg.eqCO2/mes"}
        Emisiones_Pasajero.append(dfE_pax)
        EPK.append(E_kmpax)

        P_Chofer_N = data[(data['Concepto'] == 'Chofer') & (data['Tipo'] == 'Salario')]['Valor'].iloc[0]
        N_Mantenimiento = N_Mantenimiento_V * V  # Personal promedio de mantenimiento por coche (obtenido de matriz costos 2022)
        N_Trafico = N_Trafico_V * V  # Personal promedio de tráfico por coche (obtenido de matriz costos 2022)
        N_Administrativo = N_Administrativo_V * V  # Personal promedio de administración por coche (obtenido de matriz costos 2022)

        ingreso_dia = TARIFA * paxs_dia
        ingreso_mes = TARIFA * paxs_mes

        CostoPersonalTotal_mes = (dfcostossalariales_Vmes['Valor'].sum() * V) + ((1 + CS) * N_Choferes * P_Chofer_N)

        Costo_KM = dfcostos_km['Valor'].sum()
        CostoKM_mes = Costo_KM * km_mes

        CostoCoches_Vmes = dfcostos_Vmes['Valor'].sum()
        CostoCoches_mes = V * CostoCoches_Vmes

        # CostoEstructura_Vmes  = dfcostosestructura_Vmes['Valor'].sum()
        # CostoEstructura_mes   = CostoEstructura_Vmes * V
        CostoTotal_LineaxMes = CostoPersonalTotal_mes + CostoKM_mes + CostoCoches_mes  # + CostoEstructura_mes

        Balance_mes = ingreso_mes - CostoTotal_LineaxMes
        Balance_dia = (Balance_mes / ndias_mes)
        TIPOCOSTO.append(
            {'Línea': linea, 'KM': CostoKM_mes, 'Coches': CostoCoches_mes, 'Salarios': CostoPersonalTotal_mes,
             'Total': CostoTotal_LineaxMes})
        Costo_Pasajero = CostoTotal_LineaxMes / paxs_mes
        LCosto_Pasajero.append(Costo_Pasajero)

        if (tipo == 'día') | (tipo == 'dia'):  # si elijo diario
            NPasajeros.append(paxs_dia)
            CostoPersonalTotal_dia = CostoPersonalTotal_mes / ndias_mes
            CostoKM_dia = Costo_KM * km_dia
            CostoCoches_Vdia = CostoCoches_Vmes / ndias_mes
            CostoCoches_dia = CostoCoches_Vdia * V
            CostoTotal_LineaxDia = CostoKM_dia + CostoCoches_dia + CostoPersonalTotal_dia
            SensibleKM = pd.DataFrame(
                {'Tipo': 'Costo', 'Concepto': 'SENSIBLE A LA CANTIDAD DE KM.', 'Valor': [dfcostos_km['Valor'].sum()],
                 'Unidad': '$/km'})
            SensibleV = pd.DataFrame({'Tipo': 'Costo', 'Concepto': 'SENSIBLE A LA CANTIDAD DE VEHÍCULOS.',
                                      'Valor': [dfcostos_Vmes['Valor'].sum() / ndias_mes], 'Unidad': '$/Vdia'})

            Emisiones.append(E_dia)
            dfE_dia = pd.DataFrame(E_dia)

            dfSUB = pd.DataFrame({'Concepto': ['SUBTOTAL DIARIO DEBIDO A PERSONAL CONTRATADO',
                                               'SUBTOTAL CANTIDAD DE KM POR DÍA',
                                               'SUBTOTAL DIARIO POR CANTIDAD DE COCHES'],
                                  'Valor': [CostoPersonalTotal_dia, CostoKM_dia, CostoCoches_dia],
                                  'Unidad': ['$', '$', '$']})

            dfIngreso = pd.DataFrame(
                {'Concepto': 'Ingresos estimados por ' + str(tipo), 'Valor': [ingreso_dia], 'Unidad': '$'})
            INGRESOS.append(dfIngreso)
            dictB = {'Concepto': ('BALANCE DIARIO LINEA ' + str(linea)), 'Valor': [Balance_dia],
                     'Unidad': ('$/' + str(tipo))}
            dfBalance = pd.DataFrame(dictB)
            BALANCES.append(dfBalance)

            Subsidio_L = {'Concepto': 'SUBSIDIO NECESARIO (sobre valor tarifa)',
                          'Valor': [max(0, ((-Balance_dia / paxs_dia) / TARIFA) * 100)], 'Unidad': '%'}
            Subsidio_L = pd.DataFrame(Subsidio_L)
            dfBalance = pd.concat([dfBalance, Subsidio_L])


            KMTOTAL.append(km_dia)
            nombres = ["Línea",'Número de coches (máximo por día)', 'Número de choferes', 'Personal de Mantenimiento (prom.)',
                       'Personal de Administración (prom.)', 'Personal de Tránsito (prom.)', 'Contribuciones Sociales',
                       'COSTOS SALARIALES DIARIOS', 'Kilómetros diarios recorridos',
                       'Cantidad Pasajeros por día (prom.)', 'Monto Tarifa']
            unidades = ['','V', 'P', 'P', 'P', 'P', '%', '$/' + str(tipo), 'km/' + str(tipo), ('Pasajeros/' + str(tipo)),
                        '$']
            Calculos = pd.DataFrame(data={'Concepto': nombres,
                                          'Valor': [linea,V, N_Choferes, N_Mantenimiento, N_Administrativo, N_Trafico,
                                                    CS * 100, CostoPersonalTotal_dia, km_dia, paxs_dia, TARIFA],
                                          'Unidad': unidades})
            DatosF1 = pd.DataFrame(
                {'Concepto': ['COSTO TOTAL DIARIO LINEA ' + str(linea), 'COSTO PROMEDIO POR PASAJERO', 'IPK', ''],
                 'Valor': [CostoTotal_LineaxDia, Costo_Pasajero, IPK, ''],
                 'Unidad': ['$', '$/Pasajero', '', '']})
            DatosF2 = pd.DataFrame(data={'Concepto': ['1. COSTOS SENSIBLES AL PERSONAL CONTRATADO %',
                                                      '2. COSTOS SENSIBLES A LA CANTIDAD DE COCHES %',
                                                      '3. COSTOS SENSIBLES A LA CANTIDAD DE KM. %', 'TOTAL %'],
                                         'Valor': [(CostoPersonalTotal_dia / CostoTotal_LineaxDia) * 100,
                                                   (CostoCoches_dia / CostoTotal_LineaxDia) * 100,
                                                   (CostoKM_dia / CostoTotal_LineaxDia) * 100,
                                                   (CostoTotal_LineaxDia / CostoTotal_LineaxDia) * 100],
                                         'Unidad': ['%', '%', '%', '%']})
            DatosF = pd.concat([DatosF1, DatosF2])
            DatosSalida = pd.concat([Calculos, SensibleKM, SensibleV,dfE_dia,dfE_pax,dfE_kmpax,dfSUB, dfIngreso, dfBalance, DatosF])
            DatosSalida = DatosSalida[['Concepto', 'Valor', 'Unidad']]
            DatosSalida.to_excel(writer, sheet_name='Línea ' + str(linea), index=True)
            DFSALIDAxL.append(DatosSalida)
            CostoTotal_Linea = CostoTotal_LineaxDia

        else:
            KMTOTAL.append(km_mes)
            NPasajeros.append(paxs_mes)
            SensibleKM = pd.DataFrame(
                {'Tipo': 'Costo', 'Concepto': 'SENSIBLE A LA CANTIDAD DE KM.', 'Valor': [dfcostos_km['Valor'].sum()],
                 'Unidad': '$/km'})
            SensibleV = pd.DataFrame({'Tipo': 'Costo', 'Concepto': 'SENSIBLE A LA CANTIDAD DE VEHÍCULOS.',
                                      'Valor': [dfcostos_Vmes['Valor'].sum()], 'Unidad': '$/Vmes'})
            dfE_mes=pd.DataFrame(E_mes)
            Emisiones.append(E_mes)

            dfSUB = pd.DataFrame({'Concepto': ['SUBTOTAL MENSUAL DEBIDO A PERSONAL CONTRATADO',
                                               'SUBTOTAL MENSUAL CANTIDAD DE KM POR MES',
                                               'SUBTOTAL MENSUAL POR CANTIDAD DE COCHES'],
                                  'Valor': [CostoPersonalTotal_mes, CostoKM_mes, CostoCoches_mes],
                                  'Unidad': ['$', '$', '$']})

            dfIngreso = pd.DataFrame(
                {'Concepto': 'Ingresos estimados por ' + str(tipo), 'Valor': [ingreso_mes], 'Unidad': '$'})
            INGRESOS.append(dfIngreso)
            dictB = {'Concepto': ('BALANCE MENSUAL LINEA ' + str(linea)), 'Valor': [Balance_mes],
                     'Unidad': ('$/' + str(tipo))}
            dfBalance = pd.DataFrame(dictB)
            BALANCES.append(dfBalance)
            Subsidio_L = {'Concepto': 'SUBSIDIO NECESARIO (sobre valor tarifa)',
                          'Valor': [max(0, ((-Balance_mes / paxs_mes) / TARIFA) * 100)], 'Unidad': '%'}
            Subsidio_L = pd.DataFrame(Subsidio_L)
            dfBalance = pd.concat([dfBalance, Subsidio_L])

            KMTOTAL.append(km_dia)
            nombres = ['Línea','Número de coches', 'Número de choferes', 'Personal de Mantenimiento (prom.)',
                       'Personal de Administración (prom.)', 'Personal de Tránsito (prom.)', 'Contribuciones Sociales',
                       'COSTOS SALARIALES MENSUALES', 'Kilómetros recorridos por mes',
                       'Cantidad Pasajeros por mes (prom.)', 'Monto Tarifa']
            unidades = ['','V', 'P', 'P', 'P', 'P', '%', '$/' + str(tipo), 'km/' + str(tipo), 'Pasajeros/' + str(tipo),
                        '$']
            Calculos = pd.DataFrame(data={'Concepto': nombres,
                                          'Valor': [linea,V, N_Choferes, N_Mantenimiento, N_Administrativo, N_Trafico,
                                                    CS * 100, CostoPersonalTotal_mes, km_mes, paxs_mes, TARIFA],
                                          'Unidad': unidades})
            DatosF1 = pd.DataFrame(
                {'Concepto': ['COSTO TOTAL MENSUAL LINEA ' + str(linea), 'COSTO PROMEDIO POR PASAJERO', 'IPK', ''],
                 'Valor': [CostoTotal_LineaxMes, Costo_Pasajero, IPK, ''],
                 'Unidad': ['$', '$/Pasajero', '', '']})
            DatosF2 = pd.DataFrame(data={'Concepto': ['1. COSTOS SENSIBLES AL PERSONAL CONTRATADO %',
                                                      '2. COSTOS SENSIBLES A LA CANTIDAD DE COCHES %',
                                                      '3. COSTOS SENSIBLES A LA CANTIDAD DE KM. %', 'TOTAL %'],
                                         'Valor': [(CostoPersonalTotal_mes / CostoTotal_LineaxMes) * 100,
                                                   (CostoCoches_mes / CostoTotal_LineaxMes) * 100,
                                                   (CostoKM_mes / CostoTotal_LineaxMes) * 100,
                                                   (CostoTotal_LineaxMes / CostoTotal_LineaxMes) * 100],
                                         'Unidad': ['%', '%', '%', '%']})
            DatosF = pd.concat([DatosF1, DatosF2])
            DatosSalida = pd.concat([Calculos, SensibleKM, SensibleV,dfE_mes,dfE_pax,dfE_kmpax, dfSUB, dfIngreso, dfBalance, DatosF])
            DatosSalida = DatosSalida[['Concepto', 'Valor', 'Unidad']].reset_index()
            DatosSalida.to_excel(writer, sheet_name='Línea ' + str(linea), index=False)
            DFSALIDAxL.append(DatosSalida)
            CostoTotal_Linea = CostoTotal_LineaxMes
        return CostoTotal_Linea

    # Listas con datos de interés para el reporte.
    DFSALIDAxL=[] # lista con los df de salidas de cada línea
    IPKs = [] # lista con IPK por línea
    EPK = [] # lista con emision por línea x pax x km
    Emisiones = [] # lista con emision de GEI en kg por por mes/dia de cada linea
    Emisiones_Pasajero = [] # lista con emision de GEI en kg por pasajero de cada linea
    NPasajeros = []  # lista que guarda cantidad de pasajeros por línea por mes/dia
    TIPOCOSTO = []  # voy a ir sumando los tipos de costos
    LCosto_Pasajero = []  # lista con costos de cada pasajero de la línea
    COSTOS = []  # lista con los costos totales de cada linea
    KMTOTAL = []  # mensual o diario según el horizonte
    INGRESOS = [] # lista con ingresos por línea por mes/dia
    BALANCES = [] # idem con balances
    # abro el writer
    writer = pd.ExcelWriter(direcciondestinodf, engine='xlsxwriter')

    for linea in lineas:  # Corro la función para cada linea
        COSTOS.append(Costeo(linea))  # la sumo a mi lista
    COSTOSTxL = pd.DataFrame(data={'Línea': lineas, 'Costo': COSTOS})
    COSTOTOTAL = sum(COSTOS)  # Armo el costo total de la política

    # Diferencio los tipos de costos:
    TIPOCOSTOaux = pd.DataFrame(TIPOCOSTO)
    TIPOCOSTOaux = TIPOCOSTOaux[['Salarios', 'Coches', 'KM']]
    TIPOCOSTOaux = TIPOCOSTOaux.apply(lambda x: x.sum(), axis=0)
    TIPOCOSTO = pd.DataFrame({'Concepto': ['Costos debidos a los salarios', 'Costos debidos a la cantidad de coches',
                                           'Costos debidos a la cantidad de KM'], 'Valor': TIPOCOSTOaux,
                              'Unidad': ['$/' + str(tipo), '$/' + str(tipo), '$/' + str(tipo)]})

    costosconsiderados.to_excel(writer, sheet_name='Costos Preprocesados', index=False)
    # Armo costos por pasajero, y costo de la línea más cara:
    COSTOPROMxL = np.mean(COSTOS)
    COSTOPROMPASAJEROS = np.mean(LCosto_Pasajero)
    LineaCara = lineas[COSTOS.index(max(COSTOS))]
    COSTOMAXPASAJERO = max(LCosto_Pasajero)
    LineaCara_Pasajero = lineas[LCosto_Pasajero.index(COSTOMAXPASAJERO)]

    xPASAJERO = pd.DataFrame(
        {'Línea': lineas, 'Cantidad de pasajeros': NPasajeros, 'Costo por pasajero': LCosto_Pasajero})
    LCosto_Pasajero2 = LCosto_Pasajero.copy()  # copio la lista
    LCosto_Pasajero2.remove(max(LCosto_Pasajero))  # le saco el mas caro x pasajero
    COSTOPROMPASAJEROSsinmax = np.mean(LCosto_Pasajero2)
    NPASAJEROST = sum(NPasajeros)

    # CALCULO INGRESOS Y BALANCE
    INGRESOS = pd.concat(INGRESOS)
    INGRESOTOTAL = INGRESOS['Valor'].sum()
    BALANCES = pd.concat(BALANCES)
    BALANCETOTAL = BALANCES['Valor'].sum()

    # Calculo subsidio necesario para lograr balance 0
    SUBSIDIO = max((-BALANCETOTAL / NPASAJEROST), 0)  # Subsidio por cada pasaje
    PorcentajeSUBSIDIO = (SUBSIDIO / TARIFA) * 100  # Porcentaje del pasaje subsidiado

    # Calculo Huella de carbono total

    dfHuella = pd.DataFrame(EPK)
    dfHuella["Linea"] = lineas
    dfHuella = dfHuella.set_index("Linea").sort_values(by=["Valor"]).T
    dfHuella.to_excel(writer, sheet_name='Impacto Ambiental', index=False)

    IPKprom = np.mean(IPKs)

    EMISION = {"Concepto":"Huella de CO2 total","Valor":sum([dic["Valor"][0] for dic in Emisiones])/1000,"Unidad":("ton.eqCO2/"+str(tipo))}
    EMISION_PASAJEROprom = {"Concepto":"Huella de CO2 por pasajero promedio","Valor":np.mean([dic["Valor"][0] for dic in Emisiones_Pasajero]),"Unidad":("gr.eqCO2/pasajero")}
    EMISION_PASAJEROmax =  {"Concepto":"Huella de CO2 por pasajero máxima","Valor":max([dic["Valor"][0] for dic in Emisiones_Pasajero]),"Unidad":"gr.eqCO2/pasajero"}
    dfEPK= pd.DataFrame(EPK)
    E_PKprom={"Concepto":"Huella de CO2 x pasajero x km","Valor":[(dfEPK["Valor"].mean())],"Unidad":"gr.eqCO2/pasajero.km"}
    dfE_PKprom=pd.DataFrame(E_PKprom)
    EPKMmax = dfEPK[dfEPK["Valor"] == dfEPK["Valor"].max()].reset_index()


    # Genero sheet con comparación entre líneas
    conceptos,unidades= [c for c in DFSALIDAxL[0]["Concepto"]],[c for c in DFSALIDAxL[0]["Unidad"]]
    DFSALIDAxL1 = [df[["Valor"]].reset_index() for df in DFSALIDAxL]
    def ordenarporcosto(df):
        return df.at[21,"Valor"]
    DFSALIDAxL1=sorted(DFSALIDAxL1,key= ordenarporcosto,reverse=True)
    DFSALIDAxL2 = pd.concat(DFSALIDAxL1,axis=1)
    DFSALIDAxL2["Concepto"] = conceptos
    DFSALIDAxL2["Unidad"] = unidades
    DFSALIDAxL2 = DFSALIDAxL2[["Concepto","Valor","Unidad"]]
    DFSALIDAxL2.to_excel(writer,sheet_name="Comparación líneas",index=False)


    KMTOTAL = sum(KMTOTAL)  # Incluye el N de Vehículos por cómo está dado el dato de km recorridos.
    COSTOPROM_VKM = COSTOTOTAL / KMTOTAL
    NDatosReporte = ['Costo total' + ' por ' + str(tipo), 'Tarifa', 'Ingreso total', 'Balance total política',
                     'Subsidio necesario de la tarifa', 'Costo promedio por línea', 'Cantidad estimada de pasajeros',
                     'Costo por pasajero promedio','IPK promedio' ,EMISION["Concepto"],EMISION_PASAJEROprom['Concepto'],dfE_PKprom["Concepto"][0],EMISION_PASAJEROmax["Concepto"], 'Costo promedio por km de cada coche', 'Línea más cara',
                     'Costo línea más cara', 'Línea más cara x pasajero', 'Mayor huella de CO2 x pasajero x km ', "Línea con mayor emisión x pax x km",
                     'Costo x pasajero (linea más cara x pasajero)',
                     'Costo x pasajero promedio (sin línea más cara x pasajero)']
    DatosReporte = [COSTOTOTAL, TARIFA, INGRESOTOTAL, BALANCETOTAL, PorcentajeSUBSIDIO, COSTOPROMxL, NPASAJEROST,
                    COSTOPROMPASAJEROS,IPKprom,EMISION["Valor"],EMISION_PASAJEROprom['Valor'],dfE_PKprom["Valor"][0],EMISION_PASAJEROmax["Valor"],
                    COSTOPROM_VKM, LineaCara, max(COSTOS), LineaCara_Pasajero,EPKMmax["Valor"][0],EPKMmax["Linea"][0], COSTOMAXPASAJERO,
                    COSTOPROMPASAJEROSsinmax]

    UDSReporte = ['$/' + str(tipo), '$/Pasajero', '$/' + str(tipo), '$/' + str(tipo), '%', '$/' + str(tipo),
                  'Pasajeros', '$/Pasajero','Pasajero/km',EMISION["Unidad"],EMISION_PASAJEROprom['Unidad'],
                  dfE_PKprom["Unidad"][0], EMISION_PASAJEROmax["Unidad"],'$/V.km', '', '$/' + str(tipo), '',EPKMmax["Unidad"][0],"", '$/Pasajero', '$/Pasajero']
    dfReporte = pd.DataFrame({'Concepto': NDatosReporte, 'Valor': DatosReporte, 'Unidad': UDSReporte})
    dfReporte = pd.concat([TIPOCOSTO, dfReporte])
    dfReporte.to_excel(writer, sheet_name='Reporte', index=False)
    writer.close()

    return COSTOTOTAL, COSTOSTxL

