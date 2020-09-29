# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

lon = []
lat = []
results = []
lpos = 0
longpos = 0
respos = 0

# Klobuchar parameters
A = 0                                                                    #azimuth
E = 90                                                                   #elevation angle
#Phi_u = lat
#Lambda_u = lon
#Phi_u = 40                                                              #geographic lat
#Lambda_u = -100                                                         #geographic longitude

#Time_GPS = 74700                                                         #GPS time
Time_GPS = 43200
Alpha = [3.82E-08,  1.49E-08, -1.79E-07, 0]
Beta = [1.43E+05, 0, -3.28E+05, 1.13E+05]

for i in range(-90, 91, 1):
    for j in range(-180, 181, 1):
        lat.insert(lpos ,i)
        lpos = lpos +1
        lon.insert(longpos, j)
        longpos = lpos +1
        # Earth-centered angle
        psi_i = (0.0137/ ((E/180) + 0.11)) - 0.022                        #in semicircles

        # Subionosphere latitude
        phi_i = i/180 + psi_i*np.cos (A*np.pi/180)/np.pi
        if phi_i > 0.416:
            phi_i = 0.416
        elif phi_i < -0.416:
            phi_i = -0.416

        # Subionosphere longitude
        lambdal = j/180 + psi_i*np.sin(A*np.pi/180)/np.cos(phi_i)

        # Geomagnetic latitude
        phi_m = phi_i + 0.064*np.cos((lambdal-1.617)*np.pi)/np.pi

        # Local time
        t = 4.32E+04*lambdal + Time_GPS
        if t > 86400:
            t = t - 86400
        elif t < 0:
            t = t + 86400

        # Slant factor
        F = 1.0 + 16*((0.53-E/180)**3)

        # Compute x
        Ai = Alpha[0] + Alpha[1] * phi_m + Alpha[2] * phi_m**2 + Alpha[3] * phi_m**3
        Bi = Beta[0] + Beta[1] * phi_m + Beta[2] * phi_m**2 + Beta[3] * phi_m**3

        x = 2*np.pi*(t-50400)/Bi

        # Compute Tiono for L1 frequency

        Tiono = F*(5E-09+Ai*(1-((x*x)/2)+((x**4)/24)))
        results.insert(respos, Tiono*299792458)
        respos = respos + 1
        #print("Tiono in seconds =", Tiono)
        #print("Tiono in meters =", Tiono*299792458)

#Geração de um data Frame para visualizar os dados em tabela
lista = [lat, lon, results]
df = DataFrame (lista).transpose()
df.columns = ['latitude','longitude','erro']
df.to_csv(r'dadosKlobuchar.csv')


fig = go.Figure()
fig.add_trace(go.Scatter(
    x=lon,
    y=lat,
    marker=dict(size=5, cmax=max(results), cmin=min(results), color=results, colorbar=dict(title="Tiono(meters)"),colorscale="Viridis"), mode="markers"))

fig.show()
