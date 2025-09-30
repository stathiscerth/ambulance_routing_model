# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 13:19:52 2023

@author: user
"""

import requests
import json
import math
import numpy as np
'''
mci1=[40.61578,22.74594 ]
mci2=[40.72410,23.12958]
'''
mci1 = [40.61578, 22.74594]
mci2 = [40.72410, 23.12958]

mci3 = [40.63265, 22.94120]
mci4 = [40.59920, 22.94969]


mci=[mci1,mci2,mci3,mci4]
h1=[40.57758,22.97037]
h2=[40.63050,23.04415]
h3=[40.63736,22.95888]
h4=[40.62965,22.96056]
h5=[40.63394,22.95572]
h6=[40.61283,22.96277]
h7=[40.67564,22.96160]


h=[h1,h2,h3,h4,h5,h6,h7]
#ax=numpy.random.standard_t(4,10000)
def distance_between(lon1,lat1,lon2,lat2):
    r = requests.get(f"http://router.project-osrm.org/route/v1/car/{lon1},{lat1};{lon2},{lat2}?overview=false""")
    #print(r.content)
    routes = json.loads(r.content)
    route_1 = routes.get("routes")[0]
    driving_mins=math.ceil(route_1['duration']/60) 
    distances=route_1['distance']/1000
    #print(distances)
    if distances>4.13:
        driving_m=math.ceil(2.46+0.596*distances)
    else:
        driving_m=math.ceil(2.42*math.sqrt(distances))
    return driving_m

def distance_km(lon1,lat1,lon2,lat2):
    r = requests.get(f"http://router.project-osrm.org/route/v1/car/{lon1},{lat1};{lon2},{lat2}?overview=false""")
    #print(r.content)
    routes = json.loads(r.content)
    route_1 = routes.get("routes")[0]
    driving_mins=math.ceil(route_1['duration']/60) 
    distances=route_1['distance']/1000
    #print(distances)
    
    return distances


def random_distance(lon1,lat1,lon2,lat2):
    r = requests.get(f"http://router.project-osrm.org/route/v1/car/{lon1},{lat1};{lon2},{lat2}?overview=false""")
    #print(r.content)
    routes = json.loads(r.content)
    route_1 = routes.get("routes")[0]
    #driving_mins=math.ceil(route_1['duration']/60) 
    distances=route_1['distance']/1000
    #print(distances)
    if distances>4.13:
        distance_m=math.ceil(2.46+0.596*distances)
    else:
        distance_m=math.ceil(2.42*math.sqrt(distances))
    c_d=math.sqrt(0.349+0.0006*distance_m+0.0388*distance_m*distance_m)/distance_m
    rand_time=distance_m*math.exp(c_d*np.random.standard_t(4,None))
    rand_time_c=math.ceil(rand_time)
    return rand_time_c
