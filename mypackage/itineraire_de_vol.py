from xml.dom.expatbuilder import parseString
import cartopy.crs as ccrs
from matplotlib import markers
import matplotlib.pyplot as plt
import csv
import pandas as pd
import math
import flagpy as fp
from PIL import Image


def distance_plot(dep, ariv):
    depart = []
    arrive = []

    f2 = open('mypackage\\airports.txt', encoding="utf8")
    airports = csv.reader(f2)

    for airport in airports:

        if airport[4] == dep:
            depart = [float(airport[6]), float(airport[7])]
        if airport[4] == ariv:
            arrive = [float(airport[6]), float(airport[7])]

    R = 6373.0

    lat1 = math.radians(depart[0])
    lon1 = math.radians(depart[1])
    lat2 = math.radians(arrive[0])
    lon2 = math.radians(arrive[1])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = math.sin(dlat / 2)**2 + math.cos(lat1) * \
        math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c

    fig = plt.figure(facecolor='#b0b29e')
    ax = plt.axes(projection=ccrs.Robinson())
    ax.stock_img()
    fig.set_size_inches(10, 5.5)

    plt.plot([depart[1], arrive[1]], [depart[0], arrive[0]],
             color='blue', linewidth=0.5, marker='*',
             transform=ccrs.Geodetic(),
             )

    plt.text(depart[1], depart[0], dep,
             horizontalalignment='right',
             transform=ccrs.Geodetic())

    plt.text(arrive[1], arrive[0], ariv,
             horizontalalignment='left',
             transform=ccrs.Geodetic())

    plt.savefig('static/resultats/route.png', dpi=1200)

    return math.ceil(distance)


def direct(dep, ariv):

    f = open('mypackage\\routes.txt')
    routescsv = csv.reader(f)

    vol = ''
    nb_airlines = 0

    for ligne in routescsv:

        if ligne[2] == dep:
            if ligne[4] == ariv:
                vol = "C'est un vol Direct"
                nb_airlines = nb_airlines+1

    if vol == '':
        vol = "Ce n'est pas un vol Direct "

    print(vol)
    print(nb_airlines)

    return vol, nb_airlines


def drapeau(dep, ariv):

    f2 = open('mypackage\\airports.txt', encoding="utf8")
    airports = csv.reader(f2)

    for airport in airports:
        if airport[4] == dep:

            airport_dep = airport[1]
            pays_dep = airport[3]
            Depart = fp.get_flag_img(pays_dep)
            Depart.save('static/resultats/depart.png')

        if airport[4] == ariv:

            airport_arr = airport[1]
            pays_arr = airport[3]
            Arrive = fp.get_flag_img(pays_arr)
            Arrive.save('static/resultats/arrive.png')

    print(pays_dep, airport_dep, pays_arr, airport_arr)

    return pays_dep, airport_dep, pays_arr, airport_arr


def company_plot(airline):

    MEA = []
    ar = []
    dep = []
    airo = []
    varhi = []
    new_airline_liste = []

    f = open('mypackage\\routes.txt')
    routescsv = csv.reader(f)

    f2 = open('mypackage\\airports.txt', encoding="utf8")
    airports = csv.reader(f2)

    f3 = open('mypackage\\airlines.txt', encoding="utf8")
    airlines = csv.reader(f3)

    def list_split(listA, n):
        for x in range(0, len(listA), n):
            every_chunk = listA[x: n+x]

            if len(every_chunk) < n:
                every_chunk = every_chunk + \
                    [None for y in range(n-len(every_chunk))]
            yield every_chunk

    for company in airlines:
        if company[3] != "":
            new_airline_liste = new_airline_liste+[company]

    for company in new_airline_liste:
        if company[1] == airline:
            company_name = company[3]

    for airport in airports:
        airo = airo+[[airport[4], airport[6], airport[7]]]
        if airport[4] != '\\N':
            varhi = varhi+[airport[4]]

    for ligne in routescsv:
        if ligne[0] == company_name:
            MEA = MEA+[[ligne[2], ligne[4]]]

    for listes in MEA:
        for airport in airo:

            if listes[0] == airport[0]:
                dep = dep+[float(airport[1]), float(airport[2])]
            if listes[1] == airport[0]:
                ar = ar+[float(airport[1]), float(airport[2])]

    dep = list(list_split(dep, 2))
    ar = list(list_split(ar, 2))

    i = 0
    long_routes = []
    lat_routes = []

    while i < len(dep):
        lat_routes = lat_routes+[[dep[i][0], ar[i][0]]]
        long_routes = long_routes+[[dep[i][1], ar[i][1]]]
        i = i+1

    fig = plt.figure(facecolor='#b0b29e')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    #fig.set_size_inches(10, 5.5)

    j = 0

    while j < len(long_routes):
        plt.plot(long_routes[j], lat_routes[j],
                 color=(204/255.0, 0, 153/255.0, 0.7), linewidth=0.25, marker='.', transform=ccrs.Geodetic()

                 )
        j = j+1

    plt.savefig('static/resultats/routeairline.png', dpi=1200)

    fig = plt.figure(facecolor='#b0b29e')
    ax = plt.axes(projection=ccrs.Robinson())
    ax.stock_img()
    #fig.set_size_inches(10, 5.5)

    j = 0

    while j < len(long_routes):
        plt.plot(long_routes[j], lat_routes[j],
                 color=(204/255.0, 0, 153/255.0, 0.7), linewidth=0.25, marker='.', transform=ccrs.Geodetic()

                 )
        j = j+1

    plt.savefig('static/resultats/routeZoomairline.png', dpi=1200)
