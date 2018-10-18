from myutils import *

my_stations=['R9F1B','RE4C0','R57C7','R0B29','RA3B7']
my_stations=['RA007']
my_stations=['R573E','RD51C','RA007','RB707']
my_stations=['R9F1B','RDF31','R8F35','RE4C0','R57C7','R0B29','RA3B7']

NODATA=('R6E96','R052F','R1FBA','RB511','SFD6B', 'RA70D','R9DA3','R5661','RE7F9','R9DBD','R3F3B','ILLF')

HOST1="gmazet.freeboxos.fr"
HOST2="rtserve.resif.fr"
HOST3="rtserve.iris.washington.edu"
HOST5="geofon.gfz-potsdam.de"
HOST6="geosrt1.ipgp.fr"
HOST6="rtserver.ipgp.fr"
HOST6="eida.ipgp.fr"
HOST7="caps.raspberryshakedata.com"
HOST8="fdsnws.raspberryshakedata.com"

BEFORE=10
myradius=10.0

#-------------------------------------------------------------------
def counts2Amp (c): #Only for Raspberryshake data (not reliable)
    # Av=470*c ( 470 counts ~ 1 um/s according to sensor sensentivity provided by sensor technical specs)
    # log(Ad/T)=Av/2pi and  T ~ 1.0 sec
    Av=float(c)/470
    Ad=pow(10,(Av/(2*pi)))
    return Av,Ad
    

# Build station book
# ------------------
def build_station_list(ev,provider="resif"):
    i=0
    ALLSTATIONS=[]
    if (provider=="rasp"):
        NET="AM"
        LOCCODE="00"
        CHANNEL="SHZ"
        HOST=HOST8
        RASP_STATIONBOOK="%s/stations_raspberryshake.txt" % DATADIR
        print "RASP_STATIONBOOK=%s" % RASP_STATIONBOOK
        fsta=open(RASP_STATIONBOOK,'r')
        for line in fsta:
            sta=req_station()
            sta.lat,sta.lon,sta.elevation=float(line.split()[1]),float(line.split()[2]),float(line.split()[3])
            sta.get_epidist(ev.lat,ev.lon)
            sta.network, sta.name, sta.location, sta.channel, sta.slserver = NET,line.split()[0], LOCCODE,CHANNEL,HOST
            ALLSTATIONS.append((sta.network, sta.name, sta.location, sta.channel, sta.lat, sta.lon, sta.elevation, sta.slserver,sta.epidist_deg,sta.azimuth))
            i+=1
        fsta.close()

    else:
        from obspy.clients.fdsn import RoutingClient,Client
        #client=RoutingClient("iris-federator")
        #client=RoutingClient("eida-routing")
        client=Client(provider)
        inv = client.get_stations( channel="HHZ", starttime=ev.OTutc, endtime=ev.OTutc, latitude=ev.lat, level="channel", longitude=ev.lon, maxradius=myradius, includerestricted=False)
        #print inv

        ALLCHAN=inv.get_contents()['channels']
        for c in ALLCHAN:
            sta=req_station()
            sta.network,sta.name,sta.location,sta.channel=c.split(".")
            sta.coord= inv.get_coordinates(c,ev.OTutc)
            sta.lat=sta.coord["latitude"]
            sta.lon=sta.coord["longitude"]
            sta.elevation=sta.coord["elevation"]
            sta.get_epidist(ev.lat,ev.lon)
            print sta.network,sta.name,sta.location,sta.channel, sta.coord, sta.epidist_deg
            sta.slserver=""
            ALLSTATIONS.append((sta.network, sta.name, sta.location, sta.channel, sta.lat, sta.lon, sta.elevation, sta.slserver,sta.epidist_deg,sta.azimuth))
            i+=1

    ALLSTATIONS=array(ALLSTATIONS,dtype=([('network', 'S5'),('name', 'S10'), ('location', 'S5'), ('channel', 'S10'), ('lat', 'f4'), ('lon', 'f4'), ('elevation,', 'f4'), ('slserver', 'S30'), ('epidist_deg', 'f4'), ('azimuth', 'f4') ]))
    isortbydist=argsort(ALLSTATIONS['epidist_deg'])
    CLOSE_STATIONS= ALLSTATIONS[isortbydist]

    MY_STATIONS=[]
    for sta in ALLSTATIONS:
        if (sta['name'] in my_stations):
            MY_STATIONS.append(sta)

    return CLOSE_STATIONS,MY_STATIONS
