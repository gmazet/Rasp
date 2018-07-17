########################
THISEVID=633319 #Delaware
THISEVID=634492 # S. Calif
THISEVID=635957 # Java
THISEVID=630830 # St Jean de Maurienne
THISEVID=639347 # Montenegro
THISEVID=639343 # San Francisco
########################

french_stations=['R9F1B','RE4C0','R57C7','R0B29','RA3B7']
french_stations=['R9F1B','RE4C0']


from optparse import OptionParser
from sys import exit

#----------------------------------------------------------------
def MyOptParser(parser):
    parser.add_option("--evid", action="store", dest="evid", help="Earthquake EVID")
    parser.add_option("--otime", action="store", dest="otime", help="Earthquake origin time '%Y-%m-%dT%H%M%S.%fZ' (if no EVID)")
    parser.add_option("--lat", action="store", dest="latitude", help="Epicenter latitude (if no EVID)")
    parser.add_option("--lon", action="store", dest="longitude", help="Epicenter longitude (if no EVID)")
    parser.add_option("--depth", action="store", dest="depth", default=2.0, help="Hypocenter depth (if no evid, optional)")
    parser.add_option("--mag", action="store", dest="magnitude", default=2.0, help="Earthquake magnitude (if no evid, optional)")
    parser.add_option("--length", action="store", dest="sig_length", default=120, help="Length of the signal to show in seconds [120]")
    parser.add_option("--fmin", action="store", dest="freqmin", help="Bandpass low frequency [0.5]")
    parser.add_option("--fmax", action="store", dest="freqmax", help="Bandpass high frequency [10.0]")
    parser.add_option("--ampmax", action="store", dest="ampmax", default=5000.0, help="Maximum amplitude to show [5000]")
    parser.add_option("--nbsta", action="store", dest="maxnbsta", default=3, help="Nb stations to show [3]")
    parser.add_option("--force", action="store_true", dest="force", default=False, help="Force requesting data instead of reading local miniseed file")
    parser.add_option("--fr", action="store_true", dest="fr", default=False, help="Take French stations only")
    parser.add_option("--section", action="store_true", dest="section", default=False, help="Show section")
    parser.add_option("--counts", action="store", dest="counts", default=0, help="peak-peak amplitude read as counts on the Rasp sensor")
    (options, args) = parser.parse_args(args=None, values=None)
    return options

parser=OptionParser()
try:
    MyOptions=MyOptParser(parser)
except:
    print("Incorrect options. Try --help")
    exit(1)

import os

from numpy import log,pi
Rearth=6371.0
km2deg=180/(pi*Rearth)
deg2km=(pi*Rearth)/180

import datetime,time


DATADIR="%s" % os.environ['TMP']
STATIONBOOK="%s/stations_raspberryshake.txt" % DATADIR
print "STATIONBOOK=%s" % STATIONBOOK
NODATA=('R6E96','R052F','R1FBA','RB511','SFD6B', 'RA70D','R9DA3','R5661','RE7F9','R9DBD','R3F3B')

HOST1="gmazet.freeboxos.fr"
HOST2="manchot.emsc-csem.org"
HOST3="rtserve.iris.washington.edu"
HOST4="rtserver.ipgp.fr"
HOST5="geofon.gfz-potsdam.de"
HOST6="geosrt1.ipgp.fr"
HOST6="eida.ipgp.fr"
HOST7="caps.raspberryshakedata.com"
HOST8="fdsnws.raspberryshakedata.com"

RASP_NET="AM"
RASP_LOCALCODE="00"
RASP_CHANNEL="SHZ"
RASP_HOST=HOST8

BEFORE=10

def TIMESTAMP_TO_DATETIME(timestamp):
    return datetime.datetime.fromtimestamp(timestamp+0.05)

#----------------------------------------------------------------------
class myEVID():
    def __init__(self):
        self.yr=None
        self.month=None
        self.day=None
        self.hour=None
        self.min=None
        self.sec=None
        self.lat=None
        self.lon=None
        self.depth=None
        self.fdep=None
        self.mag=None
        self.magtype=None
        self.automan=None
        self.net2=None
        self.net34=None
        self.fid=None
        self.group=None
        self.oritime=None
        self.updtime=None
        self.evid=None
        self.region=None
        self.tzid=None
        self.color=None


# Raspberryshake station coordinates,infos
class req_station:
    def __init__(self):
        self.network=""
        self.name=""
        self.location=""
        self.channel=""
        self.lat=0
        self.lon=0
        self.elevation=0
        self.slserver="rtserve.iris.washington.edu"

    def get_epidist(self,epilat,epilon):
        D,Az1,Az2=gps2dist_azimuth(self.lat, self.lon, epilat, epilon)
        self.epidist_km=D/1000
        self.epidist_deg=self.epidist_km*km2deg
        self.azimuth=Az2 # Epi -> Sta azimuth
        return self.epidist_km,self.epidist_deg, self.azimuth
        


#-------------------------------------------------------------------
def counts2Amp (c):
    # Av=470*c ( 470 counts ~ 1 um/s according to sensor sensentivity provided by sensor technical specs)
    # log(Ad/T)=Av/2pi and  T ~ 1.0 sec
    Av=float(c)/470
    Ad=pow(10,(Av/(2*pi)))
    return Av,Ad
    
#-------------------------------------------------------------------
def get_data(listofstations,ev):
    allsta={}
    arrtimes={}
    alltraces=Stream()

    global my_phase_list
    global freqmin
    global freqmax

    ista=0
    for mysta in listofstations:
        ista+=1
        sta=req_station()
        sta.network, sta.name, sta.location, sta.channel, sta.lat, sta.lon, sta.elevation, sta.slserver, sta.epidist_deg, sta.azimuth = mysta
        if (sta.name in NODATA):
            ista-=1
            continue

        sta.epidist_km=deg2km*sta.epidist_deg
        #########################
        #if (sta.epidist_deg<2.1):
            #ista-=1
            #continue
        #########################

        outseedfile="%s/%d.%s_%s_%s_%s.%s.mseed" % (DATADIR,evid,sta.network,sta.name,sta.location,sta.channel,ev.oritime)

        print ""
        print "== Request %s (%d/%d)" % (sta.name,ista,maxnbsta)
        allsta[sta.name]=sta

        if (ev.depth<0):
            ev.depth=1
        if (sta.epidist_deg<10):
            my_phase_list=["P","Pg", "Pn", "Sn", "Sg", "PP"]
        else:
            if (ev.depth>50):
                my_phase_list=["Pn", "Pg", "pP", "P", "PP", "PKP"]
            else:
                my_phase_list=["Pn", "Pg", "Sn", "Sg", "P", "PP", "PKP"]

        print "	Distance: %.1f km ; %.2f  degrees" % (sta.epidist_km,sta.epidist_deg)

        #print "Arrival times:"
        arrivals = model.get_travel_times(source_depth_in_km=ev.depth, distance_in_degree=sta.epidist_deg, phase_list=my_phase_list)
        arrtimes[sta]=[]
        for arr in arrivals:
            #print arr
            arrtimes[sta].append((arr.name,TIMESTAMP_TO_DATETIME(ev.oritime+arr.time)))

        # Just a test
        # -----------------------
        if ((counts>0) & (sta.name=="R9F1B")):
            Av,Ad=counts2Amp(counts)
            print "	Max velocity amplitude: 		%.1f um/s" % Av
            print "	Max displacement amplitude: 	%d um" % int(Ad)

        # If if local mseed file does not exist or option "force" is True, query data and save miniseed data file
        # -----------------------

        if ((not os.path.exists(outseedfile)) | (MyOptions.force)):
            #SLC = FDSN("http://fdsnws.raspberryshakedata.com", timeout=5)
            try:
                if (sta.slserver=="fdsnws.raspberryshakedata.com"):
                    SLC = FDSN("http://fdsnws.raspberryshakedata.com", timeout=5)
                
                else:
                    print ">> SL %s" % sta.slserver
                    SLC=SL(server=sta.slserver,timeout=5)
            except:
                ista-=1
                continue

            # First arrival time
            if (MyOptions.section): # IF section, all traces must start at the same time
                AT=OT
            else:
                if (len(arrivals)>0):
                    AT=OT+arrivals[0].time
                else:
                    print ("No theoretical wave found")
                    AT=OT

            t1 = AT - BEFORE
            t2 = t1 + sig_length

            try:
                st = SLC.get_waveforms(sta.network, sta.name, sta.location, sta.channel, t1, t2)
            except:
                print "Can't find data for station %s" % sta.name
                ista-=1
                continue

            if (len(st)>0):
                try:
                    print "	Save trace in %s" % outseedfile
                    st.write(outseedfile, format='MSEED') 
                except:
                    print "	ERROR: Can't save trace in %s" % outseedfile
                    ista-=1
                    continue
            else:
                print "	!!!!!! NO DATA FOUND !!!!!!!"
                ista-=1
                continue


        # If miniseed data file already exists, read it !
        # -----------------------
        else:
            print "	Mseed file %s already exist. Don't write over" % outseedfile
            # First arrival time
            if (MyOptions.section): # IF section, all traces must start at the same time
                AT=OT
            else:
                if (len(arrivals)>0):
                    AT=OT+arrivals[0].time
                else:
                    print ("No theoretical wave found")
                    AT=OT
    
            t1 = AT - BEFORE
            t2 = t1 + sig_length
            #print "\t",ista,t1,t2
            st=read(outseedfile,starttime=t1,endtime=t2)

        # Merge and try to detrend.
        # -----------------------
        st.merge()
        print st

        try:
            st.detrend()
        except:
            ista-=1
            continue
        
        if ((MyOptions.freqmin == None) & (MyOptions.freqmax == None)):
            freqmin=0.0
            freqmax=25.0
            #st.filter("highpass", freq=freqmin)
        elif ((MyOptions.freqmin != None) & (MyOptions.freqmax == None)):
            freqmin=float(MyOptions.freqmin)
            freqmax=25
            st.filter("highpass", freq=freqmin)
        elif ((MyOptions.freqmax != None) & (MyOptions.freqmin == None)):
            freqmin=0
            freqmax=float(MyOptions.freqmax)
            st.filter("lowpass", freq=float(MyOptions.freqmax))
        elif ((MyOptions.freqmax == None) & (MyOptions.freqmin == None)):
            freqmin=0
            freqmax=25
        else:
            freqmin=float(MyOptions.freqmin)
            freqmax=float(MyOptions.freqmax)
            st.filter("bandpass", freqmin=freqmin, freqmax=freqmax) 

        if (len(st)>0):
            st[0].stats.coordinates = AttribDict()
            st[0].stats.coordinates.latitude = sta.lat
            st[0].stats.coordinates.longitude = sta.lon
            st[0].stats.coordinates.elevation = sta.elevation
            st[0].stats.distance = gps2dist_azimuth(sta.lat, sta.lon, ev.lat, ev.lon)[0]
            st[0].stats.distance_deg = st[0].stats.distance/1000+km2deg

            if (ista==1):
                alltraces = st
            else:
                alltraces += st
        else:
            ista-=1
            continue
    
        if (ista>=maxnbsta):
            return allsta, arrtimes, alltraces

#-------------------------------------------------------------------
    
if (MyOptions.evid != None):
    ev=myEVID()
    try:
        evid=int(MyOptions.evid)
    except:
        print "evid must be an integer"
        raise

    from helper_tools import *
    url="http://www.seismicportal.eu/fdsnws/event/1/query?source_id=%d&format=json" % evid
    print "get url %s ..." % url
    res = geturl(url)
    print "done"
    ev1=parsejson(res['content'])
    ev2=ev1['features'][0]['properties']
    ev.lat,ev.lon,ev.mag, ev.depth, ev.OT =  ev2['lat'], ev2['lon'], ev2['mag'],ev2['depth'],ev2['time']
    ev.region = ev2['flynn_region']

    oritime=datetime.datetime.strptime(ev.OT, '%Y-%m-%dT%H:%M:%S.%fZ')
    unixtime=time.mktime(oritime.timetuple())
    ev.oritime=unixtime

    (ev.lat,ev.lon,ev.oritime,ev.mag,ev.depth)=(float(ev.lat),float(ev.lon), float(ev.oritime), float(ev.mag), int(ev.depth))

else:
    if ((MyOptions.latitude != None) & (MyOptions.longitude != None)):
        latitude=MyOptions.latitude
        longitude=MyOptions.longitude
    else:
        print "No epicenter provided. Exit"
        exit()

    if (MyOptions.otime != None):
        otime=MyOptions.otime
    else:
        print "No origin time provided. Exit"
        exit()

    if (MyOptions.depth != None):
        depth=MyOptions.depth
    else:
        depth=0

    if (MyOptions.magnitude != None):
        magnitude=MyOptions.magnitude
    else:
        magnitude=0

    try:
        evid=0
        ev=myEVID()
        ev.lat,ev.lon,ev.depth, ev.mag, ev.OT =  latitude, longitude, depth, magnitude, otime
        ev.region = ""
        oritime=datetime.datetime.strptime(ev.OT, '%Y-%m-%dT%H:%M:%S.%fZ')
        unixtime=time.mktime(oritime.timetuple())
        ev.oritime=unixtime

        (ev.lat,ev.lon,ev.oritime,ev.mag,ev.depth)=(float(ev.lat),float(ev.lon), float(ev.oritime), float(ev.mag), int(ev.depth))

    except:
        print "Missing parameters. Try --help"
        exit()

if (MyOptions.counts != None):
    counts=MyOptions.counts
else:
    counts=0

if (MyOptions.sig_length != None):
    sig_length=float(MyOptions.sig_length)
else:
    sig_length=120.0

#if (MyOptions.freqmin != None):
    #freqmin=float(MyOptions.freqmin)
#else:
    #freqmin=0.5

#if (MyOptions.freqmax != None):
    #freqmax=float(MyOptions.freqmax)
#else:
    #freqmax=10.0

if (MyOptions.ampmax != None):
    ampmax=float(MyOptions.ampmax)
else:
    ampmax=5000.0

if (MyOptions.maxnbsta != None):
    maxnbsta=int(MyOptions.maxnbsta)
else:
    maxnbsta=3


print "import obspy ..."
from numpy import array,argsort
from obspy.core.utcdatetime import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
print "done"

# Build station book
# ------------------
fsta=open(STATIONBOOK,'r')
ALLSTATIONS=[]
i=0
for line in fsta:
    sta=req_station()
    sta.lat,sta.lon,sta.elevation=float(line.split()[1]),float(line.split()[2]),float(line.split()[3])
    sta.get_epidist(ev.lat,ev.lon)
    sta.network, sta.name,  sta.location, sta.channel, sta.slserver = RASP_NET,line.split()[0], RASP_LOCALCODE,RASP_CHANNEL,RASP_HOST
    #if (sta.name != 'R9F1B'):
    #    continue
    i+=1
    ALLSTATIONS.append((sta.network, sta.name, sta.location, sta.channel, sta.lat, sta.lon, sta.elevation, sta.slserver,sta.epidist_deg,sta.azimuth))
    
    #if (i>10):
    #    break
fsta.close()

"""
# Add CLF
lat,lon,elev=48.025790,2.26,100
D,Az1,Az2=gps2dist_azimuth(lat,lon, ev.lat, ev.lon)
ALLSTATIONS.append(("G","CLF","00", "BHZ", lat,lon,elev,"SDS",D/1000*km2deg, Az2))
# Add LOR
lat,lon,elev=47.2683,3.8589, 300.0
D,Az1,Az2=gps2dist_azimuth(lat,lon, ev.lat, ev.lon)
ALLSTATIONS.append(("RD","LOR","", "BHZ", lat,lon,elev,"SDS",D/1000*km2deg, Az2))
"""

ALLSTATIONS=array(ALLSTATIONS,dtype=([('network', 'S5'),('name', 'S10'), ('location', 'S5'), ('channel', 'S10'), ('lat', 'f4'), ('lon', 'f4'), ('elevation,', 'f4'), ('slserver', 'S30'), ('epidist_deg', 'f4'), ('azimuth', 'f4') ]))
isortbydist=argsort(ALLSTATIONS['epidist_deg'])
STATION= ALLSTATIONS[isortbydist]

FR_STATIONS=[]
for sta in ALLSTATIONS:
    if (sta['name'] in french_stations):
        FR_STATIONS.append(sta)


MAXDIST=STATION['epidist_deg'][maxnbsta-1]
print "MAXDIST=%.1f degrees" % MAXDIST

#-------------------------------------
print "===================================="
print "Magnitude: %3.1f" % ev.mag
print "Region: %s" % ev.region
print "Origin time: %s" % TIMESTAMP_TO_DATETIME(ev.oritime)
print "Coord: %.2f %.2f" % (ev.lat,ev.lon)
print "Depth: %d km" % ev.depth
print "===================================="

# Origin time
# --------------------------
OT=UTCDateTime(ev.OT)

print "OT:",OT

from obspy import Stream
from obspy.core import read, AttribDict
from obspy.taup import TauPyModel
model=TauPyModel(model='ak135')

from obspy.clients.fdsn import Client as FDSN
if (MyOptions.fr):
    allsta,arrtimes,alltraces=get_data(FR_STATIONS,ev)
else:
    allsta,arrtimes,alltraces=get_data(STATION,ev)

# ------------------------------------------------------
print "Plot ..."
print "import matplotlib"
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates

mytitle="EVID %d - M%3.1f %s on %s (Lat: %.2f; Lon: %.2f; Z: %dkm)" % (evid,ev.mag,ev.region,str(OT)[0:21],ev.lat,ev.lon,ev.depth)

fig2=plt.figure(figsize=[12,8])

if (MyOptions.section):
    from matplotlib.transforms import blended_transform_factory
    alltraces.plot(fig=fig2, type='section', ev_coord=(ev.lat,ev.lon), dist_degree=True,  time_down=True, linewidth=.35, grid_linewidth=.25, offsetmin=0, recordlength=sig_length)
    axes=fig2.get_axes()
    ax=axes[0]
    
    #YLIM=axes[0].get_ylim()
    #DY=(YLIM[1]-YLIM[0])*1440
    #NEWYLIM=((YLIM[0], YLIM[0] + sig_length/60/1440))
    #ax.set_ylim(0, sig_length)

    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in alltraces:
        #print tr.stats.distance/1e6, tr.stats.station
        ax.text(tr.stats.distance/1e6, 1.0, tr.stats.station, rotation=305, va="bottom", ha="center", transform=transform)
    fig2.suptitle(mytitle,fontsize=11)
    
    #plt.tight_layout(pad=0.5,rect=(0.05,0.05,0.95,0.95))
    fig2.tight_layout(rect=(0,0,1,0.90))

    png="%s/%d.section.png" % (DATADIR,evid)
    plt.savefig(png)
    plt.show()
    exit()

else:
    alltraces.plot(fig=fig2,automerge=False)
    axes=fig2.get_axes()
    XLIM=axes[0].get_xlim()
    DX=(XLIM[1]-XLIM[0])*1440
    NEWXLIM=((XLIM[0], XLIM[0] + sig_length/60/1440))


print
i=0
for trace in alltraces:
    print "== Process trace %s" % trace
    tr=trace.stats
    Xtime_min = trace.stats.starttime
    Xtime_max = trace.stats.endtime
    sampling_rate=trace.stats.sampling_rate
    timestep=(1/sampling_rate)*1000
    npts=int(trace.stats['npts'])
    print ("	Nb points= %d ; sampling rate= %d Hz ; time step= %f ms" % (npts,sampling_rate,timestep))

    trace.trim(Xtime_min, Xtime_min + datetime.timedelta(milliseconds=sig_length*1000))
 
    axes[i].set_xlim(NEWXLIM[0], NEWXLIM[1])
    axes[i].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    axes[i].xaxis.set_major_locator(dates.SecondLocator(interval=30))
    axes[i].xaxis.set_minor_locator(dates.SecondLocator(interval=10))
    labelx = axes[i].get_xticklabels()
    plt.setp(labelx, rotation=30, fontsize=9)

    Ymax=min (max(trace.data)*1.01,ampmax) # 5000 max
    ystep=max ((int(Ymax/1000)+1)*1000/2,100)
    Ymaxunit=len(str(int(Ymax)))
    Ymaxsize=str(int(Ymax))[0]
    ystep=int(Ymaxsize)*pow(10,(int(Ymaxunit)-1)) 
    axes[i].set_ylim(-Ymax,Ymax)
    loc = ticker.MultipleLocator(base=ystep) # this locator puts ticks at regular intervals
    axes[i].yaxis.set_major_locator(loc)
    labely = axes[i].get_yticklabels()
    plt.setp(labely, fontsize=9)

    arrivals = model.get_travel_times(source_depth_in_km=ev.depth, distance_in_degree=allsta[tr.station].epidist_deg, phase_list=my_phase_list)
    arrtimes[tr.station]=[]
    for arr in arrivals:
        arrtimes[tr.station].append((arr.name,TIMESTAMP_TO_DATETIME(ev.oritime+arr.time)))

    phases_done=[]
    for pick in arrtimes[tr.station]:
        if (pick[0] in phases_done):
            continue
        phases_done.append(pick[0])
        phase_pick=pick[1]
        x=[phase_pick,phase_pick]
        y=[-Ymax,0]
        axes[i].plot(x,y)
        offsetx= datetime.timedelta(milliseconds = (NEWXLIM[1]-NEWXLIM[0])*1000*0.01)
        axes[i].text(phase_pick+offsetx, -Ymax, pick[0], style='normal', bbox={'facecolor':'lightblue', 'alpha':0.8, 'pad':4}, fontsize=8)

    streamcode="%s_%s:%s:%s" % (tr.network, tr.station, tr.location, tr.channel)
    axes[i].set_title("%s (%.1f degrees ; %d km ; Azim %d) - Filter [%.1f-%.1f Hz]" % (streamcode,allsta[tr.station].epidist_deg, int(allsta[tr.station].epidist_km), allsta[tr.station].azimuth, freqmin,freqmax), fontsize=9)

    i=i+1

fig2.suptitle(mytitle,fontsize=11)
fig2.tight_layout(pad=0.5,rect=(0,0,1,0.95))

png="%s/%d.png" % (DATADIR,evid)
plt.savefig(png)

plt.show()
exit() ################
