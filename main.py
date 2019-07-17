from optparse import OptionParser

#----------------------------------------------------------------
def MyOptParser(parser):
    parser.add_option("--evid", action="store", dest="evid", help="Earthquake EVID")
    parser.add_option("--otime", action="store", dest="otime", help="Earthquake origin time '%Y-%m-%dT%H:%M:%S.%fZ' (if no EVID)")
    parser.add_option("--lat", action="store", dest="latitude", help="Epicenter latitude (if no EVID)")
    parser.add_option("--lon", action="store", dest="longitude", help="Epicenter longitude (if no EVID)")
    parser.add_option("--depth", action="store", dest="depth", default=2.0, help="Hypocenter depth (if no EVID, optional)")
    parser.add_option("--mag", action="store", dest="magnitude", default=2.0, help="Earthquake magnitude (if no EVID, optional)")
    parser.add_option("--length", action="store", dest="sig_length", default=120, help="Length of the signal to show in seconds [120]")
    parser.add_option("--fmin", action="store", dest="freqmin", help="Bandpass low frequency [0.5]")
    parser.add_option("--fmax", action="store", dest="freqmax", help="Bandpass high frequency [10.0]")
    parser.add_option("--ampmax", action="store", dest="ampmax", default=1000000.0, help="Maximum amplitude to show [1e+6]")
    parser.add_option("--nbsta", action="store", dest="maxnbsta", default=3, help="Nb stations to show [3]")
    parser.add_option("--force", action="store_true", dest="force", default=False, help="Force requesting data instead of reading local miniseed file")
    parser.add_option("--section", action="store_true", dest="section", default=False, help="Section plot")
    parser.add_option("--mysta", action="store_true", dest="mysta", default=False, help="Request custom list of stations")
    parser.add_option("--provider", action="store", dest="provider", default="rasp", help="Data provider [rasp]")
    parser.add_option("--acc", action="store", dest="acc", default=0, help="Remove sensor response and convert to acceleration")
    (options, args) = parser.parse_args(args=None, values=None)
    return options

parser=OptionParser()
try:
    MyOptions=MyOptParser(parser)
except:
    print("Incorrect options. Try --help")
    exit(1)

#-------------------------------------------------------------------
from fdsn import *
from plot import *
    
if (MyOptions.evid != None):
    ev=myEVID()
    try:
        evid=int(MyOptions.evid)
        ev.evid=evid
    except:
        print "evid must be an integer"
        raise

    url="http://www.seismicportal.eu/fdsnws/event/1/query?source_id=%d&format=json" % evid
    print "get url %s ..." % url
    res = geturl(url)
    print "done"
    ev1=parsejson(res['content'])
    ev2=ev1['features'][0]['properties']
    ev.lat,ev.lon,ev.mag, ev.depth, ev.OT =  ev2['lat'], ev2['lon'], ev2['mag'],ev2['depth'],ev2['time']
    ev.region = ev2['flynn_region']

    oritime=datetime.strptime(ev.OT, '%Y-%m-%dT%H:%M:%S.%fZ')
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
        ev.evid=evid
        ev.lat,ev.lon,ev.depth, ev.mag, ev.OT =  latitude, longitude, depth, magnitude, otime
        ev.region = ""
        oritime=datetime.strptime(ev.OT, '%Y-%m-%dT%H:%M:%S.%fZ')
        unixtime=time.mktime(oritime.timetuple())
        ev.oritime=unixtime

        (ev.lat,ev.lon,ev.oritime,ev.mag,ev.depth)=(float(ev.lat),float(ev.lon), float(ev.oritime), float(ev.mag), int(ev.depth))

    except:
        print "Missing parameters. Try --help"
        exit()

if (MyOptions.acc != None):
    acc=1
else:
    acc=0

if (MyOptions.sig_length != None):
    sig_length=float(MyOptions.sig_length)
else:
    sig_length=120.0

if (MyOptions.freqmin != None):
    freqmin=float(MyOptions.freqmin)
else:
    #freqmin=0.5
    MyOptions.freqmin=0.5

if (MyOptions.freqmax != None):
    freqmax=float(MyOptions.freqmax)
else:
    #freqmax=10.0
    MyOptions.freqmax=25.0

if (MyOptions.ampmax != None):
    ampmax=float(MyOptions.ampmax)
else:
    ampmax=1000000.0

if (MyOptions.maxnbsta != None):
    maxnbsta=int(MyOptions.maxnbsta)
else:
    maxnbsta=3

if (MyOptions.provider):
    provider=MyOptions.provider
if (provider not in ['rasp','resif']):
    print "Take default provider (Raspberryshake)"
    provider='rasp'

#-------------------------------------
print "===================================="
print "Magnitude: %3.1f" % ev.mag
print "Region: %s" % ev.region
print "Origin time: %s" % datetime.fromtimestamp(ev.oritime)
print "Coord: %.2f %.2f" % (ev.lat,ev.lon)
print "Depth: %d km" % ev.depth
print "===================================="

# Origin time
# --------------------------
ev.OTutc=UTCDateTime(ev.OT)

print "OTutc:",ev.OTutc

from obspy import Stream
from obspy.core import read, AttribDict
from obspy.taup import TauPyModel
model=TauPyModel(model='ak135')

CLOSE_STATIONS,MY_STATIONS=build_station_list(ev,provider)

print MY_STATIONS

if (MyOptions.mysta):
    print "Request only my selection of station..."
    maxnbsta=min(len(MY_STATIONS),maxnbsta)
    phaseslist,allsta,arrtimes,alltraces=get_data(MY_STATIONS,ev,model, MyOptions, maxnbsta, provider)
else:
    phaseslist,allsta,arrtimes,alltraces=get_data(CLOSE_STATIONS,ev,model, MyOptions, maxnbsta, provider)

matplotlib_plot(ev,phaseslist, allsta,arrtimes, alltraces, model, MyOptions, DATADIR)

