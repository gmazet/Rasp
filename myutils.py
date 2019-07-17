from optparse import OptionParser
from sys import exit
import os
import datetime,time

from helper_tools import *

print "import obspy ..."
from numpy import array,argsort
from numpy import log,pi
from obspy.core.utcdatetime import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy import Stream
from obspy.core import read, AttribDict
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client as FDSN
    
Rearth=6371.0
km2deg=180/(pi*Rearth)
deg2km=(pi*Rearth)/180

DATADIR="%s" % os.environ['TMP']

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

# Station coordinates,infos
class req_station:
    def __init__(self):
        self.network=""
        self.name=""
        self.location=""
        self.channel=""
        self.lat=0
        self.lon=0
        self.elevation=0
        self.slserver=""

    def get_epidist(self,epilat,epilon):
        D,Az1,Az2=gps2dist_azimuth(self.lat, self.lon, epilat, epilon)
        self.epidist_km=D/1000
        self.epidist_deg=self.epidist_km*km2deg
        self.azimuth=Az2 # Epi -> Sta azimuth
        return self.epidist_km,self.epidist_deg, self.azimuth

