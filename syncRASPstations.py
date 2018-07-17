import requests,sys,os
import simplejson as json

stationbook="%s/stations_raspberryshake.txt" % os.environ['TMP']
URL="http://raspberryshake.net/query/stations.json"

try:
    r=requests.get(URL, timeout=10)
except Exception as e:
    print(e)
    sys.exit(1)

content=r.text

j=json.loads(content)

net=j['Inventory']['network'][0]
#print net['code']
#print net['restricted']
#print net['shared']
stations=net['station']

#[u'code', u'publicID', u'description', u'affiliation', u'sensorLocation', u'restricted', u'longitude', u'start', u'place', u'shared', u'latitude', u'elevation', u'archive']

#R9F1B  48 41  4.9    2 18 22.9 0.040                                              
#GZT    37 21 19.1   37 33 50.8 1.446                                           DDA 
#WOOW  -25-28-31.8  152  4 47.6 0.128                                               
#GTM    34 17 40.6 -116-21-21.6 0.836                                           PAS  


f=open(stationbook,'wb')
for sta in stations:
    code=sta['code']
    lat=float(sta['latitude'])
    lon=float(sta['longitude'])
    elev=float(sta['elevation'])
    if (elev<0):
        elev=0

    latD=int(lat)
    latM=int((lat-latD)*60)
    latS=(((lat-latD)*60)-latM)*60
    lonD=int(lon)
    lonM=int((lon-lonD)*60)
    lonS=(((lon-lonD)*60)-lonM)*60
    #print sta['code'],sta['latitude'],sta['longitude'],sta['elevation']
    f.write("%s %.4f %.4f %d\n" % (code,lat,lon,elev))
    #print "%-6s%3d%3d%5.1f %4d%3d%5.1f %5.3f                                           RASP" % (code,latD,latM,latS,lonD,lonM,lonS,elev/1000)

f.close()

print len(stations),"Stations"


