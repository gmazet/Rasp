from rasp_utils import *
from myutils import *
from datetime import datetime

def get_data(listofstations, event, model, options, maxnbsta=3, provider="rasp"):
    allsta={}
    arrtimes={}
    alltraces=Stream()

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
        if ((options.provider=='resif') & (sta.network not in ['FR'])):
            ista-=1
            continue
        if ((options.provider=='rasp') & (sta.network not in ['AM'])):
            ista-=1
            continue

        sta.epidist_km=deg2km*sta.epidist_deg

        outseedfile="%s/%d.%s_%s_%s_%s.%s.mseed" % (DATADIR,event.evid,sta.network,sta.name,sta.location,sta.channel,event.oritime)

        print ""
        print "== Request %s %s %s %s (%d/%d)" % (sta.network,sta.name,sta.location,sta.channel,ista,maxnbsta)
        allsta[sta.name]=sta

        if (event.depth<0):
            event.depth=1
        if (sta.epidist_deg<10):
            my_phase_list=["P","Pg", "Pn", "Sn", "Sg", "PP"]
        else:
            if (event.depth>50):
                my_phase_list=["Pn", "Pg", "pP", "P", "PP", "PKP"]
            else:
                my_phase_list=["Pn", "Pg", "Sn", "Sg", "P", "PP", "PKP"]

        print "	Distance: %.1f km ; %.2f  degrees" % (sta.epidist_km,sta.epidist_deg)

        arrivals = model.get_travel_times(source_depth_in_km=event.depth, distance_in_degree=sta.epidist_deg, phase_list=my_phase_list)
        arrtimes[sta]=[]
        for arr in arrivals:
            arrtimes[sta].append((arr.name,datetime.fromtimestamp(event.oritime+arr.time+0.05)))

        # If if local mseed file does not exist or option "force" is True, query data and save miniseed data file
        # -----------------------

        if ((not os.path.exists(outseedfile)) | (options.force)):
            try:
                if (options.provider=="rasp"):
                    FDSN_CLIENT = FDSN("http://%s" % (sta.slserver), timeout=5)
                else:
                    FDSN_CLIENT = FDSN(provider.upper(), timeout=5)
            except:
                ista-=1
                continue

            # First arrival time
            if (options.section): # IF section, all traces must start at the same time
                AT=event.OTutc
            else:
                if (len(arrivals)>0):
                    AT=event.OTutc+arrivals[0].time
                else:
                    print ("    No theoretical wave found")
                    AT=event.OTutc

            t1 = AT - BEFORE
            t2 = t1 + float(options.sig_length)

            try:
                st = FDSN_CLIENT.get_waveforms(sta.network, sta.name, sta.location, sta.channel, t1, t2)
            except:
                print (" Can't find data for station %s" % sta.name)
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
            if (options.section): # IF section, all traces must start at the same time
                AT=event.OTutc
            else:
                if (len(arrivals)>0):
                    AT=event.OTutc+arrivals[0].time
                else:
                    print ("No theoretical wave found")
                    AT=event.OTutc
    
            t1 = AT - BEFORE
            t2 = t1 + float(options.sig_length)
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
        
        if ((options.freqmin == None) & (options.freqmax == None)):
            freqmin=0.0
            freqmax=25.0
        elif ((options.freqmin != None) & (options.freqmax == None)):
            freqmin=float(options.freqmin)
            freqmax=25
            st.filter("highpass", freq=freqmin)
        elif ((options.freqmax != None) & (options.freqmin == None)):
            freqmin=0
            freqmax=float(options.freqmax)
            st.filter("lowpass", freq=float(options.freqmax))
        else:
            freqmin=float(options.freqmin)
            freqmax=float(options.freqmax)
            st.filter("bandpass", freqmin=freqmin, freqmax=freqmax) 

        if (len(st)>0):
            st[0].stats.coordinates = AttribDict()
            st[0].stats.coordinates.latitude = sta.lat
            st[0].stats.coordinates.longitude = sta.lon
            st[0].stats.coordinates.elevation = sta.elevation
            st[0].stats.distance = gps2dist_azimuth(sta.lat, sta.lon, event.lat, event.lon)[0]
            st[0].stats.distance_deg = st[0].stats.distance/1000+km2deg

            if (ista==1):
                alltraces = st
            else:
                alltraces += st
        else:
            ista-=1
            continue
    
        if (ista>=maxnbsta):
            return my_phase_list,allsta, arrtimes, alltraces

