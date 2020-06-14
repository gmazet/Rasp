import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates
from datetime import datetime, timedelta
import os

def matplotlib_plot(ev, phaseslist, allsta, arrtimes, alltraces, model, options, DATADIR):
    print "plot with matplotlib ..."
    mytitle="EVID %d - M%3.1f %s on %s (Lat: %.2f; Lon: %.2f; Z: %dkm)" % (ev.evid,ev.mag,ev.region,str(ev.OTutc)[0:21],ev.lat,ev.lon,ev.depth)

    fig2=plt.figure(figsize=[11,5])
    #fig, (ax1,ax2,ax3) = plt.subplots(3)

    axes=fig2.get_axes()

    print
    i=0
    from obspy import read_inventory, Stream
    for trace in alltraces:
        print "== Process trace %s" % trace
        tr=trace.stats

        try:
            myinv = read_inventory(os.path.expanduser('./responses/%s.xml' % tr.station))

            trace.attach_response(myinv)

            tr_copy=trace.copy()
            mystream = Stream(traces=[tr_copy])
            mystream_acc=mystream.copy()
            mystream_vel=mystream.copy()
            mystream_disp=mystream.copy()
    
            #mystream_acc.attach_response(myinv)
            mystream_acc.remove_response(output='ACC')
            #mystream_vel.attach_response(myinv)
            mystream_vel.remove_response(output='VEL')
            #mystream_disp.attach_response(myinv)
            mystream_disp.remove_response(output='DISP')

            print "<<<<<<<"
            for mytrace in mystream_acc:
                print('Station %s PGA: %.4f m/s/s (%.3f mg)' % (mytrace.stats.station, max(abs(mytrace.data)), max(abs(mytrace.data))/9.81*1000))
            for mytrace in mystream_vel:
                print('Station %s PGV: %.3f micrometers/s' % (mytrace.stats.station, max(abs(mytrace.data))*1000000))
            for mytrace in mystream_disp:
                print('Station %s PGD: %.3f micrometers' % (mytrace.stats.station, max(abs(mytrace.data))*1000000))
            print "<<<<<<<"

        except:
            print ("WARNING: no response file found for station %s" % tr.station)
            print ("\tDownload response file via associated webservice")
            print ("\tExample with Raspberryshake webservice: curl -k 'https://fdsnws.raspberryshakedata.com/fdsnws/station/1/query?network=AM&station=RDF31&level=resp&format=sc3ml' -o RDF31.xml")
            pass
           
        if (options.output=='vel'):
            print 'Plot velocity'
            trace=mystream_vel[0]
        elif (options.output=='disp'):
            print 'Plot displacement'
            trace=mystream_disp[0]
        else:
            print 'Plot raw data'

        trace.plot(fig=fig2,automerge=False)

        print 'Plot spectrogram'
        clip=[0.0,0.4]
        ###pre_filt=[0.01,0.05,25,30]
        trace.spectrogram(log=False, mult=8.0, wlen=5.0, per_lap=0.9, clip=clip)

        XLIM=axes[0].get_xlim()
        DX=(XLIM[1]-XLIM[0])*1440
        NEWXLIM=((XLIM[0], XLIM[0] + float(options.sig_length)/60/1440))

        Xtime_min = trace.stats.starttime
        Xtime_max = trace.stats.endtime
        sampling_rate=trace.stats.sampling_rate
        timestep=(1/sampling_rate)*1000
        npts=int(trace.stats['npts'])
        print ("	Nb points= %d ; sampling rate= %d Hz ; time step= %f ms" % (npts,sampling_rate,timestep))

        trace.trim(Xtime_min, Xtime_min + timedelta(milliseconds=float(options.sig_length)*1000))
 
        axes[i].set_xlim(NEWXLIM[0], NEWXLIM[1])
        axes[i].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
        axes[i].xaxis.set_major_locator(dates.SecondLocator(interval=30))
        axes[i].xaxis.set_minor_locator(dates.SecondLocator(interval=10))
        labelx = axes[i].get_xticklabels()
        plt.setp(labelx, rotation=30, fontsize=9)

        Ymax=min (max(trace.data)*1.01,int(options.ampmax)) # 5000 max
        print "Ymax=",Ymax
        ##ystep=max ((int(Ymax/1000)+1)*1000/2,100)
        ##print "Ystep=",ystep
        ##Ymaxunit=len(str(int(Ymax)))
        ##Ymaxsize=str(int(Ymax))[0]
        #ystep=int(Ymaxsize)*pow(10,(int(Ymaxunit)-1)) 
        ##print "Ystep=",ystep
        #axes[i].set_ylim(-Ymax,Ymax)
        #loc = ticker.MultipleLocator(base=ystep) # this locator puts ticks at regular intervals
        #axes[i].yaxis.set_major_locator(loc)
        labely = axes[i].get_yticklabels()
        plt.setp(labely, fontsize=9)

        #Ymax1=min (max(alltraces[0].data)*1.01,int(options.ampmax)) # 5000 max
        #Ymax2=max(alltraces[1].data)*1.01
        #Ymax3=min (max(alltraces[2].data)*1.01,int(options.ampmax)) # 5000 max
        #print Ymax1, Ymax2, Ymax3
        #Ymax=Ymax2

        #ystep=max ((int(Ymax/1000)+1)*1000/2,100)
        #Ymaxunit=len(str(int(Ymax)))
        #Ymaxsize=str(int(Ymax))[0]
        #ystep=int(Ymaxsize)*pow(10,(int(Ymaxunit)-1)) 
        #axes[i].set_ylim(-Ymax,Ymax)
        #loc = ticker.MultipleLocator(base=ystep) # this locator puts ticks at regular intervals
        #axes[i].yaxis.set_major_locator(loc)
        #labely = axes[i].get_yticklabels()
        #plt.setp(labely, fontsize=9)

 
        arrivals = model.get_travel_times(source_depth_in_km=ev.depth, distance_in_degree=allsta[tr.station].epidist_deg, phase_list=phaseslist)
        arrtimes[tr.station]=[]
        for arr in arrivals:
            arrtimes[tr.station].append((arr.name,datetime.fromtimestamp(ev.oritime+arr.time+0.05)))

        phases_done=[]
        for pick in arrtimes[tr.station]:
            if (pick[0] in phases_done):
                continue
            phases_done.append(pick[0])
            phase_pick=pick[1]
            x=[phase_pick,phase_pick]
            y=[-Ymax,0]
            axes[i].plot(x,y)
            offsetx= timedelta(milliseconds = (NEWXLIM[1]-NEWXLIM[0])*1000*0.01)
            axes[i].text(phase_pick+offsetx, -Ymax, pick[0], style='normal', bbox={'facecolor':'lightblue', 'alpha':0.8, 'pad':4}, fontsize=8)

        streamcode="%s_%s:%s:%s" % (tr.network, tr.station, tr.location, tr.channel)
        axes[i].set_title("%s (%.1f degrees ; %d km ; Azim %d) - Filter [%.1f-%.1f Hz]" % (streamcode,allsta[tr.station].epidist_deg, int(allsta[tr.station].epidist_km), allsta[tr.station].azimuth, float(options.freqmin),float(options.freqmax)), fontsize=9)
    
        i=i+1

    fig2.suptitle(mytitle,fontsize=11)
    fig2.tight_layout(pad=0.5,rect=(0,0,1,0.95))
    
    png="%s/%d.png" % (DATADIR,ev.evid)
    plt.savefig(png)

    plt.show()

