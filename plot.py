import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates
from datetime import datetime, timedelta
import os

def matplotlib_plot(ev, phaseslist, allsta, arrtimes, alltraces, model, options, DATADIR):
    print "plot with matplotlib ..."
    mytitle="EVID %d - M%3.1f %s on %s (Lat: %.2f; Lon: %.2f; Z: %dkm)" % (ev.evid,ev.mag,ev.region,str(ev.OTutc)[0:21],ev.lat,ev.lon,ev.depth)

    fig2=plt.figure(figsize=[12,8])

    if (options.section):
        from matplotlib.transforms import blended_transform_factory
        alltraces.plot(fig=fig2, type='section', ev_coord=(ev.lat,ev.lon), dist_degree=True,  time_down=True, linewidth=.35, grid_linewidth=.25, offsetmin=0, recordlength=float(options.sig_length))
        axes=fig2.get_axes()
        ax=axes[0]
        
        transform = blended_transform_factory(ax.transData, ax.transAxes)
        for tr in alltraces:
            #print tr.stats.distance/1e6, tr.stats.station
            ax.text(tr.stats.distance/1e6, 1.0, tr.stats.station, rotation=305, va="bottom", ha="center", transform=transform)
        fig2.suptitle(mytitle,fontsize=11)
    
        fig2.tight_layout(pad=0.5,rect=(0,0,1,0.90))
    
        png="%s/%d.section.png" % (DATADIR,ev.evid)
        plt.savefig(png)
        plt.show()
        exit()

    else:
        alltraces.plot(fig=fig2,automerge=False)
        axes=fig2.get_axes()
        XLIM=axes[0].get_xlim()
        DX=(XLIM[1]-XLIM[0])*1440
        NEWXLIM=((XLIM[0], XLIM[0] + float(options.sig_length)/60/1440))

    
    print
    i=0
    from obspy import read_inventory, Stream
    for trace in alltraces:
        print "== Process trace %s" % trace
        tr=trace.stats

        try:
            myinv = read_inventory(os.path.expanduser('./responses/%s.xml' % tr.station))

            tr_copy=trace.copy()
            mystream = Stream(traces=[tr_copy])
            mystream1=mystream.copy()
            mystream2=mystream.copy()
            mystream3=mystream.copy()
    
            mystream1.attach_response(myinv)
            mystream1.remove_response(output='ACC')
            mystream2.attach_response(myinv)
            mystream2.remove_response(output='VEL')
            mystream3.attach_response(myinv)
            mystream3.remove_response(output='DISP')

            print "<<<<<<<"
            for mytrace in mystream1:
                print('Station %s PGA: %.4f m/s/s (%.3f mg)' % (mytrace.stats.station, max(abs(mytrace.data)), max(abs(mytrace.data))/9.81*1000))
            for mytrace in mystream2:
                print('Station %s PGV: %.3f micrometers/s' % (mytrace.stats.station, max(abs(mytrace.data))*1000000))
            for mytrace in mystream3:
                print('Station %s PGD: %.3f micrometers' % (mytrace.stats.station, max(abs(mytrace.data))*1000000))
            print "<<<<<<<"

        except:
            print ("WARNING: no response file found for station %s" % tr.station)
            print ("\tDownload response file via Rasp webservice")
            print ("\tExample:  curl -k 'https://fdsnws.raspberryshakedata.com/fdsnws/station/1/query?network=AM&station=RDF31&level=resp&format=sc3ml' -o RDF31.xml")
            pass
           

        print
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
        ystep=max ((int(Ymax/1000)+1)*1000/2,100)
        Ymaxunit=len(str(int(Ymax)))
        Ymaxsize=str(int(Ymax))[0]
        ystep=int(Ymaxsize)*pow(10,(int(Ymaxunit)-1)) 
        axes[i].set_ylim(-Ymax,Ymax)
        loc = ticker.MultipleLocator(base=ystep) # this locator puts ticks at regular intervals
        axes[i].yaxis.set_major_locator(loc)
        labely = axes[i].get_yticklabels()
        plt.setp(labely, fontsize=9)

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

