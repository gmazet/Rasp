import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import datetime
from myutils import DATADIR, TIMESTAMP_TO_DATETIME

print

def matplotlib_plot(ev, phaseslist, allsta, arrtimes, alltraces, model, options):
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
    for trace in alltraces:
        print "== Process trace %s" % trace
        tr=trace.stats
        Xtime_min = trace.stats.starttime
        Xtime_max = trace.stats.endtime
        sampling_rate=trace.stats.sampling_rate
        timestep=(1/sampling_rate)*1000
        npts=int(trace.stats['npts'])
        print ("	Nb points= %d ; sampling rate= %d Hz ; time step= %f ms" % (npts,sampling_rate,timestep))

        trace.trim(Xtime_min, Xtime_min + datetime.timedelta(milliseconds=float(options.sig_length)*1000))
 
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
        axes[i].set_title("%s (%.1f degrees ; %d km ; Azim %d) - Filter [%.1f-%.1f Hz]" % (streamcode,allsta[tr.station].epidist_deg, int(allsta[tr.station].epidist_km), allsta[tr.station].azimuth, float(options.freqmin),float(options.freqmax)), fontsize=9)
    
        i=i+1

    fig2.suptitle(mytitle,fontsize=11)
    fig2.tight_layout(pad=0.5,rect=(0,0,1,0.95))
    
    png="%s/%d.png" % (DATADIR,ev.evid)
    plt.savefig(png)

    plt.show()

def bokeh_plot(ev, phaseslist, allsta, arrtimes, alltraces, model, options):
    print "plot with Bokeh ..."
    from bokeh.plotting import figure, output_file, show, save
    from bokeh.layouts import column
    from bokeh.models import DatetimeTickFormatter
    import pandas as pd
    from numpy import asarray,pi
    import time
    from calendar import timegm

    figs=[]

    TOOLS="hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,"
    TOOLS="pan,wheel_zoom,box_zoom,reset,save,"

    i=1
    for st in alltraces:
        starttime=st.stats.starttime# 2018-02-27T06:29:57.960000Z
        endtime=st.stats.endtime# 2018-02-27T06:31:02.880000Z

        start_utc_time = time.strptime(str(starttime), "%Y-%m-%dT%H:%M:%S.%fZ")
        start_epoch_time = timegm(start_utc_time)

        mydatetime=[]
        mytunix=[]
        for tt in st.times():
            tunix=start_epoch_time + tt
            mytunix.append(tunix)
            mydatetime.append(pd.to_datetime(tunix, unit='s'))

        mydatetime=asarray(mydatetime)
        mytunix=asarray(mytunix)

        #print type(st.times()[0]), type(mydatetime[0])
    
        npts=float(st.stats.npts)
        #print "%d samples" % npts

        df = pd.DataFrame(st.data,mydatetime,columns=['counts'])
        df['tunix']=mytunix

        #title=allsta[st.stats.station]
        title=st.stats.station
        title="EVID %d - M%3.1f %s on %s (Lat: %.2f; Lon: %.2f; Z: %dkm)" % (ev.evid,ev.mag,ev.region,str(ev.OTutc)[0:21],ev.lat,ev.lon,ev.depth)
        streamcode="%s_%s:%s:%s" % (st.stats.network, st.stats.station, st.stats.location, st.stats.channel)
        title="%s (%.1f degrees ; %d km ; Azim %d) - Filter [%.1f-%.1f Hz]" % (streamcode,allsta[st.stats.station].epidist_deg, int(allsta[st.stats.station].epidist_km), allsta[st.stats.station].azimuth, float(options.freqmin),float(options.freqmax))
        print "Title: %s" % title
        if (i>1):
            p = figure(title=title,tools=TOOLS,plot_width=1200, plot_height=300, x_axis_type="datetime",x_range=x_range)
        else:
            p = figure(title=title,tools=TOOLS,plot_width=1200, plot_height=300, x_axis_type="datetime")
        #p.line(df.index, df['counts'], color='navy', alpha=0.5)
        p.line(mydatetime, df['counts'], color='navy', alpha=0.5)
        figs.append(p)

        p.xaxis.formatter=DatetimeTickFormatter(
            milliseconds=["%H:%M:%S.%3N"],
            seconds=["%d/%m/%Y %T"],
            minsec=["%d/%m/%Y %T"],
            minutes=["%d/%m/%Y %T"],
            hourmin=["%d/%m/%Y %T"],
            hours=["%d/%m/%Y %T"],
            days=["%d/%m/%Y"],
            months=["%d/%m/%Y"],
            years=["%d/%m/%Y"],
        )
        p.xaxis.major_label_orientation = pi/4
        if (i==1):
            x_range=p.x_range
            y_range=p.y_range

        i+=1

    output_file("rasp.html")
    
    save (column(figs))

