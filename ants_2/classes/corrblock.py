# Correlation block object:
from __future__ import print_function
from obspy import Stream, UTCDateTime, read
from obspy.signal.invsim import cosine_taper
from scipy.signal import sosfilt, detrend
import numpy as np
import os, time
import re


from ants_2.tools.bookkeep import name_correlation_file
from ants_2.tools.util import get_geoinf
from ants_2.classes.corrtrace import CorrTrace
from ants_2.tools.correlations import cross_covar_fd
from ants_2.tools.treatment import ram_norm, whiten, cap, bandpass
# list of possible channels combinations indicating that the data needs to be rotated.
horizontals = ['RR','RT','TR','TT','TZ','ZT','RZ','ZR']
import matplotlib.pyplot as plt


class CorrBlock(object):


# - initialize with station pairs
    def __init__(self,block,cfg):


        self.cfg = cfg
        self._correlations = {}


        self.inv = block.inventory
        self.channels = block.channels
        self.station_pairs = block.station_pairs
        self.channel_pairs = block.channel_pairs
        self.tstart = UTCDateTime(cfg.time_begin)
        self.tend = UTCDateTime(cfg.time_end)


        if any([i in cfg.corr_tensorcomponents for i in horizontals]):

            self.baz1 = []
            self.baz2 = []

            for pair in self.station_pairs:
                cha1 = [cha for cha in self.channels if pair[0] in cha]
                cha2 = [cha for cha in self.channels if pair[1] in cha]
                geoinf = get_geoinf(cha1[0],cha2[0])
                # - find azimuth, backazimuth
                self.baz2.append(geoinf[5])
                self.baz1.append(geoinf[6])
                    

        self.newdata = Stream()
        self.data = []
        
        for sp_i in range(len(self.station_pairs)):          
            self.data.append(Stream())
            
        self.read_data(self.tend)
        self.rotate_data()
        self.sampling_rate = self.data[0][0].stats.sampling_rate
        self.delta = self.data[0][0].stats.delta
        

        for cp in block.channel_pairs:
            for pair in cp:

                cp_name = '{}--{}'.format(*pair)
                #preprstring = self.get_prepstring()

                
                self._correlations[cp_name] = CorrTrace(pair[0],pair[1],
                self.sampling_rate,stck_int=cfg.interm_stack,
                prepstring=self.get_prepstring(),
                window_length=cfg.time_window_length,
                overlap=cfg.time_overlap,corr_params=None)
                
     
        
        
    def run(self,output_file=None):


        print('Working on station pairs:')
        for sta in self.station_pairs:
            print("{}--{}".format(sta[0],sta[1]))

        t_0 = UTCDateTime(self.cfg.time_begin)
        t_end = UTCDateTime(self.cfg.time_end)
        win_len_seconds = self.cfg.time_window_length
        win_len_samples = int(round(win_len_seconds*self.sampling_rate))
        max_lag_samples = int(round(self.cfg.corr_maxlag * self.sampling_rate))
        
        
        if self.cfg.bandpass is not None:
            fmin = self.cfg.bandpass[0]
            fmax = self.cfg.bandpass[1]
            if fmax <= fmin:
                msg = "Bandpass upper corner frequency must be above lower corner frequency."
                raise ValueError(msg)

            order = self.cfg.bandpass[2]
            sos = bandpass(freqmin=fmin,freqmax=fmax,
                df=self.sampling_rate,corners=order)

                    
        # Time loop
        t = t_0
        while t <= t_end - (win_len_seconds - self.delta):
            
            
            print(t,file=output_file)
            
            # checking and reading new files
            if self.read_data(t): #True if new file was read
                self.rotate_data()
                for i in range(len(self.data)):
                    for tr in self.data[i]:
                        if tr.stats.endtime < t:
                            self.data[i].remove(tr)
                        
            
            
            # - station pair loop
            for sp_i in range(len(self.station_pairs)):

                pair = self.station_pairs[sp_i]
                
                # - select traces
                [net1, sta1] = pair[0].split('.')
                [net2, sta2] = pair[1].split('.')
                
                stream = self.data[sp_i].slice(t,t+win_len_seconds)
                    
                str1 = stream.select(network=net1, station=sta1)
                str2 = stream.select(network=net2, station=sta2)
                
                # - channel loop                
                for cpair in self.channel_pairs[sp_i]:
                
                    cp_name = '{}--{}'.format(*cpair)
                    print(cp_name,file=output_file)
                    

                    loc1, cha1 = cpair[0].split('.')[2:4]
                    loc2, cha2 = cpair[1].split('.')[2:4]

                    try:
                        tr1 = str1.select(location=loc1,channel=cha1)[0]
                        tr2 = str2.select(location=loc2,channel=cha2)[0]
                        
                    except IndexError:
                        print("Channel not found",file=output_file)
                        continue
                    
                    if tr1.stats.npts != tr2.stats.npts:
                        print("Traces do not have equal length\n",file=output_file)
                        continue
                    # - check minimum length requirement
                    # - Quite often not fulfilled due to data gaps
                    if tr1.stats.npts < win_len_samples:
                        print("Trace length < window samples\n",file=output_file)
                        continue
                        
                    if tr2.stats.npts < win_len_samples:
                        print("Trace length < window samples\n",file=output_file)
                        continue

                    if True in np.isnan(tr1.data):
                        print("Trace contains nan\n",file=output_file)
                        continue

                    if True in np.isnan(tr2.data):
                        print("Trace contains nan\n",file=output_file)
                        continue

                    if True in np.isinf(tr1.data):
                        print("Trace contains inf\n",file=output_file)
                        continue

                    if True in np.isinf(tr2.data):
                        print("Trace contains inf\n",file=output_file)
                        continue

                    
                    # - correlate
                    tracedata = [detrend(tr1.data,type='constant'),
                                 detrend(tr2.data,type='constant')]

                    taper = cosine_taper(len(tracedata[0]),p=0.05)

                    for i in range(2):
                        #print("working on iteration",cpair)
                        if self.cfg.bandpass is not None:
                            temp = sosfilt(sos,tracedata[i])
                            tracedata[i] = sosfilt(sos,temp[::-1])[::-1]

                        if self.cfg.cap_glitch:
                            tracedata[i] = cap(tracedata[i],self.cfg.cap_thresh)

                        if self.cfg.onebit:
                            tracedata[i] = np.sign(tracedata[i])

                        if self.cfg.ram_norm:
                            tracedata[i] = ram_norm(tracedata[i],self.sampling_rate,self.cfg.ram_window,self.cfg.ram_prefilt)
                    
                        tracedata[i] *= taper
                        tracedata[i] = np.fft.rfft(tracedata[i])
                        
                        if self.cfg.whiten:               
                            tracedata[i] = whiten(tracedata[i],self.sampling_rate,
                                          self.cfg.white_freqmin,
                                          self.cfg.white_freqmax,
                                          self.cfg.white_taper_samples,
                                          self.cfg.white_waterlevel,
                                          self.cfg.white_taper) 
                    
                    correlation = np.conj(tracedata[0])*tracedata[1]
                    #correlation = cross_covar_fd(tracedata[0],tracedata[1],
                    #                             max_lag_samples)#[0]

                    if self.cfg.corr_normalize:
                        ren1 = np.correlate(tracedata[0],tracedata[0],mode='valid')[0]
                        ren2 = np.correlate(tracedata[1],tracedata[1],mode='valid')[0]
                        correlation /= np.sqrt(ren1 * ren2)
                      
                    # - add to stack
                    if np.isnan(np.sum(correlation)) or np.isinf(np.sum(correlation)):
                        print("invalid data in correlation",file=output_file)
                        continue
                    #if len(correlation) == 2 * max_lag_samples + 1:
                        #self._correlations[cp_name]._add_corr(correlation,t)
                            
                    self._correlations[cp_name]._add_corr_fd(correlation,t,output_format=self.cfg.format_output)
                    #print(self._correlations[cp_name].stack_fd)
                    #else:
                    #    print('Empty window or all values zero in window.',
                    #        file=output_file)

    		# - update time
            t += self.cfg.time_window_length - self.cfg.time_overlap
#            if t.julday > jday:
#                print(t.julday,flush=True)
#                jday=t.julday
            #print("iteration time:",time.time()-t0)
        # - Write results
        for corr in self._correlations.values():
            
            corr._ifft_fd_stack()
            
            corr.write_stack(output_format=self.cfg.format_output)
            
        print(time.time())
        print('Finished a correlation block.')




    def rotate_data(self):            
    
        # rotate data
        for sp_i in range(len(self.station_pairs)):
            
            pair = self.station_pairs[sp_i]
                
            # - select traces
            [net1, sta1] = pair[0].split('.')
            [net2, sta2] = pair[1].split('.')
            
            str1 = self.newdata.select(network=net1, station=sta1)
            str2 = self.newdata.select(network=net2, station=sta2)
            
            self.data[sp_i] += str1.select(component="*Z")
            self.data[sp_i] += str2.select(component="*Z")
            
            
            str1_horizontals = (str1.select(component="*N") +
                                str1.select(component="*E"))
            str2_horizontals = (str2.select(component="*N") +
                                str2.select(component="*E"))
            
            if len(str1_horizontals+str2_horizontals) > 4:
                print("error")
                print(str1_horizontals,str2_horizontals)
            
            start=UTCDateTime(0)
            end=UTCDateTime(9e9)
            for tr in (str1_horizontals+str2_horizontals):
                if tr.stats.starttime>start:
                    start=tr.stats.starttime
                if tr.stats.endtime<end:
                    end=tr.stats.endtime                    
        
            if end > start+tr.stats.sampling_rate:
                try:
                    self.data[sp_i] += str1_horizontals.slice(starttime=start,endtime=end).rotate('NE->RT',self.baz1[sp_i])
                    self.data[sp_i] += str2_horizontals.slice(starttime=start,endtime=end).rotate('NE->RT',(self.baz2[sp_i]+180)%360.)
                except:
                    print("rotation failed")
                    print(str1_horizontals)
                    print(str2_horizontals)
                    print("slicing between:",start,end)
                self.data[sp_i]._cleanup()
            
        return()


    def read_data(self,t):
            
        new_data = False
        for tr in self.newdata:
            if tr.stats.endtime<t+1./tr.stats.sampling_rate:
                self.newdata.remove(tr)
                
        for channel in self.channels:
            
            if len(self.newdata.select(id=channel))==0:
            
                try:
                    f = self.inv[channel].pop(0)
                    self.newdata += read(f)
                    new_data = True
                except:
                    continue
                                            
        return new_data               
                     


    def read_and_rotate_data(self,t,win_len):
        
        # throw away already processed data
        #print(self.raw_data)
        for sp_i in range(len(self.station_pairs)):
            self.data[sp_i].trim(starttime=t)
            for tr in self.data[sp_i]:
                if tr.stats.endtime < t-win_len:
                    self.data[sp_i].remove(tr)
                    
        for tr in self.raw_data:
            if tr.stats.endtime < t-win_len:
                try:
                    self.data.remove(tr)
                except:
                    pass
       
        # read new data
        for channel in self.channels:
                
            st = self.raw_data.select(id=channel)

            if st.sort(['endtime'])[-1].stats.endtime < (t + win_len): 
                
                try:
                    f = self.inv[channel].pop(0)
                    self.raw_data += read(f)
                except IOError:
                    print('** problems reading file %s; or empty list' %f)
                    continue
          
        self.raw_data._cleanup()
        
        # rotate data
        for sp_i in range(len(self.station_pairs)):
            
            pair = self.station_pairs[sp_i]
                
            # - select traces
            [net1, sta1] = pair[0].split('.')
            [net2, sta2] = pair[1].split('.')
            
            str1 = self.raw_data.select(network=net1, station=sta1)
            str2 = self.raw_data.select(network=net2, station=sta2)
            
            self.data[sp_i] += str1.select(component="*Z")
            self.data[sp_i] += str2.select(component="*Z")
            if any([i in self.cfg.corr_tensorcomponents for i in horizontals]):
                self.data[sp_i] += self.rotate_traces(str1,self.baz1[sp_i])
                self.data[sp_i] += self.rotate_traces(str2,(self.baz2[sp_i]+180)%360.)
                    
            self.data[sp_i]._cleanup()
                
        return () 
    
            


    def get_prepstring(self):

        prepstring = ''
        if self.cfg.bandpass: 
            prepstring +='b'
        else:
            prepstring += '-'
        if self.cfg.cap_glitch: 
            prepstring += 'g'
        else:
            prepstring += '-'    
        if self.cfg.whiten: 
            prepstring += 'w'
        else:
            prepstring += '-'
        if self.cfg.onebit: 
            prepstring += 'o'
        else:
            prepstring += '-'
        if self.cfg.ram_norm: 
            prepstring += 'r'
        else:
            prepstring += '-'

        return prepstring
        
        

