# Correlation block object:
from obspy import Stream, Trace, read_inventory, UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth

import numpy as np
import os
import re
from glob import glob




class CorrTrace(object):

	"""
	Object holds correlation data along with metainformation (station id, geographic location) in preparation for correlation.
	Context manager: Upon exiting, the binary file holding intermediate correlation traces will be closed properly and the stack will be written.
	"""

	def __init__(self,cha1,cha2,corr_type,nlag,t0=None,t1=None):

		self.stack = np.zeros(nlag)
		self.pstak = np.zeros(nlag)

		self.begin = t0
		self.end   = t1
		self.id1   = cha1
		self.id2   = cha2
		self.id    = cha1 + '--' + cha2
		self.ctype = corr_type


		self.cnt_tot = 0
		self.cnt_int = 0

		geo_inf = self.get_geoinf()
		
			
		

		self.lat1 = geo_inf[0]
		self.lat2 = geo_inf[2]
		self.lon1 = geo_inf[1]
		self.lon2 = geo_inf[3]
		self.az   = geo_inf[5]
		self.baz  = geo_inf[6]
		self.dist = geo_inf[4]


		# open the file to dump intermediate stack results
		int_file = '{}.{}.windows.bin'.format(self.id,self.ctype)
		int_file = os.path.join('data','correlations',int_file)
		self.int_file = open(int_file,'wb')

	def get_geoinf(self):


		inv1 = '{}.{}.xml'.format(*self.id1.split('.')[0:2])
		inv2 = '{}.{}.xml'.format(*self.id2.split('.')[0:2])

		inv1 = read_inventory(os.path.join('meta','stationxml',inv1))
		inv2 = read_inventory(os.path.join('meta','stationxml',inv2))

		# Replace 'radial' and 'transverse' by 'N' and 'E'
		id1 = re.sub('\.??R$','N',self.id1)
		id2 = re.sub('\.??R$','N',self.id2)
		id1 = re.sub('\.??T$','E',id1)
		id2 = re.sub('\.??T$','E',id2)
		

		c1 = inv1.get_coordinates(id1)
		c2 = inv2.get_coordinates(id2)

		lat1, lon1, lat2, lon2 = (
			c1['latitude'],
			c1['longitude'],
			c2['latitude'],
			c2['longitude'])

		dist, az, baz = gps2dist_azimuth(lat1,lon1,lat2,lon2)

		return lat1, lon1, lat2, lon2, dist, az, baz


	def _add_corr(self):
		pass

	def write_tot(self):
		pass

	def write_int(self):
		pass

	def __enter__(self):
		return self

	def __exit__(self,type,value,traceback):

		self.write_tot()
		self.int_file.close()

class CorrBlock(object):


# - initialize with station pairs
	def __init__(self,block,inv,cfg):

		self.inv = inv
		self.correlations = []
		self.channels = []

		n_lag = 13

		for pair in block:
			try:
				self.correlations.append(CorrTrace(pair[0],pair[1],
				cfg.corr_type,n_lag))
			except:
				print('** Could not initialize correlation for %s,%s: check metadata'
					%(pair[0],pair[1]))
			
			self.channels.append(pair[0])
			self.channels.append(pair[1])


		self.channels = list(set(self.channels))
		print(self.channels)
		self.data = Stream()
		self.initialize_data()

		self.correlate(cfg)

	def correlate(self,cfg):

		t_0 = UTCDateTime(cfg.time_begin)
		t_end = UTCDateTime(cfg.time_end)

		

		self.initialize_data()

		t = min(t0,self.data)

		while t <= t_end - cfg.time_window_length:

			#self.update_data(t, cfg.win_len_sec)

			windows = self.data.slice(t, t + cfg.time_window_length)
			print(windows)

			t += cfg.time_window_length - cfg.time_overlap



	def initialize_data(self):
		
		for channel in self.channels:
			
			print(channel)
			f = self.inv.data[channel].pop(0)
			try:
				self.data += read(f)
			except IOError:
				print('** problems reading file %s' 
				%self.inv.data[channel])





			#if 
# - loop over the station pairs with n 'datatraces' (one per ID):

	# - 
	# - check starttime, if necessary, add data from 'later' file
	# - slice the traces
	# - if horizontal components are involved, rotate them
	# - check minimum length requirement
	# - transform everything to FD?
	# - correlate each relevant pair
	# - add to stack
	# - if window counter reaches n_intermediate_stack: save intermediate