import os, re
from glob import glob
from obspy import UTCDateTime
from copy import deepcopy

def find_files(indirs, format):
        
    
    content=list()

    for indir in indirs:

        if "*" in format:
            content.extend(glob(os.path.join(indir,format),recursive=True))
        else:       
            content.extend(glob(os.path.join(indir,'*'+format.lower())))
            content.extend(glob(os.path.join(indir,'*'+format.upper())))
        
    content.sort()
    
    return content

def find_files_irene(indirs,format,masterlist='./AlpArray_MasterStationList.txt'):
    statlist = ['DAVOX']#['A282A','A303A','A101B','A255A','DAVOX','BRMO','BFO']
    masterdict = {}
    with open(masterlist,'r') as f:
        lines = f.readlines()
    for line in lines:
        entry = line.split()
        network = str(entry[1])
        station = str(entry[2])
        channel = str(entry[5].split("_")[0])
        if 'EH' in channel:
            continue
        if 'HN' in channel:
            continue
        try:
            location = str(entry[5].split("_")[1])
        except:
            location = ''
        try:
            masterdict[network]
        except:
            masterdict[network] = {}
        masterdict[network][station] = {}
        masterdict[network][station]['channel'] = channel
        masterdict[network][station]['location'] = location
    content = list()
    print('creating a list of available files...')
    for indir in indirs:
        pathlist_new = [os.path.join(dir_,f) for dir_,_,files in os.walk(indir) if ('HH' in dir_ or 'BH' in dir_) for f in files]
        for i,filepath in enumerate(pathlist_new):
            try:
                network,statname,location,channel,code,year,jday = filepath.split('/')[-1].split('.')
            except:
                print('filename does not follow convention:',filepath)
                continue
            if not(statname in statlist):
                continue
            year = int(year)
            jday = int(jday)
            # check alparray master station list
            try:
                if not (masterdict[network][statname]['channel'] == channel[:2] and masterdict[network][statname]['location'] == location):
                    continue
            except: # if there is no entry in Irene's list, we rather skip this station
                #ignored_statlist.append((network,statname))
                continue
            content.append(filepath)
    content.sort()
    return content

def name_processed_file(stats,startonly=False):
    
    inf = [
        stats.network,
        stats.station,
        stats.location,
        stats.channel
    ]
    
        
    t1=stats.starttime.strftime('%Y.%j.%H.%M.%S')
    t2=stats.endtime.strftime('%Y.%j.%H.%M.%S')
    if startonly: 
        t2 = '*'
    
    inf.append(t1)
    inf.append(t2)

    inf.append(stats._format)
    
    filenew = '{}.{}.{}.{}.{}.{}.{}'.format(*inf)
    
    return filenew
            

def name_correlation_file(sta1,sta2,corr_type,fmt='SAC'):

    #name = '{}--{}.{}.{}'.format(sta1,sta2,corr_type,fmt)
    name = '{}--{}.{}'.format(sta1,sta2,fmt)
   
    return(name)



def file_inventory(cfg):
    
    stations = {}
    data = {}

    # start- and endtime specified in configuration
    t0 = UTCDateTime(cfg.time_begin) 
    t1 = UTCDateTime(cfg.time_end)

    # input directories and format (MSEED, SAC etc)
    indirs = cfg.indirs
    filefmt = cfg.input_format


    # list all files in input directories
    if False:#"*" in filefmt:
        files = find_files_irene(indirs,filefmt)
    else:
        files = find_files(indirs,filefmt)
    
    for f in files:

        # decide whether file fits time range
        fn = os.path.basename(f).split('.')
        st = UTCDateTime('{}-{}T{}:{}:{}'.format(*fn[4:9]))
        et = UTCDateTime('{}-{}T{}:{}:{}'.format(*fn[9:14]))
       
        if st > t1 or et < t0:
            continue
        else:

            station = '{}.{}'.format(*fn[0:2])
            channel = '{}.{}.{}.{}'.format(*fn[0:4])

            # - stations dictionary: What stations exist and what channels do they have
            if station not in stations.keys():
                stations[station] = []

            # - channels dictionary: Inventory of files for each channel
            if channel not in stations[station]:
                stations[station].append(channel)

            if channel not in data.keys():
                data[channel] = []

            data[channel].append(f)
    
    return(stations, data)


def station_pairs(staids,n,autocorr):
   
    #staids = self.stations.keys()
    # sort alphabetically
    staids = list(staids) 
    staids.sort()
    blcks_stations = []
    #blcks_channels = []
    idprs = []

    n_ids = len(staids)
    n_auto = 0 if autocorr else 1
    #n_blk = cfg.n_stationpairs

    for i in range(n_ids):
        for j in range(i+n_auto,n_ids):

            if len(idprs) == n:
                blcks_stations.append(idprs)
                idprs = []

            idprs.append((staids[i],staids[j]))

    if len(idprs) <= n:
        blcks_stations.append(idprs)


    # idprs = []
# 
    # for blck in blcks_stations:
        # idprs = []
        # for pair in blck:
# 
            # idprs.extend(self._channel_pairs(pair[0],pair[1],cfg))  
# 
        # if idprs != []:
            # blcks_channels.append(idprs)

    return blcks_stations


def channel_pairs(channels1,channels2,cfg):


    channels = []
    comp_pairs = []
    cha_pairs = []
    tensor_comp = cfg.corr_tensorcomponents


    for c1 in channels1:
        for c2 in channels2:

            loc1 = c1.split('.')[2]
            loc2 = c2.split('.')[2]

            if loc1 not in cfg.locations:
                continue

            if loc2 not in cfg.locations:
                continue

            if loc1 != loc2 and not cfg.locations_mix:
                continue

            comp_pairs.append([c1[-1],c2[-1]])
            
            cha1 = re.sub('E$','T',c1)
            cha1 = re.sub('N$','R',cha1)
            cha2 = re.sub('E$','T',c2)
            cha2 = re.sub('N$','R',cha2)  
            
            if cfg.update:
                f = name_correlation_file(cha1,cha2,cfg.corr_type)
                f = os.path.join('data','correlations',f)
                if os.path.exists(f):
                    continue
                
            cha_pairs.append((cha1,cha2))
   
     
    for comp in tensor_comp:
        
        if comp == 'RR' or comp == 'TR' or comp == 'RT' or comp == 'TT':
                        
            if (['E','E'] in comp_pairs and ['N','N'] in comp_pairs and
                ['E','N'] in comp_pairs and ['N','E'] in comp_pairs):
                                
                for cha_pair in cha_pairs:
                    
                    if cha_pair[0][-1] == comp[0] and cha_pair[1][-1] == comp[1]:
                        
                        channels.append(cha_pair)
        
                
        elif comp == 'ZR' or comp == 'ZT':
            
            if ['Z','E'] in comp_pairs and ['Z','N'] in comp_pairs:
                
                for cha_pair in cha_pairs:
                    
                    if cha_pair[0][-1] == comp[0] and cha_pair[1][-1] == comp[1]:
                        
                        channels.append(cha_pair)
        
        elif comp == 'RZ' or comp == 'TZ':
            
            if ['E','Z'] in comp_pairs and ['N','Z'] in comp_pairs:
                
                for cha_pair in cha_pairs:
                    
                    if cha_pair[0][-1] == comp[0] and cha_pair[1][-1] == comp[1]:
                        
                        channels.append(cha_pair)
                        
        elif comp == 'ZZ':
            
            if ['Z','Z'] in comp_pairs:
                
                for cha_pair in cha_pairs:
                    
                    if cha_pair[0][-1] == comp[0] and cha_pair[1][-1] == comp[1]:
                        
                        channels.append(cha_pair)

    return(channels)


class _block(object):

    def __init__(self):
        
        self.stations = []
        self.channels = []
        self.inventory = {}
        self.station_pairs = []
        self.channel_pairs = []

    def __repr__(self):

        return "Block containing %g channel pairs" %len(self.channel_pairs)

    def __str__(self):

        return "Block containing %g channel pairs" %len(self.channel_pairs)


class correlation_inventory(object):

    def __init__(self,cfg):

        self.cfg = cfg

        # - Find available data
        self.stations, self.files = file_inventory(cfg)
        all_stations = self.stations.keys()

        # - Determine station pairs
        # - station pairs are grouped into blocks
        self.station_blocks = station_pairs(all_stations,
            cfg.n_stationpairs, cfg.corr_autocorr)


        self.blocks = []

        # - Determine channel pairs for each station pair
        # - station pairs are grouped into blocks
        for block in self.station_blocks:
            self._add_corrblock(block)

    def __repr__(self):
        return "Correlation inventory"

    def __str__(self):
        return "%g blocks in correlation inventory" %len(self.blocks)

    def _add_corrblock(self,station_block):

        block = _block()
        

        # Station pairs and channel pairs should be at the same index.
        block.station_pairs = station_block[:]
        block.channel_pairs = []

        # Make a unique station list for this block
        for p in station_block:
            block.stations.append(p[0])
            block.stations.append(p[1])
        block.stations = list(set(block.stations))

        # Find the relevant channel combinations
        for p in station_block:
            
            sta1 = p[0]
            sta2 = p[1]
            cpairs = channel_pairs(self.stations[sta1],
                    self.stations[sta2],self.cfg)
            block.channel_pairs.append(cpairs)

        # Make a unique channel list for this block
        # for an R or T channel, both horizontal components are included               
        for cp in block.channel_pairs:
            for c in cp:
                for cha in c:
                    if cha[-1] == 'R':
                        block.channels.append(re.sub('R$','N',cha))
                        block.channels.append(re.sub('R$','E',cha))
                    elif cha[-1] == 'T':
                        block.channels.append(re.sub('T$','N',cha))
                        block.channels.append(re.sub('T$','E',cha))
                    else:
                        block.channels.append(cha)
        block.channels = list(set(block.channels))



        # Add file inventory for the channels in question
        inventory = {c: self.files[c] for c in block.channels}
        block.inventory = deepcopy(inventory)

        # No channels are found if updating, and those pairs have already been computed.
        if block.channels != []: 
            self.blocks.append(block)
