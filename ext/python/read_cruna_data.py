#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:14:12 2020

@author: arne_hoelter
"""

# import required packages
from   os.path import splitext, isfile
from   glob    import glob
from   h5py    import File
from   pandas  import DataFrame
import numpy   as     np

#%% Read CRUNA Data
# reads image data files created by the cruna code 
# and creates a mono-image data array
#
# (in) fname       : file name of the first image
# (in) component   : (optional) component to be read, default = all
#                    0: density
#                    1: velocity in x-direction
#                    2: velocity in y-direction
#                    3: velocity in z-direction
#                    4: pressure
# (in) images      : (optional) list of images to be read, default = all
#                    must be a list containing the numbers of the images to read
#                    use only to read different samples!!!
# (in) read_params : (optional) list
#                    if first entry is False, second entry must be dict of params
#
# (out) data       : data set of fname (all images)
#                    if more than one component --> [x, y, z, components, times steps]
#                    if only one component      --> [x, y, z, time steps]
# (out) params     : parameter set of fname (given image)
def read_cruna_data(fname, **kwargs):
    # set default values
    if 'component'   not in kwargs.keys(): kwargs['component']   = 'all'
    if 'images'      not in kwargs.keys(): kwargs['images']      = 'all'
    if 'read_params' not in kwargs.keys(): kwargs['read_params'] = [True, False]
    if kwargs['images'] != 'all':
        images = np.array(kwargs['images'])
           
    #%% read params of fname
    if kwargs['read_params'][0]:
        params = read_cruna_params(fname)
    else:
        params = kwargs['read_params'][1]
    #%% determine data names
    # numbers of images
    if not 'images' in locals():
        images = np.arange(1,params['parallelism']['i1']*params['parallelism']['i2'] * \
                           params['parallelism']['i3']+1, dtype=int)
            
    rawfname, ext = splitext(fname) # get file extension and raw filename
    ftype         = rawfname.split('__') # get type and number
    fimage        = ftype[1].split('_') # split images
    
    dnames = []
    for n in range(len(images)):
        block_nbr = fimage[0]
        image_nbr = str(images[n]).zfill(5)
        
        file_idents = ''
        for m in range(2,len(fimage)):
            file_idents = f'{file_idents}_{fimage[m]}'
        dnames.append(f'{ftype[0]}__{block_nbr}_{image_nbr}{file_idents}.h5')
    
    #%% determine gnames (index files)
    gnames = []
    for n in range(len(images)):
        block_nbr = fimage[0]
        image_nbr = str(images[n]).zfill(5)
        gnames.append(f'geometry_indices__{block_nbr}_{image_nbr}.h5')
    
    #%% create data field
    # default layout 3d field
    ns = [int(params['geom']['n1']),int(params['geom']['n2']),int(params['geom']['n3']),1,1]
    
    # get dimension of dataset data
    _ , DatasetDims, _ = give_DatasetInfo(dnames[0], 'data')
    
    # extent data layout corresponding to dataset dimension in file
    for d in range(3,len(DatasetDims)):
        ns[d] =  DatasetDims[d]
    
    # reduce number of dimension if only a component is read
    if kwargs['component'] != 'all':
        ns[3] = 1
        
        
    # if all samples shall be read, overwrite dnames
    # if 'sample' in dnames[0] and kwargs['images'] == 'all':
    #     raise NotImplementedError('to be implemented')
        #dnames = glob('*sample*')
    # init field
    print(f'Estimated variable size: {np.array(ns).prod()*8/(1000**3)} GB')
    if 'sample' in dnames[0]:
        data = np.zeros(list(DatasetDims))
    else:
        data = np.zeros(ns)
    
    #%% read data/params
    for r in range(len(dnames)):
        print(f'read data file: {dnames[r]}')
        
        params_image = read_cruna_params(dnames[r])
        
        _, DatasetDims, DatasetName = give_DatasetInfo(dnames[r], 'data')
        start = np.zeros(len(DatasetDims), dtype=int)
        count = DatasetDims
        
        if len(DatasetDims) <= 3 or 'geom' in dnames[r]:
            data_image = np.array(File(dnames[r], 'r')['data']).transpose()
        elif 'pena' in dnames[r]:
            data_image = np.array(File(dnames[r], 'r')['data']).transpose()
        elif 'sample' in dnames[r]:
            data_image = np.array(File(dnames[r], 'r')['data']).transpose()
        else:
            if kwargs['component'] != 'all':
                start[3] = kwargs['component']
                count[3] = kwargs['component']+1
            if 'snap' in dnames[r]:
                data_image = File(dnames[r], 'r')['data'][start[3]:count[3],:,:,:].transpose()
            else:
                data_image = File(dnames[r], 'r')['data'][:,start[3]:count[3],:,:,:].transpose()
        
        # assignment
        if 'xi10i' in params_image['geom']:
            min_xi = int(params_image['geom']['xi10i']-1)
            max_xi = int(params_image['geom']['xi11i'])
            min_yi = int(params_image['geom']['xi20i']-1)
            max_yi = int(params_image['geom']['xi21i'])
            min_zi = int(params_image['geom']['xi30i']-1)
            max_zi = int(params_image['geom']['xi31i'])
        elif isfile(gnames[r]):
            raise NotImplementedError("gnames to implement")
        else:
            raise NotImplementedError("No global block assignment available: please code")
            
        for ax in range(np.ndim(data_image), 5):
            if 'sample' not in dnames[r]:
                data_image = np.expand_dims(data_image, -1)
            
        if 'sample' in dnames[r]:
            for n in range(np.shape(data_image)[0]):
                if 'pos' in dnames[r]:
                    if np.any(data_image[n,:]):
                        data[n,:] = data_image[n,:]
                else:
                    if np.any(data_image[n,:,:]):#before data_image[n,4,:]
                        data[n,:,:] = data_image[n,:,:]
        else:
            data[min_xi:max_xi, min_yi:max_yi, min_zi:max_zi,:,:] = data_image
    
    data = np.squeeze(data)
    return data, params


#%% read cruna parameter
# (in) fname     : file name of the first image
# (out) params   : dictionary with parameters
def read_cruna_params(fname):
    rawfname, ext = splitext(fname) # get file extension
    params = dict() # create dict of keys and values (output)
    if ext=='.h5': # if data is stored in a hdf5 file
        paramsraw = np.array(DataFrame(File(fname, 'r')['params']))
    else: # assume ASCII file
        paramsraw = list()
        with open(fname) as params_obj:
            for line in params_obj:
                #line=line.replace(' ', '')  
            # fid = open(fname)
            # row=fid.readline().rstrip()
                if not line:
                    pass
                elif line=='\n':
                    pass
                elif np.any(line[0]=='%'):
                    pass
                else:
                    line=line.replace(' ', '')   
                    if np.any(line[0]=='%'):
                        pass
                    else:
                        paramsraw.append(line.replace('\n','%').split('%')[0])
                        
        paramsraw = np.array(DataFrame(paramsraw))
        
    for n in range(np.shape(paramsraw)[0]):
        # split key and value
        if ext=='.h5':
            key_value_pair = replace_all({"[":"", "]":"", "'":"", " ":"", "\\t":""}, \
                                         str(paramsraw[n]))[1:].split('=')
        else:
            key_value_pair = replace_all({"[":"", "]":"", "'":"", " ":"", "\\t":""}, \
                                         str(paramsraw[n]))[0:].split('=')
        # check if key string is empty
        if not key_value_pair[0]:
            continue
        # # store parameter name in key list
        # key2val[key_value_pair[0]] = ''
        # set value
        if np.shape(key_value_pair)[0] > 1:
            try:
                value = [float(m) for m in list(filter(None, key_value_pair[1].split(';')))]
            except:
                value = list(filter(None, key_value_pair[1].split(';')))
            
        if len(value) == 1:
            value = value[0]

        # get struct levels
        splitkeys = key_value_pair[0].split('.')
        levels = len(splitkeys)
        
          
        # write keys and values into params dict
        if levels == 1:
            params[splitkeys[0]] = value
        elif levels == 2:
            if splitkeys[0] not in params:
                params[splitkeys[0]] = dict()
            params[splitkeys[0]][splitkeys[1]] = value
        elif levels == 3:
            if splitkeys[0] not in params:
                params[splitkeys[0]] = dict()
            if splitkeys[1] not in params[splitkeys[0]]:
                params[splitkeys[0]][splitkeys[1]] = dict()
            # elif splitkeys[0] in params and if splitkeys[1] not in params[splitkeys[0]]:
            #     params[splitkeys[0]][splitkeys[1]] = dict()
            params[splitkeys[0]][splitkeys[1]][splitkeys[2]] = value
        else:
            break
            
    return params

#%% returns Dataset information of a given h5-file and DatasetName
def give_DatasetInfo(fname, DatasetName):    
    h5file = File(fname, 'r')
    for ind, key in enumerate(h5file.keys()):
        if key == DatasetName:
            DatasetIndex = ind
            DatasetName  = key
    DatasetDims = np.flip(np.shape(h5file[DatasetName]))
    return DatasetIndex, DatasetDims, DatasetName

   
#%% function to loop the replace-function
def replace_all(dict, str):
    for key in dict:
        str = str.replace(key, dict[key])
    return str
