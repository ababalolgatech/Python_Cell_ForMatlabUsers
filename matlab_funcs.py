"""
Created on Mon Aug 19 17:26:09 2019
@author: Dr. Ayodeji Babalola
"""
import matplotlib.pyplot as plt
import pickle 
import json
import h5py
import pandas as pd
import numpy as np
#import hickle as hkl
from scipy import spatial
import segyio as segyio
import Filter as ft
import os
import cell
cell = getattr(cell,"cell")
from progressbar import Percentage,Bar,ProgressBar
import sincinterpol as sincinterpol
import crewes_resample as cresamp
import math

# https://deepdish.readthedocs.io/en/latest/io.html
# To do : cell arrays not to throw errors with nr>nc indexing 
#------------------------------------------------------------------------------
def numpy_save(data,fp): 
    fp = fp + '.npy'
    np.save(fp,data) 
#------------------------------------------------------------------------------
def numpy_load(fp):
   fp = fp + '.npy'
   file_out = np.load(fp)
   return file_out
#------------------------------------------------------------------------------
def numpy_map(fp,flag = None,ind = None):
    fp = fp + '.npy'
    dat_map = np.load(fp, mmap_mode='r')
    
    if (flag is None):
        flag = 'allcols'
    if (ind is None):
        ind = 0
    
    if (flag == 'allrows'):
        dat = dat_map[:,ind]
    elif(flag  == 'allcols'):
       dat = dat_map[ind,:] 
    return dat
#------------------------------------------------------------------------------
def create_obj(obj_in,n):
   obj = [obj_in() for i in range(n)] 
   return obj  
""" 
#------------------------------------------------------------------------------
# https://pypi.org/project/hickle/
def save_obj(fname,obj):
    fname = fname + '.hkl'
    hkl.dump(obj,fname,mode='w')
#-----------------------------------------------------------------------------        
def load_obj(fname):
    fname = fname + '.hkl'
    obj=hkl.load(fname)
    return 
"""    
#-----------------------------------------------------------------------------  
def save_obj(fname,obj):
    fname = fname + '.obj'
    with open(fname, "wb") as fp:
        pickle.dump(obj, fp, pickle.HIGHEST_PROTOCOL)
        fp.close()
#-----------------------------------------------------------------------------        
def load_obj(fname):
    fname = fname + '.obj'
    with open(fname, "rb") as fp:
        obj = pickle.load(fp)
        fp.close()
    return obj
#-----------------------------------------------------------------------------  
def save_json(fname,obj):
    fname = fname+'.json'
    with open(fname, "w") as fp:
        json.dump(obj,fp) 
        fp.close()
#-----------------------------------------------------------------------------        
def load_json(fname):
    fname = fname +'.json'
    f = open(fname)
    obj = json.load(f)
    
    return obj
#-----------------------------------------------------------------------------  
def save_h5py(fname,obj,Type = None):
    # https://stackoverflow.com/questions/44049838/h5py-write-object-dynamically-to-file
    fname = fname+'.h5py'
    with h5py.File(fname, 'w') as f:
        for item in vars(obj).items():
            f.create_dataset(item[0], data = item[1])
#-----------------------------------------------------------------------------  
def save_h5py2(fname,obj):
    # https://www.pythonforthelab.com/blog/how-to-use-hdf5-files-in-python/
    fname = fname+'.h5py'
    with h5py.File(fname, "w") as fp:
        dset  = fp.create_dataset('default',data = obj)
        fp.close()        
#-----------------------------------------------------------------------------        
def load_h5py(fname,Type = None):
    if (Type is None):
        fname = fname +'.json'
    else:
        fname = fname + Type
        
    fp = open(fname)
    obj = json.load(fp)
    fp.close()
    return obj   
#-----------------------------------------------------------------------------        
def mat_load(fname,data,corr_indx = None):
    fname = fname+'.mat'
    f = h5py.File(fname,'r')
    data = f.get(data)
    data = np.array(data)
    if (corr_indx is None):
        data = data.T
    else:
        data = data.T[:,corr_indx:-1]
    return data
         
#-----------------------------------------------------------------------------
def load_ascii(fp,nheader=None,field =None):
    """
    if (nheader is not None and nheader !=0):
        nheader = nheader-1    
    """
    if (nheader is None):
        nheader = 0
        
    df = pd.read_csv(fp,skiprows =nheader,delimiter= '\s+', index_col=False)    
    matrix = df.values
    

        
    if (field is None):
        mat = matrix
    else:
        if (type(field) == np.ndarray) :
            nfield = field.size
        else:
            nfield = len(field)
        
        if(nfield != matrix.shape[1]):
            print('Number of rows is not equal to the field')
            ncol = matrix.shape[1]
            index = np.arange(ncol)
            field = np.array([])
            for ii in range(ncol):
                field[ii] = index[ii+1]
            mat = matrix
        else:       
            mat = {}        
            for i in range (0,nfield):
                mat[field[i]] = matrix[:,i]        
 
    return mat
#----------------------------------------------------------------------------
def cell2mat(dat):
    nr = dat.nr
    nc = dat.nc
    ndim = dat.dim
    tmp = np.array([])
    if (dat.type == '1D'):
        if (nr>nc):
            for i in range(nr):
                if (i==0):
                    tmp = dat[i]
                    mat = np.zeros((nr,tmp.size))     # BUG.. concatenating matrix   
                    mat[i,:] = tmp 
                else:                             
                    mat[i,:] = tmp
        elif(nc>nr):
            for i in range(nr): # BUG
                tmp = np.hstack((tmp,dat[i]))
                if (i==0):
                    mat = np.zeros((tmp.size,nc))        
                    mat[:,i] = tmp    
        elif(nc==1 and nr==1):          
            mat = dat[0]                     
    elif(dat.type == '2D'): 
        """       
        for ii in range(nc):
            for i in range (nr):
                tmp = np.concatenate((tmp,dat[i,ii]),axis = None)
            if (ii==0):
                mat = np.zeros((tmp.size,nc))        
            mat[:,ii] = tmp
        """
        nr,nc = dat.size
        for ii in range(nr):
            tmp = []
            for i in range(nc): 
                if (i==0):
                    tmp = dat[ii,i]
                else:
               # mat[(self.nsamp*ii) + ind: (self.nsamp*ii)+ind] = wavdata[ii,i] 
                   tmp = np.hstack((tmp,dat[ii,i]))
            if(ii == 0):
                mat = tmp
            else:
                mat = np.vstack((mat,tmp))
    
        tmp = np.array([])   
    return  mat
#----------------------------------------------------------------------------
def cell2cube(dat):
    print(' BUG in cell3D')
    nr,nc = dat.size()
    ndim = dat.dim
    mat = np.zeros((nr,nc,ndim))
    
    for i in range(ndim):
        mat[:,:,i] = cell2mat(dat[:,:,i])
  
    return  mat
         
 #----------------------------------------------------------------------------
def matlab_any(arr):
    if (arr.size == 0):
        out = False
    else:
        tmp = np.asarray(arr.nonzero())
        tmp = tmp.sum()
        
        if(tmp ==0):
            out = True
        else:
            out = False    
    return out
#-----------------------------------------------------------------
def var(arr,axis_=0):
    out = np.var(arr,axis = axis_,ddof=1,dtype=np.float64)
    # ddof is degree of freedom  
    return out
#-----------------------------------------------------------------
def std(arr,axis_=0):
    out = np.std(arr,axis = axis_,ddof=1,dtype=np.float64)
    # ddof is degree of freedom     
    return out
#-----------------------------------------------------------------
def corr(arr1,arr2):
    mtx = np.corrcoef(arr1,arr2)
    out = mtx[0,1]
    # ddof is degree of freedom     
    return out
#-----------------------------------------------------------------
def find(arr,val):
    if (np.isnan(val) == True):
        tmp = np.argwhere(np.isnan(arr))
        ind = tmp.astype(int)
        ind = np.asarray(ind)
    else:        
        tmp = np.nonzero(arr == val)
        #ind = np.int(tmp[0])
        #ind = np.transpose(ind)
        ind = np.asarray(tmp)
        ind = np.ravel(ind)
        
        if (ind.size ==1):
            ind = ind[0]           
    return ind

#-----------------------------------------------------------------
def findnot(arr,val):
    if (np.isnan(val) == True):
        tmp = np.argwhere(~np.isnan(arr))
        ind = tmp.astype(int)
        ind = np.asarray(ind)
    else:        
        tmp = np.nonzero(arr != val)
        #ind = np.int(tmp[0])
        #ind = np.transpose(ind)
        ind = np.asarray(tmp)
        ind = np.ravel(ind)
        
        if (ind.size ==1):
            ind = ind[0]            
    return ind 
#-----------------------------------------------------------------   
def findarr(arr,val):
    if (val.size==0):
        nval=1
        val = np.array([val])
    else:
        nval = val.size
        
    ind = np.zeros(nval)
    
    for i in range(nval):
        ind[i] = find(arr,val[i])
    return ind.astype(int)
#-----------------------------------------------------------------
def gt(arr,val):
    tmp = np.nonzero(arr > val)
        #ind = np.int(tmp[0])
        #ind = np.transpose(ind)
    ind = np.asarray(tmp)
    ind = np.ravel(ind)        
    if (ind.size ==1):
        ind = ind[0]            
    return ind
#-----------------------------------------------------------------
def lt(arr,val):        
    tmp = np.nonzero(arr < val)
        #ind = np.int(tmp[0])
        #ind = np.transpose(ind)
    ind = np.asarray(tmp)
    ind = np.ravel(ind)        
    if (ind.size ==1):
        ind = ind[0]            
    return ind
#-----------------------------------------------------------------
def ge(arr,val):
    tmp = np.nonzero(arr >= val)
        #ind = np.int(tmp[0])
        #ind = np.transpose(ind)
    ind = np.asarray(tmp)
    ind = np.ravel(ind)        
    if (ind.size ==1):
        ind = ind[0]            
    return ind
#-----------------------------------------------------------------
def le(arr,val):        
    tmp = np.nonzero(arr <= val)
        #ind = np.int(tmp[0])
        #ind = np.transpose(ind)
    ind = np.asarray(tmp)
    ind = np.ravel(ind)        
    if (ind.size ==1):
        ind = ind[0]            
    return ind
#-----------------------------------------------------------------
def minn(arr):        
   val = np.min(arr)            
   return int(val)
#-----------------------------------------------------------------
def maxx(arr):        
   val = np.max(arr)            
   return int(val)
#-----------------------------------------------------------------
def replace(dat,numm,replc):
    ind = find(dat,np.nan)
    dat[ind] = replc
    return dat

#-----------------------------------------------------------------
def mround(dat):
    ind1 = dat % 1 >= 0.5
    ind2 = np.logical_not(ind1)
    data = np.zeros(dat.size)
    data[ind1] = np.ceil(dat[ind1])  
    data[ind2] = np.round(dat[ind2])
    return data
#-----------------------------------------------------------------
def seis_resamp(seis_in,dtold,dtnew,options = None):
    # this resampling codes increases the sample rates linearly
    # if ns_old is 101 at dt=4, ns_new will be 404
    if (options is None):
        options = 'default'
        
    nr = seis_in.nr
    nc = seis_in.nc            
    seis_out = cell(nr,nc)
    scal = np.int(dtold/dtnew)
    if (options == 'default'):
        for ii in range(nc):
            for i in range(nr):
                seis_out[i,ii] = sincinterpol.upsample3(seis_in[i,ii],scal)            
    else:
        for ii in range(nc):
            for i in range(nr):
                seis_out[i,ii] = ft.resamp_by_interp(seis_in[i,ii],dtold,dtnew)           
            
    return seis_out
#-----------------------------------------------------------------
def seis_resamp2(seis_in,dtold,dtnew,options = None):
    # this resampling codes deos not increases the sample rates linearly
    # if ns_old is 101 at dt=4, ns_new will be 401
    
    if (options is None):
        options = 'sinc'
        
    nr = seis_in.nr
    nc = seis_in.nc            
    seis_out = cell(nr,nc)
    if(seis_in.type == '1D'):
      old_nsamp =seis_in[0].size   
    else:
        old_nsamp =seis_in[0,0].size 
    tmax = old_nsamp*dtold - dtold    
    new_nsamp =  np.int(dtold/dtnew) * (old_nsamp-1) + 1
    old_tseis = np.linspace(0,tmax,old_nsamp,dtype = int)
    new_tseis = np.linspace(0,tmax,new_nsamp,dtype = int)
     
    if (options == 'sinc'):  # better
        for ii in range(nc):
            for i in range(nr):
                seis_out[i,ii] = ft.sinc_interp(seis_in[i,ii],old_tseis,new_tseis)
    elif(options == 'crewes'):
        for ii in range(nc):
            for i in range(nr):
                seis_out[i,ii],tmp = cresamp.resamp_crewes(seis_in[i,ii],old_tseis,dtnew)     
    else:
        for ii in range(nc):
            for i in range(nr):
                seis_out[i,ii] = ft.scipy_interp(seis_in[i,ii],old_tseis,new_tseis)         
            
    return seis_out
#-----------------------------------------------------------------
def m_size(arr):    
    if(type(arr) == tuple):
        sz = (1,1)   
    else:
        try:
            sz =( arr.shape[0],arr.shape[1])
        except :
            try:
                sz =( 1,arr.shape[1])
            except:
               try :  
                   sz =( arr.shape[0],1)
               except:
                   try:
                       tmp = arr.size
                       if (tmp==1):
                           sz = (1,1)
                       else:
                           sz=(tmp,1)
                   except :
                        pass
    return sz
#-----------------------------------------------------------------
def shape(arr):
    if (arr.ndim == 1):
        nr = arr.size
        nc = 1
    else:
        nr,nc = arr.shape
    return nr,nc
#----------------------------         
def process_arr(dat):
    # check for all data and return an array as matrix
    nr,nc = m_size(dat)
    
    if (nr == 1 or nc==1):
        mat = np.zeros((nr,nc))
        if (nr >1 and nc==1):
            mat[:,0] = dat
        elif (nc>1 and nr==1):
            mat[0,:] = dat
        else:
            raise Exception('BUG: Inside matlab_func-> process_arr')
    else:
        mat = dat

    return mat

#----------------------------         
def convert2_arr(dat):
    # convvert 1d matrix to array
    nr,nc = m_size(dat)
    if (nr == 1 or nc==1):
        mat = dat.flatten()  
    return mat
#-----------------------------------------------------------------    
def padd(gth,n):
   pf = np.zeros(n)
   pb = np.zeros(n)
   pf[0:n] = gth[0]
   pb[0:n] = gth[-1]
   
   pdd = np.concatenate((pf,gth,pb),axis = None)
   return pdd
#-----------------------------------------------------------------
def upadd(gth,n):
   a = gth.shape[0]
   pdd = gth[n:(a-n)]
   return pdd  
   
#-----------------------------------------------------------------
def kde_tree(arr,loc):
   # kdtree = spatial.cKDTree(arr,k=1)
    kdtree = spatial.cKDTree(arr)
    if (loc.size ==2):
         dist,ind = kdtree.query(loc,k=1)
    else:
        nr,nc = loc.shape
        dist = np.zeros(nr)
        ind = np.zeros(nr,dtype = int)
        for i in range(nr):
            dist[i],ind[i] = kdtree.query(loc[i,:],k=1)

    return ind,dist
#---------------------------- 
def progress(miniter,maxiter,msg = None):   
    itr = round((miniter/maxiter)*100)
    rem = np.remainder(itr,10)
    if (rem == 0):
        if (msg is None):
            print (str(itr) + '%')
        else:
            print( msg + str (itr) + '%')
        
    if (itr == maxiter):
        print('done')

#------------------------------------------------------------------------------   
def seglog(zones,z_in,Vp,Vs,Rho):
    zones = np.fix(zones)
    aa = find(z_in, zones[0])  # np.nonzero(z_in == zones[0])
    bb = find(z_in, zones[1])
    if (aa.size==0):
        raise Exception ('lower zone exceeded')
    elif(bb.size==0):
        raise Exception ('upper zone exceeded') 
    z_out = z_in[aa:bb+1]  
    vp = Vp[aa:bb+1]
    vs = Vs[aa:bb+1] 
    rho = Rho[aa:bb+1] 
        
    return vp,vs,rho,z_out
#------------------------------------------------------------------------------
def segdat(zones,tseis,seis):
    zones = np.fix(zones)
    aa = find(tseis, zones[0])  # np.nonzero(z_in == zones[0])
    bb = find(tseis, zones[1])    
    if (aa.size==0):
        raise Exception ('lower zone exceeded')
    elif(bb.size==0):
        raise Exception ('upper zone exceeded')  
        
    time = tseis[aa:bb] 
    
    if (isinstance(seis,cell) == False):
        raise Exception ('seismic inpute must be cell arrays')
    
    ndat = seis.nr
    nang = seis.nc
    data = cell(ndat,nang)
    
    for ii in range(nang):
        for i in range(ndat):
            tmp = seis[i,ii]
            data[i,ii] = tmp[aa:bb+1]
    return time,data

#------------------------------------------------------------------------------
def segdat_trc(zones,tseis,seis):
    zones = np.fix(zones)
    aa = find(tseis, zones[0])  # np.nonzero(z_in == zones[0])
    bb = find(tseis, zones[1])    
    if (aa.size==0):
        raise Exception ('lower zone exceeded')
    elif(bb.size==0):
        raise Exception ('upper zone exceeded')   
    time = tseis[aa:bb+1] 
    data = seis[aa:bb+1]
    return time,data

#------------------------------------------------------------------------------
def tlcord(time,dt,vp1,vs1,rho1):
# function to recoordinate log sampling with seismic
# should be applied prior to resampling with seismic sampling rate    
    if(time[1]-time[0] == dt):
        tcord = time
        vp = vp1
        vs = vs1
        rho = rho1
    else:
        minT = np.amin(time)
        maxT = np.amax(time)
        num1 = minT - (dt - np.round(np.remainder(minT,dt)))
        num2 = maxT - (dt - np.round(np.remainder(maxT,dt)))
        
        t1 = find(time,num1)
        t2 = find(time,num2)
        vp = vp1[t1:t2]
        vs = vs1[t1:t2]
        rho = rho1[t1:t2]
        tcord = time[t1:t2]
    return rho,vs,vp,tcord
        
#------------------------------------------------------------------------------
def segment_logs(zones,z_in,dat_in):
    zones = np.round(zones)
    aa = find(z_in, zones[0])  # np.nonzero(z_in == zones[0])
    bb = find(z_in, zones[1])
    if (aa.size==0):
        raise Exception ('lower zone exceeded')
    elif(bb.size==0):
        raise Exception ('upper zone exceeded')        
    dat_out = dat_in[aa:bb+1]
    z_out = z_in[aa:bb+1]        
    return dat_out,z_out
#----------------------------
def reverse_array(arr):
    dat = arr[::-1]
    return dat
#-----------------------------------------------------------------------------
def rms (sig):
    rms = np.sqrt(np.mean(sig**2))
    return rms
#-----------------------------------------------------------------------------
def file_size (fp,Type = None):
    if (Type is None):
        Type = 'numpy'
    
    if (Type == 'numpy'):
        fp = fp + '.npy'
    else:
        fp = fp + '.obj'
        
    statinfo = os.stat(fp)
    sz = statinfo.st_size * 1e-6
    return sz

#-----------------------------------------------------------------------------
def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]
#-----------------------------------------------------------------------------
def remove_nan_interp(dat):
    """
    remove np.nan values by interpolation
    """
    nans, x= nan_helper(dat)
    dat[nans]= np.interp(x(nans), x(~nans), dat[~nans])
    return dat
#-----------------------------------------------------------------------------
def remove_zero_interp(dat):
    nans, x= nan_helper(dat)
    dat[nans]= np.interp(x(nans), x(~nans), dat[~nans])
    return dat
#-----------------------------------------------------------------------------    
def remove_nan_topbase(time_in, dat_in):
    """
    remove np.nan values from top and base of data passed in
    """
    testnan = find(dat_in,np.nan)
    if (testnan.size != 0):            
       #   indnan = np.argwhere(np.isnan(dat_in))      
       ns = dat_in.size
       indnan = np.isnan(dat_in)   
       indx1 = np.nonzero(indnan==0)
       indx2 = np.nonzero(indnan==1)
    #        indx1 = np.argwhere(np.isnan(dat_in))  
            
       if (indx1 is None):
            time_out = time_in
            dat_out = dat_in
       else:
            min_zone = np.min(indx2)
            max_zone = np.max(indx2) 
            if(max_zone +1 < ns): # no nan value at the base
                min_zone = max_zone + 1
                max_zone = ns
            elif (max_zone+1 ==ns and min_zone==0):
                min_zone = np.min(indx1)
                max_zone = np.max(indx1) + 1  # : loses one index
            else:
                min_zone = np.min(indx1)
                max_zone = np.max(indx1)
       time_out = time_in[min_zone:max_zone]  # : loses one index
       dat_out = dat_in[min_zone:max_zone]
    else:
        dat_out = dat_in
        time_out = time_in
        
    return time_out,dat_out    
#-----------------------------------------------------------------------------
def plotd(dat_in,time=None)  :
    fig, (ax) = plt.subplots(1,1)
    if (time is None):
        ax.plot(dat_in)
    else:
        ax.plot(dat_in,time)
    plt.show()
    plt.pause(1)
    return ax
#-----------------------------------------------------------------------------
def scatterd(dat1,dat2=None)  :
    plt.figure()
    if (dat2 is None):
        raise Exception('input second data')
    else:
        plt.scatter(dat1,dat2)
    plt.show()
    plt.pause(1)
#-----------------------------------------------------------------------------
def imshow(dat_in,time=None)  :
    plt.figure()
    if (time is None):
        #plt.imshow(dat_in,interpolation='bicubic',cmap = 'nipy_spectral', aspect='auto')
        plt.imshow(dat_in,interpolation='nearest',cmap = 'nipy_spectral', aspect='auto')
        plt.colorbar()
    else:
        plt.imshow(dat_in,interpolation='bicubic',cmap = 'nipy_spectral', aspect='auto')
        # plt.imshow(dat_in,interpolation='bicubic',cmap = 'gist_ncar', aspect='auto')
        plt.colorbar()
    #plt.show()
    #plt.pause(1)
    
#-----------------------------------------------------------------------------
def imshowd(dat_in,time=None)  :
    plt.figure()
    if (time is None):
        plt.imshow(dat_in,interpolation='bicubic',cmap = 'nipy_spectral', aspect='auto')
        plt.colorbar()
    else:
        #plt.imshow(dat_in,interpolation='bicubic',cmap = 'nipy_spectral', aspect='auto')
        plt.imshow(dat_in,interpolation='nearest',cmap = 'nipy_spectral', aspect='auto')
        # plt.imshow(dat_in,interpolation='bicubic',cmap = 'gist_ncar', aspect='auto')
        plt.colorbar()
    plt.show()
    plt.pause(1)    
 #-----------------------------------------------------------------------------
def plot(dat_in,time=None)  :
    fig, (ax) = plt.subplots(1,1)
    if (time is None):
        ax.plot(dat_in)
    else:
        ax.plot(dat_in,time)
    plt.show()
    return ax
#-----------------------------------------------------------------------------
def histd(dat_in,time=None)  :
    plt.figure()
    plt.hist(dat_in)
    plt.show()
    plt.pause(1)   
    #-------------------------------    
def get_indx(indx): 
    if (type(indx) ==  int) :
        n =1
        indx = np.array([indx])
    elif (type(indx) ==  str) :
        n =1
        indx = np.array([indx]) 
    elif (type(indx) ==  np.int32) :
        n =1
        indx = np.array([indx])         
    else:        
        n = indx.size    
    return indx,n
#----------------------------         
def read_numpy_seismic(fp,indx,flag = None):
    if (flag is None):
       flag = 'allcols'
    indx,nXL = get_indx(indx)
    fp,nseismic = get_indx(fp)    
  
    mat_cell = cell(nseismic)   
    for ii in range(nseismic):
        nfp = fp[ii] + '.npy'  
        tmp = cell(nXL)
        dat_map =  np.load(nfp, mmap_mode='r') 
        if (flag == 'allcols'):
            tmp = dat_map[indx,:] # each trace in a column
        else:
             tmp = dat_map[:,indx]             
        mat_cell[ii] = tmp
    
    if (nseismic >1):    
        aa,ns = tmp.shape
        cube = np.zeros((nXL,ns,nseismic))
        
        for i in range(nseismic):
            cube[:,:,i] = mat_cell[i]
    else:
        cube = tmp
        
       
    return cube
#----------------------------         
def read_numpy_seismic_cell(fp,indx,flag = None):
    if (flag is None):
       flag = 'allcols'    
    indx,nXL = get_indx(indx)
    fp,nseismic = get_indx(fp)

    mat_cell = cell(nXL,nseismic)   
    for ii in range(nseismic):
        nfp = fp[ii] + '.npy'
        dat_map =  np.load(nfp, mmap_mode='r') 
        if (flag == 'allcols'):
            tmp = dat_map[indx,:] # each trace in a column
        else:
             tmp = dat_map[:,indx]
        for i in range(nXL):
              mat_cell[i,ii] = tmp[i,:]
        tmp = []
       
    return mat_cell
 
#----------------------------         
def read_seismic(fp,options):
    #https://github.com/equinor/segyio/blob/master/python/segyio/tracefield.py
     SDataCoord = {}
     if options == 'petrel':
         with segyio.open(fp, ignore_geometry=True) as segyfile:
            segyfile.mmap()
            SDataCoord['Lines'] = segyfile.attributes(segyio.TraceField.TRACE_SEQUENCE_FILE)[:] # byte 5 - Petrel & SMT
            SDataCoord['Traces'] = segyfile.attributes(segyio.TraceField.CDP)[:]
            SDataCoord['X'] = segyfile.attributes(segyio.TraceField.SourceX)[:]
            SDataCoord['Y'] = segyfile.attributes(segyio.TraceField.SourceY)[:]
            SDataCoord['start_time'] = segyfile.attributes(segyio.TraceField.DelayRecordingTime)[:]
            SData = segyfile.trace.raw[:]
     elif options =='fp':
        with segyio.open(fp, ignore_geometry=True) as segyfile:
            segyfile.mmap()
            SDataCoord['Lines'] = segyfile.attributes(segyio.TraceField.FieldRecord)[:] # byte 9 -FP
            SDataCoord['Traces'] = segyfile.attributes(segyio.TraceField.CDP)[:]
            SDataCoord['X'] = segyfile.attributes(segyio.TraceField.SourceX)[:]
            SDataCoord['Y'] = segyfile.attributes(segyio.TraceField.SourceY)[:]
            SDataCoord['start_time'] = segyfile.attributes(segyio.TraceField.DelayRecordingTime)[:]
            SData = segyfile.trace.raw[:]
     return SData, SDataCoord       
 
#----------------------------     
def write_segy(fp,SDataFile,SDataCoordFile,dt = None, start_time = None,options = None):
# https://segyio.readthedocs.io/en/latest/segyio.html?highlight=segyio%20su#segyio.su.words.sut
# The link shows the byte locations
# make sure you load with trace-headers only
# time(dt and start_time) must be in milliseconds
    
    txt  = { ' 1. Make sure you load with trace-headers only,\
               2.time(dt and start_time) must be in milliseconds,\
               inline = byte 9, crossline = byte 21, X = byte 73, Y,byte 77,\
               sart_time = byte 109, dt = byte 115\
               Seygio modified by Dr. Ayodeji Babalola' }
               
# time(dt and start_time) must be in milliseconds}
    if (options is None):
        options = 'loaded'
    if (dt is None):
        dt = 1
    if (start_time is None):
        start_time  = 0
    
    SDataCoord = load_obj(SDataCoordFile)
    
    ns = SDataCoord['Lines'].size
    Lines = SDataCoord['Lines']
    Traces = SDataCoord['Traces']
    line_nos = np.unique(Lines)
    trc_nos  = np.unique(Traces)
    dlines = line_nos[1] -line_nos[0]
    dtraces = trc_nos[1]  -trc_nos[0]
    
    Inlines = np.array([np.min(line_nos), np.max(line_nos)+dlines])
    Xlines = np.array([np.min(trc_nos), np.max(trc_nos)+dtraces])
    
   
    trc = read_numpy_seismic (SDataFile,0)
    samples = trc.size
    
    
    spec = segyio.spec()
    spec.sorting = 2
    spec.format = 1
    spec.samples = range(int(samples))
    spec.ilines = range(*map(int, Inlines))
    spec.xlines = range(*map(int, Xlines))
   
    
    #spec.tracecount=ilines[1]-ilines[0]
    spec.tracecount = Inlines[1] - Inlines[0]    
    pbar = ProgressBar(maxval = ns).start()

    if (options == 'loaded'): 
        seis = numpy_load(SDataFile)
        # Add flags for ns later
        with segyio.create(fp, spec) as f:
            for i in range(ns):                
                pbar.update(i)
                #f.text[0] = txt  do it later
                f.header[i] = {                    
                    segyio.su.offset : 1,
                    segyio.su.fldr  : SDataCoord['Lines'][i].astype(int),   # byte 9
                    segyio.su.cdp   : SDataCoord['Traces'][i].astype(int),  # byte 21
                    segyio.su.sx    : SDataCoord['X'][i].astype(int),       # byte 73
                    segyio.su.sy    : SDataCoord['Y'][i].astype(int),       # byte 77
                    segyio.su.dt    : dt,                            # byte 115                  
                    segyio.su.delrt : start_time,                    # byte 109
                    segyio.su.ns    : samples                        # byte 115
                }
                f.trace[i] = seis[i,:]                       
            f.bin.update(tsort=segyio.TraceSortingFormat.INLINE_SORTING) 
    else:
        with segyio.create(fp, spec) as f:
            for i in range(ns):
                pbar.update(i)
                #f.text = txt
                f.header[i] = {                
                    segyio.su.fldr  : SDataCoord['Lines'][i].astype(int),   # byte 9
                    segyio.su.cdp   : SDataCoord['Traces'][i].astype(int),  # byte 21
                    segyio.su.sx    : SDataCoord['X'][i].astype(int),       # byte 73
                    segyio.su.sy    : SDataCoord['Y'][i].astype(int),       # byte 77
                    segyio.su.dt    : dt,                       # byte 115                  
                    segyio.su.delrt : start_time,               # byte 109
                    segyio.su.ns    : samples                   # byte 115 
                }
                f.trace[i] = np.array(read_numpy_seismic(SDataFile,i))    
            f.bin.update(tsort=segyio.TraceSortingFormat.INLINE_SORTING)         
        
    pbar.finish()


#----------------------------         
def hist(dat,bins):
    #data = dat[:-1] + np.gradient(dat)/2
    data = dat + np.gradient(dat)/2
    aa,bb = np.histogram(data,bins)
    aa = np.ones(data.shape)
    return aa,bb
 
#---------------------------------------            
def calc_dist(x1,y1,x2,y2):        
    distance = np.sqrt((y2 - y1)**2 + (x2 -x1)**2) 
    if(math.isnan(distance)):
        distance = 0
    return distance
#----------------------------         
def rand(num):
    return  np.random.uniform(0,1,num)
#----------------------------         
def bxfun_add(a,b):
    n_b = b.size    
    tmp = np.tile(a,(n_b,1))
    tmp = tmp.T
    
    for i in range(n_b):
        tmp[:,i] = tmp[:,i] + b[i]
        
    return tmp
#---------------------------
def bxfun_minus(a,b):     
    n_b = b.size    
    tmp = np.tile(a,(n_b,1))
    tmp = tmp.T
    
    for i in range(n_b):
        tmp[:,i] = tmp[:,i] - b[i]
        
    return tmp

#---------------------------
def bxfun_mult(a,b,): 
    """
    a mult b
    """
    n_b = b.size    
    tmp = np.tile(a,(n_b,1))
    tmp = tmp.T
    
    for i in range(n_b):
        tmp[:,i] = tmp[:,i] * b[i]
        
    return tmp

#---------------------------
def bxfun_div(a,b,):     
    n_b = b.size    
    tmp = np.tile(a,(n_b,1))
    tmp = tmp.T
    
    for i in range(n_b):
        tmp[:,i] = tmp[:,i] / b[i]
        
    return tmp  
#---------------------------
def bxfun_hypot(a,b,): 
    n_a = a.size    
    n_b = b.size    
    tmp = np.tile(a,(n_b,1))
    tmp = tmp.T
    tmp2 = np.zeros((n_a,n_b))
    
    for i in range(n_b):
        tmp2[:,i] = np.hypot(tmp[:,i],b[i])
        
    return tmp2         
#******************************************************************************
if __name__ == '__main__':
    xx = np.array([1,1.5,2,2.5])
    xx_out = mround(xx)
    """
    SDataCoordFile = 'models\\rand_update\\LF\\seis_coord'
    SDataFile = 'models\\rand_update\\LF\\AI3D_mean'
    fp_out = 'AI3D_mean_rand_update_lf.sgy'
    write_segy(fp_out,SDataFile,SDataCoordFile,1000,2148)
    
    SDataCoordFile = 'models\\rand_update\\LF\\seis_coord'
    SDataFile = 'models\\rand_update\\LF\\AI3D_rand'
    fp_out = 'AI3D_rand_rand_update_LF.sgy'
    write_segy(fp_out,SDataFile,SDataCoordFile,1000,2148)    
    """
    """
    SDataCoordFile = 'models\\mean_update\\LF\\seis_coord'
    SDataFile = 'models\\mean_update\\LF\\AI3D_rand'
    fp_out = 'AI3D_rand_mean_update_LF.sgy'
    write_segy(fp_out,SDataFile,SDataCoordFile,1000,2148)
    
    SDataCoordFile = 'models\\mean_update\\HF\\seis_coord'
    SDataFile = 'models\\mean_update\\HF\\AI3D_mean'
    fp_out = 'AI3D_mean_mean_update_HF.sgy'
    write_segy(fp_out,SDataFile,SDataCoordFile,1000,2148)  
    SDataFile = 'models\\mean_update\\HF\\AI3D_rand'
    fp_out = 'AI3D_rand_mean_update_HF.sgy'
    write_segy(fp_out,SDataFile,SDataCoordFile,1000,2148)
    """
    
    
    SDataCoordFile = 'models\\seis_coord'
    SDataFile = 'models\\AI3D_mean'
    fp_out = 'AI3D_mean_dec6.sgy'
    write_segy(fp_out,SDataFile,SDataCoordFile,1000,2148)  
    SDataFile = 'models\\AI3D_rand'
    fp_out = 'AI3D_rand_dec6.sgy'
    write_segy(fp_out,SDataFile,SDataCoordFile,1000,2148)    
