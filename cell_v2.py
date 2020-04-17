# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 13:26:22 2019

@author: Dr. Ayodeji Babalola
"""
import numpy as np
#----------------------------------------------------------------------------
class cell():
    def __init__(self,nr,nc=None,dim =None):
        self.nr  = nr
        self.nc  = nc
        self.dim  = dim
        self.mat = None
        self.mat_out = None
        self.type = None
        self.info = None
        self.init()
#----------------------------        
    def init(self):
        if (self.nr is not None and self.nc is None and self.dim is None):
            self.type = '1D'
        elif(self.nc is not None and self.dim is None):
            self.type = '2D'
        elif (self.nr is not None and self.nc is not None and self.dim is None):
            self.type = '3D'
            
        if(self.nc is None):
            self.nc = 1
            self.dim = None
            self.mat=   [[] for _ in range(self.nr)]
        elif(self.nc is not None and self.dim is None):
            #self.mat = [ [[] for _ in range(self.nc)] for _ in range(self.nr)] # prefered
            self.mat = [[[] for i in range(self.nc)] for j in range(self.nr)]
        elif(self.dim is not None):
            self.mat = [[[ [] for col in range(self.nr)] for col in range(self.nc)] for row in range(self.dim)]
#----------------------------            
    def __repr__(self):
        return repr(self.mat)
#----------------------------      
    def __setitem__(self,ind,value):
        if (type(ind) == int):
            nr = ind
            if(nr > self.nr):
                raise Exception ("number of rows is out of bounds") 
            else:
                self.mat[nr] = value       
            
        elif(len(ind) == 2):
            nr = ind[0]
            nc = ind[1]
            if(nr > self.nr):
                raise Exception ("number of rows is out of bounds")
            elif(nc>self.nc):
                raise Exception ("number of columns is out of bounds")
            else:
                self.mat[nr][nc] = value
            
        elif(len(ind) == 3):
            nr = ind[0]
            nc = ind[1]
            dim = ind[2]
            if(nr > self.nr):
                raise Exception ("number of rows is out of bounds")
            elif(nc>self.nc):
                raise Exception ("number of columns is out of bounds")
            elif(dim>self.dim):
                raise Exception ("number of dimensions is out of bounds")
            else:
                self.mat[nr][nc][dim] = value     
        
        return self.mat_out
#----------------------------    
    def __getitem__(self,ind):
        if (type(ind) == int):
            nr = ind
            if(nr > self.nr):
                raise Exception ("number of rows is out of bounds") 
            else:
                mat_out = self.mat[nr]         
            
        elif(len(ind) == 2):
            nr = ind[0]
            nc = ind[1]
            if(nr > self.nr):
                raise Exception ("number of rows is out of bounds")
            elif(nc>self.nc):
                raise Exception ("number of columns is out of bounds")
            else:
                mat_out = self.mat[nr][nc]
            
        elif(len(ind) == 3):
            nr = ind[0]
            nc = ind[1]
            dim = ind[2]
            if(nr > self.nr):
                raise Exception ("number of rows is out of bounds")
            elif(nc>self.nc):
                raise Exception ("number of columns is out of bounds")
            else:
                mat_out = self.mat[nr][nc][dim]   
        return mat_out
              
        def __delitem__(self, ind):
            if (type(ind) == int):
                nr = ind
                if(nr > self.nr):
                    raise Exception ("number of rows is out of bounds") 
                else:
                   del self.mat[nr]         
            
            elif(len(ind) == 2):
                nr = ind[0]
                nc = ind[1]
                if(nr > self.nr):
                    raise Exception ("number of rows is out of bounds")
                elif(nc>self.nc):
                    raise Exception ("number of columns is out of bounds")
                else:
                    del self.mat[nr][nc]
            
            elif(len(ind) == 3):
                nr = ind[0]
                nc = ind[1]
                dim = ind[2]
                if(nr > self.nr):
                    raise Exception ("number of rows is out of bounds")
                elif(nc>self.nc):
                    raise Exception ("number of columns is out of bounds")
                else:
                    del self.mat[dim][nc][dim]
#----------------------------                    
    def size(self,dim = None):
        if (self.type == '1D'):
            sz = np.array([self.nr,1])        
        if (self.type == '2D'):
            sz = np.array([self.nr,self.nc])
        elif(self.type == '3D'):
            sz = np.array([self.nr,self.nc,self.dim])   
                
        return sz
    """
    #----------------------------                    
    def size(self,dim = None):
        if (self.type == '1D'):
            sz = np.array([self.nr,1,1])        
        if (self.type == '2D'):
            sz = np.array([self.nr,self.nc,1])
        elif(self.type == '3D'):
            sz = np.array([self.nr,self.nc,self.dim])   
                
        return sz[0],sz[1],sz[2]
    """
#----------------------------    
    def size_data(self): 
# mimicking cellfub=n@size
# works only for 1d arrays 
# buggy
        if (self.nc is None and self.nc is None):
            sz = np.array(self.nr)
            for i in range(self.nr):
                sz[i] = self.mat[i].size
        elif(self.nc is not None and self.dim is None):
            sz = np.zeros((self.nr,self.nc))
            for ii in range(self.nc):
                for i in range(self.nr):
                    tmp = self.mat[i][ii]
                    if (tmp == []):
                        sz[i,ii] = 0
                    else:
                        sz[i,ii] = tmp.size
                    
        elif(self.dim is not None):
            sz = np.zeros(((self.nr,self.nc,self.dim)))
            for iii in range(self.dim):
                tmp2d = self.mat[iii]  # slice of 2D
                for ii in range(self.nc):
                    tmp1d  = tmp2d[ii]  # extract col
                    for i in range(self.nr): 
                        tmp = tmp1d[i]   # extract row
                        sz[i,ii] = tmp.size
        return sz
#----------------------------    
    def extr_allrows(self,col):
        tmp  = np.array([])
        for i in range(self.nr):
            if (i==0):
                tmp = self.mat[i][col]  
            else:
                tmp = np.concatenate((tmp,self.mat[i][col]),axis = None)
        return tmp
        
#----------------------------    
    def extr_allcols(self,rows):
        tmp  = np.array([])
        for i in range(self.nc):
            if (i==0):
                tmp = self.mat[rows][i]
            else:
                tmp = np.concatenate((tmp,self.mat[rows][i]),axis  = None)   
        return tmp
#----------------------------    
    def extr_allrows_cell(self,col):
        tmp  = cell(self.nr)
        for i in range(self.nr):
            tmp[i] = self.mat[i][col]

        return tmp
        
#----------------------------    
    def extr_allcols_cell(self,rows):
        tmp  = cell(self.nc)
        for i in range(self.nc):
            tmp[i] = self.mat[rows][i]

        return tmp    
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:49:07 2019

@author: Ayo
"""

