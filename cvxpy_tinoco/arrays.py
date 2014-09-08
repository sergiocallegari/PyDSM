#***********************************************************************#
# Copyright (C) 2010-2012 Tomas Tinoco De Rubira                        #
#                                                                       #
# This file is part of CVXPY                                            #     
#                                                                       #
# CVXPY is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #   
# (at your option) any later version.                                   # 
#                                                                       #
# CVXPY is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program.  If not, see <http://www.gnu.org/licenses/>. #
#***********************************************************************#

import numpy as np
from cvxpy.defs import *
import scipy.sparse as sp

#***********************************************************************#
# Class definition: cvxpy_array                                         #
#***********************************************************************#
class cvxpy_array(object):

    # Method: __init__
    def __init__(self,m,n):
        """
        Class constructor.

        :param m: Number of rows.
        :param n: Number of columns.
        """

        self.shape = (m,n)
        temp1 = {}
        for i in range(0,m,1):
            temp2 = {}
            for j in range(0,n,1):
                temp2[j] = 0.0
            temp1[i] = temp2
        self.data = temp1
        self.type = ARRAY
    
    # Method:__getattribute__
    def __getattribute__(self,name):

        # Transpose
        if name == 'T':
            new_ar = cvxpy_array(self.shape[1],self.shape[0])
            for i in range(0,self.shape[0],1):
                for j in range(0,self.shape[1],1):
                    new_ar.data[j][i] = self.data[i][j]
            return new_ar

        # Variables
        elif name == 'variables':
            l = []
            for i in range(0,self.shape[0],1):
                for j in range(0,self.shape[1],1):
                    if not np.isscalar(self.data[i][j]):
                        l += self.data[i][j].variables
            return cvxpy_list(set(l))

        # Parameters
        elif name == 'parameters':
            l = []
            for i in range(0,self.shape[0],1):
                for j in range(0,self.shape[1],1):
                    if not np.isscalar(self.data[i][j]):
                        l += self.data[i][j].parameters
            return cvxpy_list(set(l))

        # Value
        elif name == 'value':
            mat = cvxpy_matrix(np.zeros((self.shape[0],self.shape[1])),np.float64)
            for i in range(0,self.shape[0],1):
                for j in range(0,self.shape[1],1):
                    if np.isscalar(self.data[i][j]):
                        mat[i,j] = self.data[i][j]
                    else:
                        mat[i,j] = self.data[i][j].value
            return mat
        
        # Other attribute
        else:
            return object.__getattribute__(self,name)
    
    # Method: __setitem__
    def __setitem__(self,key,value):

        if type(key) != tuple:
            raise TypeError('Invalid Key')
        if len(key) != 2:
            raise TypeError('Invalid Key')
        if (not (issubclass(type(key[0]),int) or 
                 issubclass(type(key[0]),np.integer)) or
            not (issubclass(type(key[1]),int) or 
                 issubclass(type(key[1]),np.integer))):
            raise TypeError('Invalid Key')
        if key[0] < 0 or key[0] >= self.shape[0]:
            raise ValueError('Index out of range')
        if key[1] < 0 or key[1] >= self.shape[1]:
            raise ValueError('Index out of range')

        self.data[key[0]][key[1]] = value

    # Method: __getitem__
    def __getitem__(self,key):

        if type(key) != tuple:
            raise TypeError('Invalid Key')
        if len(key) != 2:
            raise TypeError('Invalid Key')
        if (not (issubclass(type(key[0]),int) or 
                 issubclass(type(key[0]),np.integer)) and
            type(key[0]) is not slice):
            raise TypeError('Invalid Key')
        if (not (issubclass(type(key[1]),int) or 
                 issubclass(type(key[1]),np.integer)) and
            type(key[1]) is not slice):
            raise TypeError('Invalid Key')

        # Prepare left slice
        sl = [None,None,None]
        if np.isscalar(key[0]):
            sl = [key[0],key[0]+type(key[0])(1),1]
        else:
            if key[0].start == None:
                sl[0] = 0
            else:
                sl[0] = key[0].start
            if key[0].stop == None:
                sl[1] = self.shape[0]
            else:
                sl[1] = key[0].stop
            if key[0].step == None:
                sl[2] = 1
            else:
                sl[2] = key[0].step

        # Prepare right slices
        sr = [None,None,None]
        if np.isscalar(key[1]):
            sr = [key[1],key[1]+type(key[1])(1),1]
        else:
            if key[1].start == None:
                sr[0] = 0
            else:
                sr[0] = key[1].start
            if key[1].stop == None:
                sr[1] = self.shape[1]
            else:
                sr[1] = key[1].stop
            if key[1].step == None:
                sr[2] = 1
            else:
                sr[2] = key[1].step
                    
        # Check type
        if not all(list(map(lambda x: (issubclass(type(x),int) or 
                                       issubclass(type(x),np.integer)),
                            [sl[0],sl[1],sl[2],sr[0],sr[1],sr[2]]))):
            raise TypeError('Invalid key')            

        # Calculate new size            
        ll = range(sl[0],sl[1],sl[2])
        lr = range(sr[0],sr[1],sr[2])
        new_m = len(ll)
        new_n = len(lr)
        
        # Scalar
        if new_m == 1 and new_n == 1:
            if sl[0] < 0 or sl[0] >= self.shape[0]:
                raise ValueError('Index out of range')
            if sr[0] < 0 or sr[0] >= self.shape[1]:
                raise ValueError('Index out of range')
            return self.data[sl[0]][sr[0]]

        # Array
        new_ar = cvxpy_array(new_m,new_n)
        for i in range(0,new_m,1):
            
            if ll[i] < 0 or ll[i] >= self.shape[0]:
                raise ValueError('Index out of range')
            
            for j in range(0,new_n,1):
                
                if lr[j] < 0 or lr[j] >= self.shape[1]:
                    raise ValueError('Index out of range')
                
                new_ar.data[i][j] = self.data[ll[i]][lr[j]]
        return new_ar

    # Method: __str__
    def __str__(self):

        output = '['
        for i in range(0,self.shape[0],1):
            if i:
                output += '\n ['
            else:
                output += '['
            for j in range(0,self.shape[1],1):
                output += ' ' + str(self.data[i][j]) + ' '
            output += ']'
        output += '] '
        return output

    # Method: __add__
    def __add__(self,other):

        return self.__addsub__(other,LEFT_ADD)

    # Method: __sub__
    def __sub__(self,other):

        return self.__addsub__(other,LEFT_SUBTRACT)

    # Method: __radd__
    def __radd__(self,other):

        return self.__addsub__(other,RIGHT_ADD)

    # Method: __rsub__
    def __rsub__(self,other):

        return self.__addsub__(other,RIGHT_SUBTRACT)

    # Method: __addsub__
    def __addsub__(self,other,action):
        """
        Handles addition and subraction.

        :param other: Other operand.
        :param action: Keyword (See cvxpy.defs).
        """
        
        # New array
        (m,n) = self.shape
        new_ar = cvxpy_array(m,n)

        # Scalar
        if (np.isscalar(other) or 
            type(other).__name__ in SCALAR_OBJS):
            for i in range(0,m,1):
                for j in range(0,n,1):
                    if action == LEFT_ADD:
                        new_ar.data[i][j] = self.data[i][j] + other
                    elif action == LEFT_SUBTRACT:
                        new_ar.data[i][j] = self.data[i][j] - other
                    elif action == RIGHT_ADD:
                        new_ar.data[i][j] = other + self.data[i][j]
                    else:
                        new_ar.data[i][j] = other - self.data[i][j]
            return new_ar

        # Array or matrix
        elif (type(other).__name__ in MATRIX_OBJS or 
              type(other).__name__ in ARRAY_OBJS):
            if other.shape != self.shape:
                raise ValueError('Invalid Dimensions')
            else:
                for i in range(0,m,1):
                    for j in range(0,n,1):
                        if action == LEFT_ADD:
                            new_ar.data[i][j] = self.data[i][j] + other[i,j]
                        elif action == LEFT_SUBTRACT:
                            new_ar.data[i][j] = self.data[i][j] - other[i,j]
                        elif action == RIGHT_ADD:
                            new_ar.data[i][j] = other[i,j] + self.data[i][j]
                        else:
                            new_ar.data[i][j] = other[i,j] - self.data[i][j]
                return new_ar

        # Sparray or spmatrix
        elif (type(other) is cvxpy_spmatrix or 
              type(other) is cvxpy_sparray):
            if other.shape != self.shape:
                raise ValueError('Invalid Dimensions')
            else:
                for i in range(0,m):
                    for j in range(0,n):
                        new_ar.data[i][j] = self.data[i][j]
                for i in range(0,m,1):
                    rowi_indeces = other.rows[i]
                    rowi_values = other.data[i]
                    for k in range(0,len(rowi_indeces)):
                        j = rowi_indeces[k]
                        value = rowi_values[k]
                        if action == LEFT_ADD:
                            new_ar.data[i][j] = new_ar.data[i][j] + value
                        elif action == LEFT_SUBTRACT:
                            new_ar.data[i][j] = new_ar.data[i][j] - value
                        elif action == RIGHT_ADD:
                            new_ar.data[i][j] = value + new_ar.data[i][j]
                        else:
                            new_ar.data[i][j] = value - new_ar.data[i][j]
                return new_ar

        # Not implemented
        else:
            return NotImplemented
    
    # Method: __mul__
    def __mul__(self,other):

        return self.__mulhandle__(other,LEFT_MULTIPLY)
    
    # Method: __rmul__
    def __rmul__(self,other):

        return self.__mulhandle__(other,RIGHT_MULTIPLY)
   
    # Method: __mulhandle__
    def __mulhandle__(self,other,action):
        """
        Handles multiplication.

        :param other: Other operand.
        :param action: Keyword (See cvxpy.defs).
        """

        # Dimensions
        (m,n) = self.shape

        # Scalar
        if (np.isscalar(other) or 
            type(other).__name__ in SCALAR_OBJS):
            new_ar = cvxpy_array(m,n)
            for i in range(0,m):
                for j in range(0,n):
                    if action == LEFT_MULTIPLY:
                        new_ar.data[i][j] = self.data[i][j] * other
                    else:
                        new_ar.data[i][j] = other * self.data[i][j]
            return new_ar

        # Array or matrix
        elif (type(other).__name__ in MATRIX_OBJS or 
              type(other).__name__ in ARRAY_OBJS):
            (p,q) = other.shape
            if action == LEFT_MULTIPLY:
                if p != n:
                    raise ValueError('Invalid dimensions')
                else:
                    new_ar = cvxpy_array(m,q)
                    for i in range(0,m):
                        for j in range(0,q):
                            temp = 0
                            for k in range(0,n):
                                temp += self.data[i][k]*other[k,j]
                            new_ar.data[i][j] = temp

                    # Convert to scalar if shape is (1,1)
                    if new_ar.shape == (1,1):
                        return new_ar.data[0][0]
                    else:
                        return new_ar
            else:
                if q != m:
                    raise ValueError('Invalid dimensions')
                else:
                    new_ar = cvxpy_array(p,n)
                    for i in range(0,p):
                        for j in range(0,n):
                            temp = 0
                            for k in range(0,m):
                                temp += other[i,k]*self.data[k][j]
                            new_ar.data[i][j] = temp

                    # Convert to scalar if shape = (1,1)
                    if new_ar.shape == (1,1):
                        return new_ar.data[0][0]
                    else:
                        return new_ar

        # Sparray or spmatrix
        elif (type(other) is cvxpy_spmatrix or 
              type(other) is cvxpy_sparray):
            (p,q) = other.shape
            if action == LEFT_MULTIPLY:
                if p != n:
                    raise ValueError('Invalid dimensions')
                else:
                    new_ar = cvxpy_array(m,q)
                    other._create_col_based_rep()
                    for j in range(0,q):
                        colj_indeces = other.cols[j]
                        colj_values = other.datac[j]
                        for i in range(0,m):
                            temp = 0.
                            for r in range(0,len(colj_indeces)):
                                k = colj_indeces[r]
                                temp += self.data[i][k]*colj_values[r]
                            new_ar.data[i][j] = temp    

                    # Convert to scalar if shape is (1,1)
                    if new_ar.shape == (1,1):
                        return new_ar.data[0][0]
                    else:
                        return new_ar
            else:
                if q != m:
                    raise ValueError('Invalid dimensions')
                else:
                    new_ar = cvxpy_array(p,n)
                    for i in range(0,p):
                        rowi_indeces = other.rows[i]
                        rowi_values = other.data[i]
                        for j in range(0,n):
                            temp = 0.
                            for r in range(0,len(rowi_indeces)):
                                k = rowi_indeces[r]
                                temp += rowi_values[r]*self.data[k][j]
                            new_ar.data[i][j] = temp

                    # Convert to scalar if shape = (1,1)
                    if new_ar.shape == (1,1):
                        return new_ar.data[0][0]
                    else:
                        return new_ar
        else:
            return NotImplemented
    
    # Method: __neg__
    def __neg__(self):

        # Allocate new array
        (m,n) = self.shape
        new_ar = cvxpy_array(m,n)

        # Fill and return array
        for i in range(0,m,1):
            for j in range(0,n,1):
                new_ar.data[i][j] = -self.data[i][j]
        return new_ar
    
    # Method: is_affine
    def is_affine(self):
        """
        Determines if the array contains affine
        expressions only.
        """
        
        for i in range(0,self.shape[0],1):
            for j in range(0,self.shape[1],1):
                if (not np.isscalar(self.data[i][j]) and
                    not self.data[i][j].is_affine()):
                    return False
        return True

#***********************************************************************#
# Class definition: cvxpy_var                                           #
#***********************************************************************#
class cvxpy_var(cvxpy_array):

    # Method: __init__
    def __init__(self,m,n,structure=None,name=None):
        """
        Class constructor.

        :param name: String.
        :param m: Number of rows.
        :param n: Number of columns.
        :param structure: Keyword (See cvxpy.defs).
        """
        
        # Set name
        if name == None:
            name = BASE_VAR+str(cvxpy_scalar_var.i)
            cvxpy_scalar_var.i += 1
            
        # Call parent constructor
        cvxpy_array.__init__(self,m,n)
        
        # No structure
        if structure == None:
            for i in range(0,m,1):
                for j in range(0,n,1):
                    v = cvxpy_scalar_var(name+'['+str(i)+','+str(j)+']')
                    self.data[i][j] = v
        
        # Lower triangular
        elif structure == LOWER_TRIANGULAR:
            if m!=n:
                raise ValueError('Invalid dimensions')
            for i in range(0,m,1):
                for j in range(0,i+1,1):
                    v = cvxpy_scalar_var(name+'['+str(i)+','+str(j)+']')
                    self.data[i][j] = v

        # Upper triangular
        elif structure == UPPER_TRIANGULAR:
            if m!=n:
                raise ValueError('Invalid dimensions')
            for i in range(0,m,1):
                for j in range(i,n,1):
                    v = cvxpy_scalar_var(name+'['+str(i)+','+str(j)+']')
                    self.data[i][j] = v
        
        # Symmetric
        elif structure == SYMMETRIC:
            if m!=n:
                raise ValueError('Invalid dimensions')
            for i in range(0,m,1):
                for j in range(0,i+1,1):
                    v = cvxpy_scalar_var(name+'['+str(i)+','+str(j)+']')
                    self.data[i][j] = v
                    self.data[j][i] = v
        
        # Error
        else:
            raise ValueError('Invalid structure')

#***********************************************************************#
# Class definition: cvxpy_param                                         #
#***********************************************************************#
class cvxpy_param(cvxpy_array):

    # Method: __init__
    def __init__(self,m,n,attribute=None,name=None):
        """
        Class constructor.

        :param name: String.
        :param m: Number of rows.
        :param n: Number of columns.
        :param attribute: Keyword (See cvxpy.defs).
        """

        # Set name
        if name == None:
            name = BASE_PARAM+str(cvxpy_scalar_param.i)
            cvxpy_scalar_param.i += 1

        # Call parent constructor
        cvxpy_array.__init__(self,m,n)
        
        # Fill in
        for i in range(0,m,1):
            for j in range(0,n,1):
                v = cvxpy_scalar_param(attribute,
                                       name+'['+str(i)+','+str(j)+']')
                self.data[i][j] = v

    # Method: __setattr__
    def __setattr__(self,name,value):

        # Value
        if name == 'value':
               
            # Check type
            if type(value) is not cvxpy_matrix:
                raise TypeError('Invalid value type')
    
            # Check dimensions
            if value.shape != self.shape:
                raise ValueError('Invalid value dimensions')

            # Store new values
            for i in range(0,self.shape[0],1):
                for j in range(0,self.shape[1],1):
                    self.data[i][j].value = value[i,j]
        
        # Other
        else:
            object.__setattr__(self,name,value)

#***********************************************************************#
# Class definition: cvxpy_matrix                                        #
#***********************************************************************#
class cvxpy_matrix(np.matrix):

    # Method: __init__
    def __init__(self,data,dtype=None,copy=True):
        """
        Class constructor.

        :param data: Array-like or string.
        :param dtype: Data type.
        :param copy: Boolean.

        See numpy.matrix documentation.
        """

        # Parent constructor is called automatically
        pass


    # Method:__getattr__
    def __getattribute__(self,name):

        # Inverse
        if name == 'I':
            temp = np.array(self.copy())
            temp = np.array(np.matrix(temp).I)
            return cvxpy_matrix(temp,np.dtype(np.float64))

        # Other attribute
        else:
            return np.matrix.__getattribute__(self,name)


    # Method: __add__
    def __add__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.matrix.__add__(self,other)

    # Method: __radd__
    def __radd__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.matrix.__radd__(self,other)

    # Method: __sub__
    def __sub__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.matrix.__sub__(self,other)

    # Method: __rsub__
    def __rsub__(self,other):

        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.matrix.__rsub__(self,other)

    # Method: __mul__
    def __mul__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.matrix.__mul__(self,other)
        
    # Method: __rmul__
    def __rmul__(self,other):

        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.matrix.__rmul__(self,other)

#***********************************************************************#
# Class definition: cvxpy_spmatrix                                      #
#***********************************************************************#
class cvxpy_spmatrix(sp.lil_matrix):

    # Method: __init__
    def __init__(self,arg,dtype=None,copy=False):
        """
        Class constructor.

        :param arg: Dense matrix, sparse matrix or shape tuple.
        :param dtype: Data type.
        :param copy: Boolean.

        See scipy.sparse.lil_matrix documentation.
        """
        
        # Call parent constructor
        sp.lil_matrix.__init__(self,arg,dtype=dtype)
        self.format = sp.lil_matrix.__name__[:3]

    # Method: _create_col_based_rep
    def _create_col_based_rep(self):
        self.cols = {}
        self.datac = {}
        for j in range(0,self.shape[1]):
            self.cols[j] = []
            self.datac[j] = []
        for i in range(0,self.shape[0]):
            rowi_indeces = self.rows[i]
            rowi_values = self.data[i]
            for k in range(0,len(rowi_indeces)):
                j = rowi_indeces[k]
                value = rowi_values[k]
                self.cols[j].append(i)
                self.datac[j].append(value)

    # Method:__getattr__
    def __getattribute__(self,name):

        # Inverse
        if name == 'T':
            return cvxpy_spmatrix(self.transpose())
        
        # Other attribute
        else:
            return sp.lil_matrix.__getattribute__(self,name)

    # Method: __add__
    def __add__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return sp.lil_matrix.__add__(self,other)

    # Method: __radd__
    def __radd__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.lil_matrix.__radd__(self,other)

    # Method: __sub__
    def __sub__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.lil_matrix.__sub__(self,other)

    # Method: __rsub__
    def __rsub__(self,other):

        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.lil_matrix.__rsub__(self,other)

    # Method: __mul__
    def __mul__(self,other):
        
        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.lil_matrix.__mul__(self,other)
        
    # Method: __rmul__
    def __rmul__(self,other):

        if (type(other).__name__ in SCALAR_OBJS or
            type(other).__name__ in ARRAY_OBJS):
            return NotImplemented
        else:
            return np.lil_matrix.__rmul__(self,other)

#***********************************************************************#
# Class definition: cvxpy_sparray                                       #
#***********************************************************************#
class cvxpy_sparray(object):

    # Method: __init__
    def __init__(self,m,n):
        """
        Class constructor.

        :param m: Number of rows.
        :param n: Number of columns.
        """

        self.shape = (m,n)
        self.type = SPARRAY
        self.rows = {}
        self.data = {}
        for i in range(0,m):
            self.rows[i] = []
            self.data[i] = []
        self.nnz = 0
        
    # Method: _create_col_based_rep
    def _create_col_based_rep(self):
        self.cols = {}
        self.datac = {}
        for j in range(0,self.shape[1]):
            self.cols[j] = []
            self.datac[j] = []
        for i in range(0,self.shape[0]):
            rowi_indeces = self.rows[i]
            rowi_values = self.data[i]
            for k in range(0,len(rowi_indeces)):
                j = rowi_indeces[k]
                value = rowi_values[k]
                self.cols[j].append(i)
                self.datac[j].append(value)
        
    # Method:__getattribute__
    def __getattribute__(self,name):

        # Transpose
        if name == 'T':
            new_ar = cvxpy_sparray(self.shape[1],self.shape[0])
            for i in range(0,self.shape[0],1):
                rowi_indeces = self.rows[i]
                rowi_values = self.data[i]
                for k in range(0,len(rowi_indeces)):
                    new_ar[rowi_indeces[k],i] = rowi_values[k]
            return new_ar

        # Variables
        elif name == 'variables':
            l = []
            for i in range(0,self.shape[0],1):
                for obj in self.data[i]:
                    if not np.isscalar(obj):
                        l += obj.variables
            return cvxpy_list(set(l))

        # Parameters
        elif name == 'parameters':
            l = []
            for i in range(0,self.shape[0],1):
                for obj in self.data[i]:
                    if not np.isscalar(obj):
                        l += obj.parameters
            return cvxpy_list(set(l))

        # Value
        elif name == 'value':
            mat = cvxpy_spmatrix((self.shape[0],self.shape[1]),np.float64)
            for i in range(0,self.shape[0],1):
                rowi_indeces = self.rows[i]
                rowi_values = self.data[i]
                for k in range(0,len(rowi_indeces)):
                    if np.isscalar(rowi_values[k]):
                        mat[i,rowi_indeces[k]] = rowi_values[k] 
                    else:
                        mat[i,rowi_indeces[k]] = rowi_values[k].value
            return mat
        
        # Other attribute
        else:
            return object.__getattribute__(self,name)

    # Method: __setitem__
    def __setitem__(self,key,value):

        if type(key) != tuple:
            raise TypeError('Invalid Key')
        if len(key) != 2:
            raise TypeError('Invalid Key')
        if (not (issubclass(type(key[0]),int) or 
                 issubclass(type(key[0]),np.integer)) or
            not (issubclass(type(key[1]),int) or 
                 issubclass(type(key[1]),np.integer))):
            raise TypeError('Invalid Key')
        if key[0] < 0 or key[0] >= self.shape[0]:
            raise ValueError('Index out of range')
        if key[1] < 0 or key[1] >= self.shape[1]:
            raise ValueError('Index out of range')
            
        i = key[0]
        j = key[1]
        rowi_indeces = self.rows[i]
        rowi_values = self.data[i]
        
        found = False
        k = 0
        while k < len(rowi_indeces):
            
            if j > rowi_indeces[k]:
                k += 1
                    
            elif j == rowi_indeces[k]:
                found = True
                break;
            
            else:
                break
                    
        # present
        if found:
            if value == 0.:
                self.rows[i] = rowi_indeces[:k]+rowi_indeces[k+1:]
                self.data[i] = rowi_values[:k]+rowi_values[k+1:]
                self.nnz -= 1
            else:
                rowi_values[k] = value
                
        # not present
        else:
            if value != 0.:
                rowi_indeces.insert(k,j)
                rowi_values.insert(k,value)
                self.nnz += 1

    # Method: __getitem__
    def __getitem__(self,key):

        if type(key) != tuple:
            raise TypeError('Invalid Key')
        if len(key) != 2:
            raise TypeError('Invalid Key')
        if (not (issubclass(type(key[0]),int) or 
                 issubclass(type(key[0]),np.integer)) and
            type(key[0]) is not slice):
            raise TypeError('Invalid Key')
        if (not (issubclass(type(key[1]),int) or 
                 issubclass(type(key[1]),np.integer)) and
            type(key[1]) is not slice):
            raise TypeError('Invalid Key')

        # Prepare left slice
        sl = [None,None,None]
        if np.isscalar(key[0]):
            sl = [key[0],key[0]+type(key[0])(1),1]
        else:
            if key[0].start == None:
                sl[0] = 0
            else:
                sl[0] = key[0].start
            if key[0].stop == None:
                sl[1] = self.shape[0]
            else:
                sl[1] = key[0].stop
            if key[0].step == None:
                sl[2] = 1
            else:
                sl[2] = key[0].step

        # Prepare right slices
        sr = [None,None,None]
        if np.isscalar(key[1]):
            sr = [key[1],key[1]+type(key[1])(1),1]
        else:
            if key[1].start == None:
                sr[0] = 0
            else:
                sr[0] = key[1].start
            if key[1].stop == None:
                sr[1] = self.shape[1]
            else:
                sr[1] = key[1].stop
            if key[1].step == None:
                sr[2] = 1
            else:
                sr[2] = key[1].step
                    
        # Check type
        if not all(list(map(lambda x: (issubclass(type(x),int) or 
                                       issubclass(type(x),np.integer)),
                            [sl[0],sl[1],sl[2],sr[0],sr[1],sr[2]]))):
            raise TypeError('Invalid key')            

        # Calculate new size            
        ll = range(sl[0],sl[1],sl[2])
        lr = range(sr[0],sr[1],sr[2])
        new_m = len(ll)
        new_n = len(lr)
        
        # Scalar
        if new_m == 1 and new_n == 1:
            if sl[0] < 0 or sl[0] >= self.shape[0]:
                raise ValueError('Index out of range')
            if sr[0] < 0 or sr[0] >= self.shape[1]:
                raise ValueError('Index out of range')
            
            i = sl[0]
            j = sr[0]
            rowi_indeces = self.rows[i]
            rowi_values = self.data[i]
            
            for k in range(0,len(rowi_indeces)):
                if j == rowi_indeces[k]:
                    return rowi_values[k]
            return np.float64(0.)

        # Array
        new_ar = cvxpy_sparray(new_m,new_n)
        for i in range(0,new_m,1):

            lli = ll[i]
            
            if lli < 0 or lli >= self.shape[0]:
                raise ValueError('Index out of range')
                
            rowi_indeces = self.rows[lli]
            rowi_values = self.data[lli]                
                
            for j in range(0,new_n,1):

                lrj = lr[j]

                if lrj < 0 or lrj >= self.shape[1]:
                    raise ValueError('Index out of range')
        
                for k in range(0,len(rowi_indeces)):
                    if lrj == rowi_indeces[k]:
                        new_ar[i,j] = rowi_values[k]
        return new_ar

    # Method: __str__
    def __str__(self):

        output = ''
        for i in range(0,self.shape[0],1):
            rowi_indeces = self.rows[i]
            rowi_values = self.data[i]
            for k in range(0,len(rowi_indeces)):
                output += ' (' + str(i)+','+str(rowi_indeces[k])+')\t\t'
                output += str(rowi_values[k]) + '\n'
        return output[:-1]
        
# Load modules
from cvxpy.scalars import cvxpy_scalar_var
from cvxpy.scalars import cvxpy_scalar_param
from cvxpy.constraints import cvxpy_list
