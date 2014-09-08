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
from .defs import *
from .scalars import cvxpy_obj
from .arrays import cvxpy_array
from .arrays import cvxpy_matrix
from functools import reduce

#***********************************************************************#
# Class definition: cvxpy_constr                                        #
#***********************************************************************#
class cvxpy_constr(object):

    # Method: __init__
    def __init__(self,left,constraint_type,right):
        """
        Class constructor.

        :param left: Left hand side.
        :param constraint type: Keyword (See cvxpy.defs).
        :param right: Right hand side.
        """

        self.left = left
        self.type = constraint_type
        self.right = right

    # Method: __getattribute__
    def __getattribute__(self,name):

        # Variables
        if name == 'variables':
            if self.type == BELONGS:
                return self.left.variables
            else:
                return cvxpy_list(set(self.left.variables + 
                                      self.right.variables))

        # Parameters
        elif name == 'parameters':
            if self.type == BELONGS:
                return self.left.parameters
            else:
                return cvxpy_list(set(self.left.parameters + 
                                      self.right.parameters))

        # Other
        else:
            return object.__getattribute__(self,name)

    # Method: __str__
    def __str__(self):

        l_text = str(self.left)
        type_text = self.type
        r_text = str(self.right)
        return l_text+' '+type_text+' '+r_text

    # Method: is_dcp
    def is_dcp(self):
        """
        Determines if the constraint is DCP-compliant.
        """

        if (self.type == EQUALS and
            self.left.is_affine() and
            self.right.is_affine()):
            return True
        elif (self.type == LESS_EQUALS and
              self.left.is_convex() and
              self.right.is_concave()):
            return True
        elif (self.type == GREATER_EQUALS and
              self.left.is_concave() and
              self.right.is_convex()):
            return True
        elif (self.type == BELONGS and
              self.left.is_affine()):
            return True
        else:
            return False
        
    # Method: is_affine
    def is_affine(self):
        
        if self.type == BELONGS:
            return False
        else:
            return self.left.is_affine() and self.right.is_affine()
        
#***********************************************************************#
# Class definition: cvxpy_list                                          #
#***********************************************************************#
class cvxpy_list(list):

    # Method: __getattributes__
    def __getattribute__(self,name):

        # Variables
        if name == 'variables':
            all_vars = reduce(lambda x,y: x+y,
                              list(map(lambda x: x.variables,self)),[])
            return cvxpy_list(set(all_vars))

        # Parameters
        elif name == 'parameters':
            all_params = reduce(lambda x,y: x+y,
                                list(map(lambda x: x.parameters,self)),[])
            return cvxpy_list(set(all_params))

        # Other
        else:
            return object.__getattribute__(self,name)
                    
    # Method: _get_eq
    def _get_eq(self):
        """
        Gives equality constraints.
        """
    
        return cvxpy_list([c for c in self if c.type == EQUALS])

    # Method: _get_ineq_in
    def _get_ineq_in(self):
        """
        Gives inequality and membership constraints.
        """
        
        return cvxpy_list([c for c in self if c.type != EQUALS])

    # Method: is_dcp
    def is_dcp(self):
        """
        Determines if all constraints are DCP-compliant.
        """

        return all(list(map(lambda x:x.is_dcp(),self)))

    # Method: is_affine
    def is_affine(self):
        """
        Determines if all constraints are affine.
        """
        
        return all(list(map(lambda x:x.is_affine(),self)))
        
    # Method: __add__
    def __add__(self,other):        
        
        return cvxpy_list(list(self) + other)

    # Method: __radd__
    def __radd__(self,other):
        
        return cvxpy_list(other + list(self))

    # Method: __str__
    def __str__(self):

        output = ''
        for i in range(0,len(self),1):
            output += str(self[i])
            if i != len(self)-1:
                output += '\n'
        return output

#***********************************************************************#
# Function definition: compare                                          #
#***********************************************************************#
def compare(obj1,constraint_type,obj2):
    """
    Compares obj1 with obj2.

    :param obj1: Left hand side obejct.
    :param constraint_type: Keyword (See cvxpy.defs.).
    :param obj2: Right hand side object.
    """
    
    # Both scalars 
    if ((np.isscalar(obj1) or type(obj1).__name__ in SCALAR_OBJS) and
        (np.isscalar(obj2) or type(obj2).__name__ in SCALAR_OBJS)):
        
        # Upgrade scalars to cvxpy_obj
        if np.isscalar(obj1):
            obj1 = cvxpy_obj(CONSTANT,obj1,str(obj1))
        if np.isscalar(obj2):
            obj2 = cvxpy_obj(CONSTANT,obj2,str(obj2))

        # Construct and return constraint
        return cvxpy_constr(obj1,constraint_type,obj2)

    # Upgrate scalars to arrays
    if ((type(obj1) is cvxpy_matrix or type(obj1).__name__ in ARRAY_OBJS) and
        (np.isscalar(obj2) or type(obj2).__name__ in SCALAR_OBJS)):
        (m,n) = obj1.shape
        new_ar = cvxpy_array(m,n)
        for i in range(0,m,1):
            for j in range(0,n,1):
                new_ar[i,j] = obj2
        obj2 = new_ar
    if ((type(obj2) is cvxpy_matrix or type(obj2).__name__ in ARRAY_OBJS) and
        (np.isscalar(obj1) or type(obj1).__name__ in SCALAR_OBJS)):
        (m,n) = obj2.shape
        new_ar = cvxpy_array(m,n)
        for i in range(0,m,1):
            for j in range(0,n,1):
                new_ar[i,j] = obj1
        obj1 = new_ar
    
    # Both arrays
    if ((type(obj1) is cvxpy_matrix or type(obj1).__name__ in ARRAY_OBJS) and
        (type(obj2) is cvxpy_matrix or type(obj2).__name__ in ARRAY_OBJS)):
        constr = []
        if obj1.shape != obj2.shape:
            raise ValueError('Invalid dimensions')
        (m,n) = obj1.shape
        for i in range(0,m,1):
            for j in range(0,n,1):
                constr += [compare(obj1[i,j],constraint_type,obj2[i,j])]
        return cvxpy_list(constr)

    # Invalid arguments
    raise TypeError('Objects not comparable')    
