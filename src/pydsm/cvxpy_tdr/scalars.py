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
from functools import reduce

#***********************************************************************#
# Class definition: cvxpy_obj                                           #
#***********************************************************************#
class cvxpy_obj(object):

    # Method: __init__
    def __init__(self,object_type,value,name):
        """
        Class constructor.

        :param object_type: Keyword (See cvxpy.def).
        :param value: Number.
        :param name: String.
        """

        self.type = object_type
        self.value = value
        self.name = name
        self.shape = (1,1)
        self.variables = cvxpy_list()
        self.parameters = cvxpy_list()
        self.T = self

    def __lt__(self, other):
        return self.name < other.name

    # Method: is_convex
    def is_convex(self):
        """
        Determines if object follows DCP rules of convexity.
        """

        return True

    # Method: is_concave
    def is_concave(self):
        """
        Determines if object follows DCP rules of concavity.
        """

        return True

    # Method: is_dcp
    def is_dcp(self):
        """
        Determines if object is DCP-compliant.
        """

        return True

    # Method: is_affine
    def is_affine(self):
        """
        Determines if object is affine.
        """

        return True

    # Method: is_nonnegative_constant
    def is_nonnegative_constant(self):
        """
        Determines if object is a nonnegative constant.
        """

        return self.value >= 0

    # Method: is_nonpositive_constant
    def is_nonpositive_constant(self):
        """
        Determines if object is a nonpositive constant.
        """

        return self.value <= 0

    # Method: __add__
    def __add__(self,other):

        return self.__addsub__(other,SUMMATION)

    # Method: __sub__
    def __sub__(self,other):

        return self.__addsub__(other,SUBTRACTION)

    # Method: __radd__
    def __radd__(self,other):
        return self.__raddsub__(other,SUMMATION)

    # Method: __rsub__
    def __rsub__(self,other):

        return self.__raddsub__(other,SUBTRACTION)

    # Method: __addsub__
    def __addsub__(self,other,op):
        """
        Handles left add and left subtract.

        :param other: Right hand side.
        :param op: Keyword (See cvxpy.defs)
        """

        # Create operator
        operator = cvxpy_obj(OPERATOR,np.nan,SUMMATION)

        # Create args to combine
        if (type(self) is cvxpy_tree and
            self.item.name == SUMMATION):
            args = self.children
        else:
            args = [self]

        # Number
        if np.isscalar(other):
            if other == 0.0:
                return self
            if op == SUMMATION:
                constant = cvxpy_obj(CONSTANT,other,str(other))
            else:
                constant = cvxpy_obj(CONSTANT,-other,str(-other))
            return cvxpy_tree(operator,args+[constant])

        # Scalar variable or scalar param
        elif (type(other) is cvxpy_scalar_var or
              type(other) is cvxpy_scalar_param):
            if op == SUMMATION:
                return cvxpy_tree(operator,args+[other])
            else:
                return cvxpy_tree(operator,args+[-other])

        # Tree
        elif type(other) is cvxpy_tree:
            if other.item.name == SUMMATION:
                if op == SUMMATION:
                    return cvxpy_tree(operator,
                                      args+other.children)
                else:
                    return cvxpy_tree(operator,
                                      args+[-x for x in other.children])
            else:
                if op == SUMMATION:
                    return cvxpy_tree(operator,args+[other])
                else:
                    return cvxpy_tree(operator,args+[-other])

        # Matrix
        elif type(other) is cvxpy_matrix:
            (m,n) = other.shape
            new_ar = cvxpy_array(m,n)
            for i in range(0,m,1):
                for j in range(0,n,1):
                    if op == SUMMATION:
                        new_ar[i,j] = self+other[i,j]
                    else:
                        new_ar[i,j] = self-other[i,j]
            return new_ar

        # Not implemented
        else:
            return NotImplemented

    # Method: __raddsub__
    def __raddsub__(self,other,op):
        """
        Handles right add and right subtract.

        :param other: Left hand side.
        :param op: Keyword (See cvxpy.defs)
        """

        # Create operator
        operator = cvxpy_obj(OPERATOR,np.nan,SUMMATION)

        # Create args to combine
        if (type(self) is cvxpy_tree and
            self.item.name == SUMMATION):
            if op == SUMMATION:
                args = self.children
            else:
                args = [-x for x in self.children]
        else:
            if op == SUMMATION:
                args = [self]
            else:
                args = [-self]

        # Number
        if np.isscalar(other):
            if other == 0.0:
                if op == SUMMATION:
                    return self
                else:
                    return -self
            constant = cvxpy_obj(CONSTANT,other,str(other))
            return cvxpy_tree(operator,[constant]+args)

        # Scalar variable or scalar param
        elif (type(other) is cvxpy_scalar_var or
              type(other) is cvxpy_scalar_param):
            return cvxpy_tree(operator,[other]+args)

        # Tree
        elif type(other) is cvxpy_tree:
            if other.item.name == SUMMATION:
                return cvxpy_tree(operator,other.children+args)
            else:
                return cvxpy_tree(operator,[other]+args)

        # Matrix
        elif type(other) is cvxpy_matrix:
            (m,n) = other.shape
            new_ar = cvxpy_array(m,n)
            for i in range(0,m,1):
                for j in range(0,n,1):
                    if op == SUMMATION:
                        new_ar[i,j] = other[i,j]+self
                    else:
                        new_ar[i,j] = other[i,j]-self
            return new_ar

        # Not implemented
        else:
            return NotImplemented

    # Method: __mul__
    def __mul__(self,other):

        return self.__mulhandle__(other)

    # Method: __rmul__
    def __rmul__(self,other):

        return self.__mulhandle__(other)

    # Method: __mulhandle__
    def __mulhandle__(self,other):
        """
        Handles left and right multiply.
        Important: If other is a constant, the multiplication
        tree formed has the constant object as the first child.

        :param other: Other operand.
        """

        # Create operator
        operator = cvxpy_obj(OPERATOR,np.nan,MULTIPLICATION)

        # Number
        if np.isscalar(other):

            # Multiplication by 0
            if other == 0.0:
                return 0.0

            # Multiplication by 1
            elif other == 1.0:
                return self

            # Distribution over addition
            elif (type(self) is cvxpy_tree and
                  self.item.name == SUMMATION):
                new_children = []
                for child in self.children:
                    new_children += [other*child]
                return cvxpy_tree(self.item,new_children)

            # Associativity
            elif (type(self) is cvxpy_tree and
                  self.item.name == MULTIPLICATION and
                  self.children[0].type == CONSTANT):
                new_object = other*self.children[0]
                if new_object.value == 1.0:
                    return self.children[1]
                else:
                    return cvxpy_tree(self.item,[new_object,
                                                 self.children[1]])

            # Constant times constant
            elif type(self) is cvxpy_obj:
                new_const = cvxpy_obj(CONSTANT,self.value*other,
                                      str(self.value*other))
                return new_const

            # Else
            else:
                constant = cvxpy_obj(CONSTANT,other,str(other))
                return cvxpy_tree(operator,[constant,self])

        # Matrix
        elif type(other) is cvxpy_matrix:
            (m,n) = other.shape
            new_ar = cvxpy_array(m,n)
            for i in range(0,m,1):
                for j in range(0,n,1):
                    new_ar[i,j] = other[i,j]*self
            return new_ar

        # Sparse matrix
        elif type(other) is cvxpy_spmatrix:
            (m,n) = other.shape
            new_ar = cvxpy_sparray(m,n)
            for i in range(0,m,1):
                rowi_indeces = other.rows[i]
                rowi_values = other.data[i]
                for k in range(0,len(rowi_indeces)):
                    new_ar[i,rowi_indeces[k]] = rowi_values[k]*self
            return new_ar

        # Scalar param
        elif type(other) is cvxpy_scalar_param:
            return cvxpy_tree(operator,[other,self])

        # Tree
        elif type(other) is cvxpy_tree:
            if not len(self.variables):
                return cvxpy_tree(operator,[self,other])
            elif not len(other.variables):
                return cvxpy_tree(operator,[other,self])
            else:
                return NotImplemented

        # Not implemented
        else:
            return NotImplemented

    # Method: __neg__
    def __neg__(self):

        return (-1.0)*self

    # Method: __str__
    def __str__(self):

        return self.name

#***********************************************************************#
# Class definition: cvxpy_scalar_var                                    #
#***********************************************************************#
class cvxpy_scalar_var(cvxpy_obj):

    # Variable counter
    i = 0

    # Method: __init__
    def __init__(self,name=None):
        """
        Class constructor.

        :param name: String.
        """

        # Assign a name
        if name is None:
            name = BASE_VAR+str(cvxpy_scalar_var.i)
        cvxpy_scalar_var.i += 1

        # Call parent constructor
        cvxpy_obj.__init__(self,VARIABLE,np.nan,name)

    # Method: __getattribute__
    def __getattribute__(self,name):

        # Variables
        if name == 'variables':
            return cvxpy_list([self])

        # Other
        else:
            return object.__getattribute__(self,name)

    # Method: is_nonnegative_constant
    def is_nonnegative_constant(self):
        """
        Determines if object is a nonnegative constant.
        """

        return False

    # Method: is_nonpositive_constant
    def is_nonpositive_constant(self):
        """
        Determines if object is a nonpositive constant.
        """

        return False

#***********************************************************************#
# Class definition: cvxpy_scalar_param                                  #
#***********************************************************************#
class cvxpy_scalar_param(cvxpy_obj):

    # Param counter
    i = 0

    # Method: __init__
    def __init__(self,attribute=None,name=None):
        """
        Class constructor.

        :param name: String.
        :param attribute: Keyword (See cvxpy.defs)
        """

        # Assign a name
        if name is None:
            name = BASE_PARAM+str(cvxpy_scalar_param.i)
        cvxpy_scalar_param.i += 1

        # Call parent constructor
        cvxpy_obj.__init__(self,PARAMETER,np.nan,name)

        # Store attribute
        self.attribute = attribute

    # Method: __getattribute__
    def __getattribute__(self,name):

        # Parameters
        if name == 'parameters':
            return cvxpy_list([self])

        # Other
        else:
            return object.__getattribute__(self,name)

    # Method: __setattr__
    def __setattr__(self,name,value):

        # Value
        if name == 'value':
            if not np.isscalar(value):
                raise TypeError('Invalid value type')
            if value < 0 and self.attribute == NONNEGATIVE:
                raise ValueError('Value must be nonnegative')
            if value > 0 and self.attribute == NONPOSITIVE:
                raise ValueError('Value must be nonpositive')

        # Attribute
        if name == 'attribute':
            if value not in [None,NONNEGATIVE,NONPOSITIVE]:
                raise ValueError("Invalid attribute")
            if value == NONNEGATIVE and self.value < 0:
                self.value = np.nan
            if value == NONPOSITIVE and self.value > 0:
                self.value = np.nan

        # Set attribute
        object.__setattr__(self,name,value)

    # Method: is_nonnegative_constant
    def is_nonnegative_constant(self):
        """
        Determines if object is a nonnegative constant.
        """

        return self.attribute == NONNEGATIVE

    # Method: is_nonpositive_constant
    def is_nonpositive_constant(self):
        """
        Determines if object is a nonpositive constant.
        """

        return self.attribute == NONPOSITIVE

#***********************************************************************#
# Class definition: cvxpy_tree                                          #
#***********************************************************************#
class cvxpy_tree(cvxpy_obj):

    # Method: __init__
    def __init__(self,item,children):
        """
        Class constructor.

        :param item: Operator or function.
        :param children: List of scalar objects.
        """

        self.item = item
        self.children = children
        cvxpy_obj.__init__(self,TREE,np.nan,'')

    # Method: __getattribute__
    def __getattribute__(self,name):

        # Variables
        if name == 'variables':
            l = list(map(lambda x: x.variables,self.children))
            return cvxpy_list(set(reduce(lambda x,y:x+y,l,[])))

        # Parameters
        elif name == 'parameters':
            l = list(map(lambda x: x.parameters,self.children))
            return cvxpy_list(set(reduce(lambda x,y:x+y,l,[])))

        # Value
        elif name == 'value':

            # Summation
            if (self.item.type == OPERATOR and
                self.item.name == SUMMATION):
                return np.sum(list(map(lambda x:x.value,self.children)))

            # Multiplication
            elif (self.item.type == OPERATOR and
                  self.item.name == MULTIPLICATION):
                return self.children[0].value*self.children[1].value

            # Function
            elif self.item.type == FUNCTION:
                return self.item(list(map(lambda x: x.value,self.children)))

            # Error
            else:
                raise TypeError('Invalid tree item')

        # Other
        else:
            return object.__getattribute__(self,name)

    # Method: __str__
    def __str__(self):

        # Multiplication
        if (self.item.type == OPERATOR and
            self.item.name == MULTIPLICATION):
            left_text = str(self.children[0])
            if (self.children[0].type == TREE and
                self.children[0].item.type == OPERATOR):
                left_text = '('+left_text+')'
            right_text = str(self.children[1])
            if (self.children[1].type == TREE and
                self.children[1].item.type == OPERATOR):
                    right_text = '('+right_text+')'
            return left_text+'*'+right_text

        # Summation
        elif (self.item.type == OPERATOR and
              self.item.name == SUMMATION):
            text = str(self.children[0])
            for i in range(1,len(self.children),1):
                x = self.children[i]
                text += ' + '+str(x)
            return text

        # Function
        elif self.item.type == FUNCTION:
            children_text = list(map(str,self.children))
            args = ''
            for i in range(0,len(children_text),1):
                args += children_text[i]
                if i != len(children_text)-1:
                    args += ', '
            return self.item.name+'('+args+')'

        # Error
        else:
            raise TypeError('Invalid tree item')

    # Method: is_convex
    def is_convex(self):
        """
        Determines if tree follows DCP rules of convexity.
        """

        # Summation
        if (self.item.type == OPERATOR and
            self.item.name == SUMMATION):
            return all(list(map(lambda x: x.is_convex(),self.children)))

        # Multiplication
        elif (self.item.type == OPERATOR and
              self.item.name == MULTIPLICATION):
            ob1 = self.children[0]
            ob2 = self.children[1]
            case1 = ob1.is_nonnegative_constant() and ob2.is_convex()
            case2 = ob1.is_nonpositive_constant() and ob2.is_concave()
            case3 = ob2.is_affine()
            return case1 or case2 or case3

        # Function
        elif self.item.type == FUNCTION:
            fcn = self.item
            args = self.children
            return (fcn.action == MINIMIZE) and fcn.is_dcp(args)

        # Error
        else:
            raise TypeError('Invalid tree item')

    # Method: is_concave
    def is_concave(self):
        """
        Determines if tree follows DCP rules of concavity.
        """

        # Summation
        if (self.item.type == OPERATOR and
            self.item.name == SUMMATION):
            return all(list(map(lambda x: x.is_concave(),self.children)))

        # Multiplication
        elif (self.item.type == OPERATOR and
              self.item.name == MULTIPLICATION):
            ob1 = self.children[0]
            ob2 = self.children[1]
            case1 = ob1.is_nonnegative_constant() and ob2.is_concave()
            case2 = ob1.is_nonpositive_constant() and ob2.is_convex()
            case3 = ob2.is_affine()
            return case1 or case2 or case3

        # Function
        elif self.item.type == FUNCTION:
            fcn = self.item
            args = self.children
            return (fcn.action == MAXIMIZE) and fcn.is_dcp(args)

        # Error
        else:
            raise TypeError('Invalid tree item')

    # Method: is_dcp
    def is_dcp(self):
        """
        Determines if tree is DCP-compliant.
        """

        return self.is_convex() or self.is_concave()

    # Method: is_affine
    def is_affine(self):
        """
        Determines if tree is an affine expression.
        """

        # Summation
        if (self.item.type == OPERATOR and
            self.item.name == SUMMATION):
            return all(list(map(lambda x: x.is_affine(),self.children)))

        # Multiplication
        elif (self.item.type == OPERATOR and
              self.item.name == MULTIPLICATION):
            ob2 = self.children[1]
            return ob2.is_affine()

        # Function
        elif self.item.type == FUNCTION:
            return len(cvxpy_list(self.children).variables) == 0

        # Error
        else:
            raise TypeError('Invalid tree item')

    # Method: is_nonnegative_constant
    def is_nonnegative_constant(self):
        """
        Determines if tree is guaranteed to be a
        nonnegative constant.
        """

        # Summation
        if (self.item.type == OPERATOR and
            self.item.name == SUMMATION):
            return all(list(map(lambda x: x.is_nonnegative_constant(),
                           self.children)))

        # Multiplication
        elif (self.item.type == OPERATOR and
              self.item.name == MULTIPLICATION):
            ob1 = self.children[0]
            ob2 = self.children[1]
            both_nonneg = all(list(map(lambda x:x.is_nonnegative_constant(),
                                  [ob1,ob2])))
            both_nonpos = all(list(map(lambda x:x.is_nonpositive_constant(),
                                  [ob1,ob2])))
            return both_nonneg or both_nonpos

        # Function
        elif self.item.type == FUNCTION:
            return False

        # Error
        else:
            raise TypeError('Invalid tree item')

    # Method: is_nonpositive_constant
    def is_nonpositive_constant(self):
        """
        Determines if tree is guaranteed to be a
        nonpositive constant.
        """

        # Summation
        if (self.item.type == OPERATOR and
            self.item.name == SUMMATION):
            return all(list(map(lambda x: x.is_nonpositive_constant(),
                           self.children)))

        # Multiplication
        elif (self.item.type == OPERATOR and
              self.item.name == MULTIPLICATION):
            ob1 = self.children[0]
            ob2 = self.children[1]
            p1 =(ob1.is_nonnegative_constant(),
                 ob1.is_nonpositive_constant())
            p2 =(ob2.is_nonnegative_constant(),
                 ob2.is_nonpositive_constant())
            return (p1[0] and p2[1]) or (p1[1] and p2[0])

        # Function
        elif self.item.type == FUNCTION:
            return False

        # Error
        else:
            raise TypeError('Invalid tree item')

# Load modules
from .arrays import cvxpy_array,cvxpy_sparray
from .arrays import cvxpy_matrix,cvxpy_spmatrix
from .constraints import cvxpy_list
