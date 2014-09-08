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
from cvxpy.utils import sum
from cvxpy.arrays import *
from cvxpy.scalars import *
from cvxpy.constraints import *
from cvxpy.interface import equals
from cvxpy.interface import belongs
from cvxpy.interface import less_equals
from cvxpy.interface import greater_equals

# Function: re_eval
def re_eval(arg,replace_map):
    """
    Re-evaluates argument.

    :param replace_map: Dictionary.
    """

    # Number
    if np.isscalar(arg):
        return arg
    
    # Constant object
    elif type(arg) is cvxpy_obj:
        return arg.value
    
    # Scalar variable
    elif type(arg) is cvxpy_scalar_var:
        if arg in replace_map.keys():
            return re_eval(replace_map[arg],replace_map)
        else:
            return arg
    
    # Scalar param
    elif type(arg) is cvxpy_scalar_param:
        if arg in replace_map.keys():
            return re_eval(replace_map[arg],replace_map)
        else:
            return arg

    # Summation
    elif (type(arg) is cvxpy_tree and 
          arg.item.type == OPERATOR and
          arg.item.name == SUMMATION):
        new_children = list(map(lambda x:re_eval(x,replace_map),arg.children))
        return sum(new_children)

    # Multiplication
    elif (type(arg) is cvxpy_tree and 
          arg.item.type == OPERATOR and
          arg.item.name == MULTIPLICATION):
        child1 = re_eval(arg.children[0],replace_map)
        child2 = re_eval(arg.children[1],replace_map)
        return child1*child2

    # Function
    elif (type(arg) is cvxpy_tree and 
          arg.item.type == FUNCTION):
        new_children = list(map(lambda x:re_eval(x,replace_map),arg.children))
        return arg.item(new_children)

    # Constraint
    elif type(arg) is cvxpy_constr:
        
        # Not set membership
        if arg.type != BELONGS:
            left = re_eval(arg.left ,replace_map)
            right= re_eval(arg.right,replace_map)
            if arg.type == EQUALS:
                return equals(left,right)
            elif arg.type == LESS_EQUALS:
                return less_equals(left,right)
            elif arg.type == GREATER_EQUALS:
                return greater_equals(left,right)
            else:
                raise TypeError('Invalid constraint')

        # Set membership
        else:
            left = re_eval(arg.left,replace_map)
            return belongs(left,arg.right)

    # Array
    elif (type(arg) is cvxpy_array or
          type(arg) is cvxpy_var or
          type(arg) is cvxpy_param):
        (m,n) = arg.shape
        new_ar = cvxpy_array(m,n)
        for i in range(0,m,1):
            for j in range(0,n,1):
                new_ar[i,j] = re_eval(arg[i,j],replace_map)
        return new_ar

    # List
    elif (type(arg) is list or
          type(arg) is cvxpy_list):

        new_list = []
        for c in arg:
            new_list += [re_eval(c,replace_map)]
        return cvxpy_list(new_list)

    # Invalid
    else:
        raise TypeError('Invalid argument')
