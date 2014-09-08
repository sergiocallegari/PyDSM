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
from ..defs import *
from ..scalars import *
from ..constraints import *
from ..interface import equals
from ..interface import less_equals
from ..interface import greater_equals
from ..interface import variable
from ..arrays import cvxpy_var
from ..arrays import cvxpy_array

# Function: expand
def expand(arg):
    
    # Constant
    if type(arg) is cvxpy_obj:
        return arg,[]
    
    # Scalar variable
    elif type(arg) is cvxpy_scalar_var:
        return arg,[]

    # Summation
    elif (type(arg) is cvxpy_tree and 
          arg.item.type == OPERATOR and
          arg.item.name == SUMMATION):

        # Get item and children
        item = arg.item
        children = arg.children

        # New variale
        v = variable()            
            
        # Expand children
        new_children = []
        new_constr = []
        for child in children:
            
            # Multiplication
            if (child.type == TREE and
                child.item.type == OPERATOR and
                child.item.name == MULTIPLICATION):
                child_var,child_constr = expand(child.children[1])
                new_children += [child.children[0].value*child_var]
                new_constr += child_constr
                    
            # Else
            else:
                child_var,child_constr = expand(child)
                new_children += [child_var]
                new_constr += child_constr
             
        # Return (Always right side is the new variable)
        new_tree = cvxpy_tree(item,new_children)
        return v,[equals(new_tree,v)]+new_constr
        
    # Multiplication
    elif (type(arg) is cvxpy_tree and 
          arg.item.type == OPERATOR and
          arg.item.name == MULTIPLICATION):

        # Get item and children
        item = arg.item
        children = arg.children

        # New variable
        v = variable()

        # Apply expand to second operand (first is a constant)
        child_var,child_constr = expand(children[1])

        # Return result (Always right side is the new variable)
        new_tree = cvxpy_tree(item,[children[0],child_var])
        new_eq = [equals(new_tree,v)]
        new_eq += child_constr
        return v,new_eq
    
    # Function
    elif (type(arg) is cvxpy_tree and 
          arg.item.type == FUNCTION):

        # Get item and children
        item = arg.item
        children = arg.children

        # New varibale 
        v = variable()

        # Analyze children
        new_children = []
        new_constr = []
        for child in children:
            child_var,child_constr = expand(child)
            new_children += [child_var]
            new_constr += child_constr
                
        # Return (Always right side is the new variable)
        new_tree = cvxpy_tree(item,new_children)
        if arg.is_convex():
            new_constr += [less_equals(new_tree,v)]
        else:
            new_constr += [greater_equals(new_tree,v)]
        return v,new_constr
    
    # Constraint
    elif type(arg) is cvxpy_constr:
        
        # Not set membership
        if arg.type != BELONGS:

            # Apply expand to left and right side
            obj1,constr_list1 = expand(arg.left)
            obj2,constr_list2 = expand(arg.right)
                                         
            # Return new constraints
            new_constr = cvxpy_constr(obj1,arg.type,obj2)
            new_list = [new_constr]
            new_list += constr_list1
            new_list += constr_list2
            return new_list

        # Set membership
        else:
            obj,constr_list = expand(arg.left)
            new_constr = cvxpy_constr(obj,arg.type,arg.right)
            return [new_constr]+constr_list

    # Array
    elif (type(arg) is cvxpy_array or
          type(arg) is cvxpy_var):
        (m,n) = arg.shape
        new_list = []
        new_ar = cvxpy_array(m,n)
        for i in range(0,m,1):
            for j in range(0,n,1):

                # Number: Upgrade
                if np.isscalar(arg[i,j]):
                    new_ar[i,j] = cvxpy_obj(CONSTANT,arg[i,j],str(arg[i,j]))
                    
                # Not a number
                else:
                    obj,constr_list = expand(arg[i,j])
                    new_ar[i,j] = obj
                    new_list += constr_list
        return new_ar,new_list
    
    # List of constraints
    elif (type(arg) is list or 
          type(arg) is cvxpy_list):
        return reduce(lambda x,y:x+y,list(map(expand,arg)),[])
        
    # Invalid
    else:
        raise TypeError('Invalid argument')
