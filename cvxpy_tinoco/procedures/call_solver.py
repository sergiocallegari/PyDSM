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
import cvxopt as opt
from .expand import expand
from ..defs import *
import scipy.sparse as sp
from cvxopt import solvers
from ..scalars import cvxpy_obj
from ..constraints import cvxpy_list

# Function: call_solver
def call_solver(p,quiet):
    """
    Calls solver. 

    :param p: Convex cvxpy_program. 
              Assumed to be expanded.
    :param quiet: Boolean.
    """

    # Set printing format for cvxopt sparse matrices
    opt.spmatrix_str = opt.printing.spmatrix_str_triplet 	
    
    # Expand objects defined via partial minimization
    constr_list = cvxpy_list(pm_expand(p.constraints))
    
    # Get variables
    variables = constr_list.variables
    variables.sort()

    # Count variables
    n = len(variables)

    # Create variable-index map
    var_to_index = {}
    for i in range(0,n,1):
        var_to_index[variables[i]] = i

    # Construct objective vector
    c = construct_c(p.objective,var_to_index,n,p.action)

    # Construct Ax == b
    A,b = construct_Ab(constr_list._get_eq(),var_to_index,n)

    # Construct  Gx <= h
    G,h,dim_l,dim_q,dim_s = construct_Gh(constr_list._get_ineq_in(),
                                         var_to_index,n)
    
    # Construct F
    F = construct_F(constr_list._get_ineq_in(),var_to_index,n)

    # Call cvxopt
    solvers.options['maxiters'] = p.options['maxiters']
    solvers.options['abstol'] = p.options['abstol']
    solvers.options['reltol'] = p.options['reltol']
    solvers.options['feastol'] = p.options['feastol']
    solvers.options['show_progress'] = not quiet
    dims = {'l':dim_l, 'q':dim_q, 's':dim_s}
    if F is None:
        r =  solvers.conelp(c,G,h,dims,A,b)
    else:
        r =  solvers.cpl(c,F,G,h,dims,A,b)

    # Store numerical values
    if r['status'] != PRIMAL_INFEASIBLE:
        for v in variables:
            v.value =  r['x'][var_to_index[v]]

    # Return result
    return r

# Function: construct_c
def construct_c(objective,mapping,n,action):
    """
    Creates the vector c for the objective function c.T*x

    :param objective: cvxpy_scalar_var or cvxpy_obj.
    :param mapping: Dictionary.
    :param n: Size for vector c.
    :param action: MINIMIZE OR MAXIMIZE
    """

    # Check objective
    if (objective.type != CONSTANT and 
        objective.type != VARIABLE):
        raise TypeError('Bad objective: Cannot construct c')

    # Construct c vector
    c = opt.matrix(0.0,(n,1))
    obj_vars = objective.variables
    if not len(obj_vars):
        return c
    elif len(obj_vars) == 1:
        obj_var = obj_vars[0]
        i = mapping[obj_var]
        if action == MINIMIZE:
            c[i] = 1.
        else:
            c[i] = -1.
    else:
        raise TypeError('Bad objective: Cannot construct c')

    # Return vector
    return c

# Function: construct_Ab
def construct_Ab(constr_list,mapping,n):
    """
    Constructs matrix A and vector b from a list 
    of equality constraints. 
    
    :param constr_list: List of equality constraints.
    :param mapping: Dictionary.
    :param n: Number of variables.
    """

    # Sort constraints
    def cmp(x,y):
        lx = len(x.variables)
        if ((type(x.left)  is cvxpy_obj and x.left.value != 0) or
            (type(x.right) is cvxpy_obj and x.right.value != 0)):
             lx += 1
        ly = len(y.variables)
        if ((type(y.left)  is cvxpy_obj and y.left.value != 0) or
            (type(y.right) is cvxpy_obj and y.right.value != 0)):
            ly += 1
        return 1 if lx >= ly else -1

    # Sort constraints keys
    def cmp_keys(x):
        lx = len(x.variables)
        if ((type(x.left)  is cvxpy_obj and x.left.value != 0) or
            (type(x.right) is cvxpy_obj and x.right.value != 0)):
            lx += 1
        return lx

    constr_list.sort(key=cmp_keys)
    
    # Construct matrices
    Ab = None
    for constr in constr_list:

        # Get elements
        left = constr.left
        right = constr.right

        # New row
        t = sp.lil_matrix((1,n+1))
        
        # Deal with right element
        if right.type == CONSTANT:
            t[0,n] += right.value*1.
        elif right.type == VARIABLE:
            t[0,mapping[right]] += -1.
        else:
            raise TypeError('Bad equality: Cannot construct A,b')

        # Left is a variable
        if left.type == VARIABLE:
            t[0,mapping[left]] += 1.
   
        # Left is a constant
        elif left.type == CONSTANT:
             t[0,n] += -left.value*1.

        # Left is an operation tree
        elif left.type == TREE and left.item.type == OPERATOR:
            
            # Addition
            if left.item.name == SUMMATION:
                for arg in left.children:
                    if (arg.type == TREE and
                        arg.item.type == OPERATOR and
                        arg.item.name == MULTIPLICATION):
                        ch1 = arg.children[1]
                        ch0 = arg.children[0]
                        t[0,mapping[ch1]] += ch0.value*1.
                    elif arg.type == VARIABLE:
                        t[0,mapping[arg]] += 1.
                    elif arg.type == CONSTANT:
                        t[0,n] += -arg.value*1.
                    else:
                        raise TypeError('Bad equality: Cannot construct A,b')

            # Multiplication
            elif left.item.name == MULTIPLICATION:
                op1 = left.children[0]
                op2 = left.children[1]
                t[0,mapping[op2]] += 1.*op1.value
            
            # Error
            else:
                raise TypeError('Bad equality: Cannot construct A,b')
    
        # Error
        else:
            raise TypeError('Bad equality: Cannot construct A,b')
        
        # Compute norm
        t = t.tocsr()
        t_norm = np.sqrt((t*t.T)[0,0])

        # Process
        if t_norm < EPSILON:
            continue
        if Ab is None:
            Ab = t
        else:
            t_c = t.nonzero()[1]
            Ab_c = Ab.nonzero()[1]
            c_overlap = True
            for x in t_c:
                if x not in Ab_c:
                    c_overlap = False
                    break
            valid_row = True
            if c_overlap:
                tn = t/t_norm
                for i in range(Ab.shape[0]-1,-1,-1):
                    s = Ab[i,:]
                    s_c = s.nonzero()[1]
                    if len(s_c) < len(t_c):
                        break
                    if set(s_c) == set(t_c):
                        e = s-(tn*s.T)*tn
                        if np.sqrt((e*e.T)[0,0]) < EPSILON:
                            valid_row = False
                            break
            if valid_row:
                Ab = sp.vstack((Ab,t)).tocsr()

    # Return matrices
    if Ab is None:
        Ab = opt.spmatrix(0.0,[],[],(0,n+1))
    else:
        nz = Ab.nonzero()
        Ab = opt.spmatrix(Ab.data.tolist(),nz[0].tolist(),
                          nz[1].tolist(),Ab.shape)
    return Ab[:,0:n],opt.matrix(Ab[:,n])

# Function: construct_Gh
def construct_Gh(constr_list,mapping,n):
    """
    Creates the matrix G and vector h from a list 
    of inequality and membership constraints.
    
    :param constr_list: List of inequality and membership constraints.
    :param mapping: Dictionary.
    :param n: Number of varaibles.
    """

    # Check constraints
    for c in constr_list:
        if c.left.type == ARRAY:
            if c.right.type != SET:
                raise TypeError('Bad constraint: Cannot construct G,h')
        else:
            if ((c.left.type != VARIABLE and c.left.type != CONSTANT) or
                (c.right.type != VARIABLE and c.right.type != CONSTANT)):
                raise TypeError('Bad constraint: Cannot construct G,h')

    # Initialize G,h and dimensions
    dim_l = 0
    dim_q = []
    dim_s = []
    G = opt.spmatrix(0.0,[],[],(0,n))
    h = opt.matrix(0.0,(0,1))

    # Nonnegative Orthant
    for c in constr_list:
        
        # Left is constant or variable
        if c.left.type == CONSTANT or c.left.type == VARIABLE:

            # Get constraint elements
            t = c.type
            ob1 = c.left
            ob2 = c.right

            # New row
            rowG = opt.spmatrix(0.0,[],[],(1,n))
            rowh = opt.matrix(0.0,(1,1))
            if ob1.type == VARIABLE:
                rowG[0,mapping[ob1]]+=1. if t==LESS_EQUALS else -1.
            else:
                rowh[0,0]+=-ob1.value*1. if t==LESS_EQUALS else ob1.value*1.
            if ob2.type == VARIABLE:
                rowG[0,mapping[ob2]]+=1. if t==GREATER_EQUALS else -1.
            else:
                rowh[0,0]+=-ob2.value*1. if t==GREATER_EQUALS else ob2.value*1.

            # Append to G,h
            G = opt.sparse([G,rowG])
            h = opt.matrix([h,rowh])

            # Increment size of cone
            dim_l += 1
        
    # Second order cone
    for c in constr_list:

        # Set membership that expands to SOC
        if (c.right.type == SET and
            c.right.expansion_type == SOC):
        
            # Get G,h, section
            el = c.left
            set_atom = c.right
            newG,newh,r = set_atom._construct(el,mapping,n)

            # Attach to G,h
            G = opt.sparse([G,newG])
            h = opt.matrix([h,newh])

            # Attach size of cone
            dim_q = dim_q + [r]

    # Semidefinite cone
    for c in constr_list:

        # Set membership that expands to LMI
        if (c.right.type == SET and
            c.right.expansion_type == SDC):
            
            # Get G,h, section
            el = c.left
            set_atom = c.right
            newG,newh,t = set_atom._construct(el,mapping,n)

            # Attach to G,h
            G = opt.sparse([G,newG])
            h = opt.matrix([h,newh])

            # Attach size of cone
            dim_s = dim_s + [t]
      
    # Return 
    return G,h,dim_l,dim_q,dim_s

# Function
def construct_F(constr_list,mapping,n):
    """
    Constructs function F needed by cvxopt 
    to solve nonlinear programs.
    
    :param constr_list: List of inequality and set 
                        membership constraints.
    :param mapping: Dictionary.
    :param n: Numbe of variables.
    """
    
    # Check constraints
    for c in constr_list:
        if c.left.type == ARRAY:
            if c.right.type != SET:
                raise TypeError('Bad constraint: Cannot construct F')
        else:
            if ((c.left.type != VARIABLE and c.left.type != CONSTANT) or
                (c.right.type != VARIABLE and c.right.type != CONSTANT)):
                raise TypeError('Bad constraint: Cannot construct F')

    # Lists
    fs = []
    grads = []
    hess = []
    inds = []
    
    for c in constr_list:
        
        # Set that has a (f,gradf,hessf) implementation
        if (c.right.type == SET and
            c.right.expansion_type == DIF):

            # Get objects
            el = c.left
            set_atom = c.right

            # Get functions
            n_f,n_grad,n_hess,n_ind = set_atom._construct(el,mapping,n)
            
            # Append
            fs = fs + [n_f]
            grads = grads + [n_grad]
            hess = hess + [n_hess]
            inds = inds + [n_ind]
        
    # Check if no functions
    if not len(fs):
        return None

    # Construct F
    def F(x=None,z=None):
        
        # Case 1
        if x is None and z is None:
            x0 = opt.matrix(1., (n,1))
            return len(fs),x0

        # Case 2
        elif x is not None and z is None:
            if all(list(map(lambda y: y(x),inds))):
                f = opt.matrix(0.0,(len(fs),1))
                for i in range(0,len(fs),1):
                    f[i] = fs[i](x)
                Df = opt.spmatrix(0.0,[],[],(0,n))
                for i in range(0,len(grads),1):
                    Df = opt.sparse([Df,grads[i](x).T])
                return f,Df
            else:
                return None,None

        # Case 3
        else:
            f = opt.matrix(0.0,(len(fs),1))
            for i in range(0,len(fs),1):
                f[i] = fs[i](x)
            Df = opt.spmatrix(0.0,[],[],(0,n))
            for i in range(0,len(grads),1):
                Df = opt.sparse([Df,grads[i](x).T])
            H = opt.spmatrix(0.0,[],[],(n,n))
            for i in range(0,len(hess),1):
                H = H + z[i]*hess[i](x)
            return f,Df,H

    # Return F
    return F

# Function: Partial minimization expansion
def pm_expand(constr_list):
    """
    Expands functions that are implemented
    using partial minimization descriptions.

    :param constr_list: List of constraints.
    """

    new_list = []
    for c in constr_list:
        if (c.left.type == TREE and
            c.left.item.type == FUNCTION):
            new_constr = expand(c.left.item._pm_expand(c))
            new_list += pm_expand(new_constr)
        elif (c.right.type == SET and
              c.right.expansion_type == PM):
            new_constr = expand(c.right._pm_expand(c))
            new_list += pm_expand(new_constr)
        else:
            new_list += [c]

    # Return new list
    return new_list
