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
from .arrays import *
from .scalars import *
from .constraints import *
from .procedures.expand import expand
from .procedures.re_eval import re_eval
from .procedures.solve_prog import solve_prog

#***********************************************************************#
# Class definition: cvxpy_program                                       #
#***********************************************************************#
class cvxpy_program(object):

    # Method: __init__
    def __init__(self,action,objective,constraints,
                formals,options,name):
        """
        Class constructor.

        :param action: Keyword (see cvxpy.def).
        :param objective: Scalar expression.
        :param constraints: List of constraints.
        :param formals: List of variables or parameters.
        :para options: Dictionary.
        :param name: String.
        """
        
        # Attributes
        self.name = name
        self.type = FUNCTION
        self.action = action
        self.constraints = cvxpy_list(constraints)
        
        # Objective
        if np.isscalar(objective):
            self.objective = cvxpy_obj(CONSTANT,objective,str(objective))
        elif (type(objective) is cvxpy_obj or
              type(objective).__name__ in SCALAR_OBJS):
            self.objective = objective
        else:
            raise TypeError('Invalid objective')

        # Formals
        for x in formals:
            if type(x).__name__ not in FORMAL_OBJS:
                raise TypeError('Invalid formal')
        self.formals = cvxpy_list(formals)

        # Options
        if options != None:
            self.options = options
        else:
            self.options = CONFIGURATION.copy()

    # Method: _get_expanded_objects
    def _get_expanded_objects(self,args):
        """
        Expands list of objects to list of scalar objects.

        :param args: List of objects.
        """

        args_expanded = []
        for x in args:
            if (np.isscalar(x) or 
                type(x).__name__ in SCALAR_OBJS):
                args_expanded += [x]
            elif (type(x) is cvxpy_matrix or 
                  type(x).__name__ in ARRAY_OBJS):
                for i in range(0,x.shape[0],1):
                    for j in range(0,x.shape[1],1):
                        args_expanded += [x[i,j]]
            else:
                raise TypeError('Invalid argument')
        return args_expanded

    # Method: _upgrade_scalars
    def _upgrade_scalars(self,args):
        """
        Upgrades scalars to cvxpy_obj.

        :param args: List of objects.
        """

        args_upgraded = []
        for x in args:
            if np.isscalar(x):
                args_upgraded += [cvxpy_obj(CONSTANT,x,str(x))]
            else:
                args_upgraded += [x]
        return args_upgraded

    # Method: _get_expanded_program
    def _get_expanded_program(self):
        """
        Constructs a new program by applying the
        expansion algorithm to the objective
        and constraints.
        """

        new_obj,obj_constr = expand(self.objective)
        more_constr = expand(self.constraints)
        return cvxpy_program(self.action,new_obj,more_constr+obj_constr,
                             [],self.options,self.name)
    
    # Method: __getattribute__
    def __getattribute__(self,name):

        # Variables
        if name == 'variables':
            return cvxpy_list(set(self.constraints.variables + 
                                  self.objective.variables))

        # Parameters
        elif name == 'parameters':
            return cvxpy_list(set(self.constraints.parameters + 
                                  self.objective.parameters))

        # Other
        else:
            return object.__getattribute__(self,name)

    # Method: solve
    def solve(self,quiet=False,return_status=False):
        """
        Solves optimization program using current
        parameter values.
        """

        # Check DCP
        if not self.is_dcp():
            raise ValueError('Program is not DCP')

        # Create param-value map
        replace_map = {}
        for param in self.parameters:
            if np.isnan(param.value):
                raise ValueError('Invalid parameter value: NaN')
            else:
                replace_map[param] = param.value
        
        # Create new program
        new_p = cvxpy_program(self.action,
                              re_eval(self.objective,replace_map),
                              re_eval(self.constraints,replace_map),
                              [],self.options,'')

        # Solve new program
        obj,valid = solve_prog(new_p,quiet)
        if return_status:
            return obj,valid
        else:
            return obj
    
    # Method: __str__
    def __str__(self):

        return self.name

    # Method: show
    def show(self):
        """
        Shows description of the optimization program.
        """

        if self.action == MINIMIZE:
            output = '\nminimize '
        else:
            output = '\nmaximize '
        output += str(self.objective)+'\n'
        output += 'subject to\n'
        output += str(self.constraints)
        print(output)

    # Method: __call__
    def __call__(self,*args):
        """
        Invokes program.

        :param args: List of arguments.        
        """
        
        # Get arguments reference
        args_ref = self._get_expanded_objects(self.formals)
        
        # Internal call
        if len(args) == 1 and type(args[0]) is list:
            
            # Get list of arguments
            args_list = args[0]

            # Check number of arguments
            if len(args_list) != len(args_ref):
                raise TypeError('Invalid number of arguments')

        # User call
        else:
     
            # Check number of arguments
            if len(args) != len(self.formals):
                raise TypeError('Invalid number of arguments')

            # Check syntax
            for i in range(0,len(args),1):
                if (np.isscalar(args[i]) or 
                    type(args[i]).__name__ in SCALAR_OBJS):
                    if self.formals[i].shape != (1,1):
                        raise ValueError('Invalid argument shape')
                elif (type(args[i]) is cvxpy_matrix or 
                      type(args[i]).__name__ in ARRAY_OBJS):
                    if self.formals[i].shape != args[i].shape:
                        raise ValueError('Invalid argument shape')
                else:
                    raise TypeError('Invalid argument')

            # Expand arguments
            args_list = self._get_expanded_objects(args)

        # Verify values
        for i in range(0,len(args_ref)):
            if np.isscalar(args_list[i]) and np.isnan(args_list[i]):
                raise ValueError('Invalid argument value: NaN')

        # Upgrade scalars to objects
        args_upgraded = self._upgrade_scalars(args_list)

        # Verify replacements
        for i in range(0,len(args_ref)):
            if (len(args_upgraded[i].variables) != 0 and
                type(args_ref[i]) is not cvxpy_scalar_var):
                raise TypeError('Invalid replacement')

        # All numeric
        if all(list(map(lambda x: type(x) is cvxpy_obj,args_upgraded))):

            # Create map and new program
            replace_map = {}
            for i in range(0,len(args_ref),1):
                replace_map[args_ref[i]] = args_upgraded[i]
            new_p = cvxpy_program(self.action,
                                  re_eval(self.objective,replace_map),
                                  re_eval(self.constraints,replace_map),
                                  [],self.options,'')

            # Solve
            obj,valid = new_p.solve(quiet=True,return_status=True)
            if not valid:
                raise ValueError('Unable to compute value')
            else:
                return obj
        
        # Not all numeric
        else:
            
            # Return expression tree
            return cvxpy_tree(self,args_upgraded)

    # Method: is_dcp
    def is_dcp(self,args=None):
        """
        Determines if program is DCP-compliant.

        :param args: List of scalar arguments.
        """

        # No arguments: Check body of program
        if args == None:
            if (self.action == MINIMIZE and
                not self.objective.is_convex()):
                return False
            elif (self.action == MAXIMIZE and
                  not self.objective.is_concave()):
                return False
            else:
                return self.constraints.is_dcp()

        # Arguments given: Replace arguments and then check
        else:

            # Create reference list of arguments
            args_ref = self._get_expanded_objects(self.formals)

            # Check number of arguments
            if len(args) != len(args_ref):
                raise TypeError('Invalid number of arguments')

            # Create map and new program
            replace_map = {}
            for k in range(0,len(args),1):
                replace_map[args_ref[k]] = args[k]
            new_p = cvxpy_program(self.action,
                                  re_eval(self.objective,replace_map),
                                  re_eval(self.constraints,replace_map),
                                  [],self.options,'')

            # Check dcp on resulting program
            return new_p.is_dcp()
    
    # Method: _pm_expand
    def _pm_expand(self,constr):
        """
        Returns constraints that embed the program into
        the parent program.

        :param constr: cvxpy_constr of the form 
                       self(args), comparison-type, variable.
        """

        # Get arguments
        args = constr.left.children

        # Create argument reference list
        args_ref = self._get_expanded_objects(self.formals)
        
        # Check number of arguments
        if len(args) != len(args_ref):
            raise TypeError("Invalid number of arguments")

        # Create map and new program
        replace_map = {}
        for k in range(0,len(args),1):
            replace_map[args_ref[k]] = args[k]
        new_p = cvxpy_program(self.action,
                              re_eval(self.objective,replace_map),
                              re_eval(self.constraints,replace_map),
                              [],self.options,'')

        # Embed
        right = constr.right
        new_constr = []
        if self.action == MINIMIZE:
            new_constr += [compare(new_p.objective,LESS_EQUALS,right)]
        else:
            new_constr += [compare(new_p.objective,GREATER_EQUALS,right)]
        new_constr += new_p.constraints
        
        # Return constraints
        return new_constr
