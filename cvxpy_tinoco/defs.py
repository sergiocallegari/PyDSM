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

# Names
BASE_VAR = 'v'
BASE_PARAM = 'a'

# Epsilon
EPSILON = 1e-12

# Variable structure
SYMMETRIC = 'symmetric'
LOWER_TRIANGULAR = 'lower_triangular'
UPPER_TRIANGULAR = 'upper_triangular'

# Parameter attribute
NONNEGATIVE = 'nonnegative'
NONPOSITIVE = 'nonpositive'

# Groups
VARIABLE_OBJS = ['cvxpy_scalar_var','cvxpy_var']
PARAMETER_OBJS = ['cvxpy_scalar_param','cvxpy_param']
SCALAR_OBJS = ['cvxpy_scalar_var','cvxpy_scalar_param','cvxpy_tree']
ARRAY_OBJS = ['cvxpy_var','cvxpy_param','cvxpy_array']
MATRIX_OBJS = ['cvxpy_matrix']
FORMAL_OBJS = VARIABLE_OBJS + PARAMETER_OBJS

# Object types
SET = 'set'
TREE = 'tree'
ARRAY = 'array'
OPERATOR = 'operator'
CONSTANT = 'constant'
FUNCTION = 'function'
VARIABLE = 'variable'
PARAMETER = 'parameter'
SPARRAY = 'sparray'

# Operations
MULTIPLICATION = '*'
SUMMATION = '+'
SUBTRACTION = '-'
LEFT_ADD = 'la'
LEFT_SUBTRACT = 'ls'
RIGHT_ADD = 'ra'
RIGHT_SUBTRACT = 'rs'
LEFT_MULTIPLY = 'lm'
RIGHT_MULTIPLY = 'rm'

# Expansion types
SDC = 'sdc'
SOC = 'soc'
DIF = 'dif'
PM = 'pm'

# Actions
MINIMIZE = 'minimize'
MAXIMIZE = 'maximize'
 
# Constraint types
EQUALS = '=='
LESS_EQUALS = '<='
GREATER_EQUALS = '>='
BELONGS = 'in'

# Configuration
CONFIGURATION ={'maxiters':100,
                'abstol':1e-7,
                'reltol':1e-6,
                'feastol':1e-6}

# Solver fields
OPTIMAL = 'optimal'
UNKNOWN = 'unknown'
DUAL_INFEASIBLE = 'dual infeasible'
PRIMAL_INFEASIBLE = 'primal infeasible'
PRIMAL_OBJECTIVE = 'primal objective'

