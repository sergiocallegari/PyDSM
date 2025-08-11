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
from .call_solver import call_solver

# Function
def solve_prog(p,quiet):
    """
    Solves optimization program.

    :param p: cvxpy_program
    """

    # Expand
    p_expanded = p._get_expanded_program()

    # Compute signs
    if p.action == MINIMIZE:
        sign = 1.
    else:
        sign = -1.

    # Solve convex relaxation
    if not quiet:
        print('\nCalling CVXOPT ...')
    sol = call_solver(p_expanded,quiet)

    valid = True

    # Primal infeasible
    if sol['status'] == PRIMAL_INFEASIBLE:
        obj = np.inf*sign

    # Uknown
    elif sol['status'] == UNKNOWN:
        obj = sol[PRIMAL_OBJECTIVE]*sign
        valid = False

    # Dual infeasible
    elif sol['status'] == DUAL_INFEASIBLE:
        obj = np.inf*(-sign)

    # Optimal
    else:
        obj =sol[PRIMAL_OBJECTIVE]*sign

    # Return optimal value
    return obj,valid
