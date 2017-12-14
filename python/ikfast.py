#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement # for python 2.5

__author__    = 'Rosen Diankov'
__copyright__ = 'Copyright (C) 2009-2012 Rosen Diankov <rosen.diankov@gmail.com>'
__license__   = 'Lesser GPL, Version 3'
__version__   = '0x1000004a' # hex of the version, has to be prefixed with 0x. also in ikfast.h

# ========== TGN's modification of IPython ==========
import os
def clc():
    os.system('clear')
    os.system('clear')

import traceback
def ikfast_print_stack():
    tb = traceback.extract_stack()
    pattern = '%-30s %5s %24s' 
    print( '\n'+pattern % ('        FUNCTION','LINE', 'FILE      '))
    keyword_of_interest = [ 'ikfast_IKFastSolver.py', 'ikfast_AST.py', 'ikfast.py', 'inversekinematics.py']
    print('--------------------------------------------------------------')
    for function_call in tb:
        for keyword in keyword_of_interest:
            if (keyword in function_call[0]) and (function_call[2] not in 'ikfast_print_stack'):
                print(pattern % (function_call[2], function_call[1], keyword))
                break

ipython_str = 'ikfast_print_stack(); ' + \
              'from IPython.terminal import embed; ' + \
              'ipshell = embed.InteractiveShellEmbed(banner1="", config=embed.load_default_config())(local_ns=locals())'

def print_matrix(matrices, ind=None):
    if ind is None:
        ind = range(len(matrices))
    for i in ind:
        print i, ':',
        print matrices[i]
        print '\n',

LOGGING_FORMAT = ' %(levelname)-6s [%(filename)s : %(funcName)s : %(lineno)d]\n' + \
                 '\t%(message)s\n'
# ========== End of TGN's modification ==============

# numpy for computing eigenvalues
import numpy
# check sympy version
from sympy import __version__ as sympy_version
if sympy_version < '0.7.0':
    raise ImportError('ikfast needs sympy 0.7.x or greater')
sympy_smaller_073 = sympy_version < '0.7.3'
if sympy_version > '0.7.1':
    _zeros, _ones = zeros, ones
    zeros = lambda args: _zeros(*args)
    ones = lambda args: _ones(*args)
# import rest of sympy
from sympy import *

from itertools import izip, chain, product
try:
    from itertools import combinations, permutations
except ImportError:
    def permutations(iterable, r=None):
        # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
        # permutations(range(3)) --> 012 021 102 120 201 210
        pool = tuple(iterable)
        n = len(pool)
        r = n if r is None else r
        if r > n:
            return
        indices = range(n)
        cycles = range(n, n-r, -1)
        yield tuple(pool[i] for i in indices[:r])
        while n:
            for i in reversed(range(r)):
                cycles[i] -= 1
                if cycles[i] == 0:
                    indices[i:] = indices[i+1:] + indices[i:i+1]
                    cycles[i] = n - i
                else:
                    j = cycles[i]
                    indices[i], indices[-j] = indices[-j], indices[i]
                    yield tuple(pool[i] for i in indices[:r])
                    break
            else:
                return

# changes to sympy:
# core/power.py Pow
def Pow_eval_subs(self, old, new):
    if self == old:
        return new
    
    if old.func is self.func and self.base == old.base:
        coeff1, terms1 = self.exp.as_coeff_mul()
        coeff2, terms2 = old.exp.as_coeff_mul()
        if terms1==terms2:
#             pow = coeff1/coeff2
#             if pow.is_Integer or self.base.is_commutative:
#                 return Pow(new, pow) # (x**(2*y)).subs(x**(3*y),z) -> z**(2/3)
            # only divide if coeff2 is a divisor of coeff1
            if coeff1.is_integer and coeff2.is_integer and (coeff1/coeff2).is_integer:
                return new ** (coeff1/coeff2) # (x**(2*y)).subs(x**(3*y),z) -> z**(2/3*y)
            
    if old.func is C.exp:
        coeff1, terms1 = old.args[0].as_coeff_mul()
        coeff2, terms2 = (self.exp*C.log(self.base)).as_coeff_mul()
        if terms1==terms2:
            # only divide if coeff2 is a divisor of coeff1
            if coeff1.is_integer and coeff2.is_integer and (coeff1/coeff2).is_integer:
                return new ** (coeff1/coeff2) # (x**(2*y)).subs(exp(3*y*log(x)),z) -> z**(2/3*y)
            
    return Pow(self.base._eval_subs(old, new), self.exp._eval_subs(old, new))

# simplify/simplify.py
def trigsimp_custom(self, **args):
    """
    Default trigsimp (in sympy >= 0.7.3) reduces sum of sin/cos products, for example

        trigsimp(-sin(x)⋅cos(y) + sin(y)⋅cos(x))
        >> -sin(x - y)

    We have to undo this step, which is what happens here.
    """
    from sympy.simplify import trigsimp as sympy_trigsimp 
    from sympy.simplify.fu import TR10
    return TR10(sympy_trigsimp(self, **args))

if not sympy_smaller_073:
    trigsimp = trigsimp_custom
else:
    power.Pow._eval_subs = Pow_eval_subs

from sympy.core import function # for sympy 0.7.1+

class fmod(function.Function):
    """defines floating-point mod"""
    nargs = 2
    is_real = True
    is_Function = True

class atan2check(atan2):
    """defines floating-point mod"""
    nargs = 2
    is_real = True
    is_Function = True

class GinacUtils:
    @staticmethod
    def ConvertToGinac(eq,localsymbolmap):
        if eq.is_Add:
            geq = None
            for arg in eq.args:
                geq2 = GinacUtils.ConvertToGinac(arg,localsymbolmap)
                if geq is None:
                    geq = geq2
                else:
                    geq += geq2
            return geq if geq is not None else swiginac.numeric(0)

        elif eq.is_Mul:
            geq = None
            for arg in eq.args:
                geq2 = GinacUtils.ConvertToGinac(arg,localsymbolmap)
                if geq is None:
                    geq = geq2
                else:
                    geq *= geq2
            return geq if geq is not None else swiginac.numeric(1)

        elif eq.is_Pow:
            gbase = GinacUtils.ConvertToGinac(eq.base,localsymbolmap)
            if eq.exp == S.One:
                return gbase

            elif eq.exp == -S.One:
                return 1/gbase

            else:
                return pow(gbase,GinacUtils.ConvertToGinac(eq.exp,localsymbolmap))

        elif eq.is_number:
            return swiginac.numeric(str(eq))

        elif eq.is_Symbol:
            if str(eq) in localsymbolmap:
                return localsymbolmap[str(eq)]

            else:
                gsym = swiginac.symbol(str(eq))
                localsymbolmap[str(eq)] = gsym
                return gsym

        raise ValueError('unknown equation %s'%str(eq))

    @staticmethod
    def ConvertFromGinac(geq):
        if isinstance(geq, swiginac.add):
            return Add(*[GinacUtils.ConvertFromGinac(geq.op(i)) for i in range(geq.nops())])
        
        elif isinstance(geq, swiginac.mul):
            return Mul(*[GinacUtils.ConvertFromGinac(geq.op(i)) for i in range(geq.nops())])
        
        elif isinstance(geq, swiginac.power):
            ebase = GinacUtils.ConvertFromGinac(geq.op(0))
            if geq.op(1) == 1:
                return ebase

            elif geq.op(1) == -1:
                return S.One/ebase
            
            else:
                return Pow(ebase,GinacUtils.ConvertFromGinac(geq.op(1)),evaluate=False)
            
        elif isinstance(geq, swiginac.numeric):
            if geq.is_integer():
                return Integer(str(geq))

            elif geq.is_rational():
                return Rational(str(geq))

            else:
                return geq.eval()

        elif isinstance(geq, swiginac.symbol):
            return Symbol(str(geq))

        else:
            raise ValueError('unknown equation %s'%str(eq))

    @staticmethod
    def ConvertMatrixToGinac(M,name='M', localsymbolmap={}):
        gM = swiginac.symbolic_matrix(M.shape[0],M.shape[1],'M')
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                gM[i,j] = GinacUtils.ConvertToGinac(M[i,j],localsymbolmap)
        return gM

    @staticmethod
    def GetPolyTermsFromGinac(geq, gothersymbols, othersymbols):
        """return a dict of monom:coeff items
        """
        terms = {}
        for i, gothersymbol in enumerate(gothersymbols):
            for degree in range(geq.ldegree(gothersymbol),geq.degree(gothersymbol)+1):
                monomprefix = (0,)*i + (degree,)
                gcoeff = geq.coeff(gothersymbol,degree)
                if i+1 < len(gothersymbols):
                    newterms = GinacUtils.GetPolyTermsFromGinac(gcoeff,gothersymbols[i+1:],othersymbols[i+1:])
                    for newmonom, newcoeff in newterms.iteritems():
                        assert(len(newmonom)==len(gothersymbols)-i-1)
                        terms[monomprefix+newmonom] = newcoeff
                else:
                    # ConvertFromGinac is very slow
                    terms[monomprefix] = gcoeff#GinacUtils.ConvertFromGinac(gcoeff)
        return terms

    @staticmethod
    def SolveUpperTriangular(gA, gB, name='X'):
        """solves for gA * X = gB.
        All parameters have to be ginac objects
        """
        gX = swiginac.symbolic_matrix(gB.rows(),gB.cols(),name)
        for i in reversed(xrange(gA.rows())):
            if gA[i, i] == 0:
                raise ValueError("Matrix must be non-singular.")
            
            gX[i, 0] = (gB[i, 0] - sum(gA[i, k] * gX[k, 0] for k in xrange(i+1, gA.rows()))) / gA[i, i]
            
        return gX

from ikfast_IKFastSolver import IKFastSolver

if __name__ == '__main__':
    import openravepy
    parser = OptionParser(description="""IKFast: The Robot Kinematics Compiler                                             
Software License Agreement (Lesser GPL v3). 
Copyright (C) 2009-2011 Rosen Diankov. 
IKFast is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

IKFast is part of OpenRAVE. This program can be used with robots or kinbodies defined and is independent of the OpenRAVE databases.
                                                  
Example usage for 7 DOF Barrett WAM where the 3rd joint is a free parameter:                               

python ikfast.py --robot=robots/barrettwam.robot.xml --baselink=0 --eelink=7 --savefile=ik.cpp --freeindex=2

""",version=__version__)
    parser.add_option('--robot', action='store', type='string', dest='robot',default=None,
                      help='robot file (COLLADA or OpenRAVE XML)')
    parser.add_option('--savefile', action='store', type='string', dest='savefile',default='ik.cpp',
                      help='filename where to store the generated c++ code')
    parser.add_option('--baselink', action='store', type='int', dest='baselink',
                      help='base link index to start extraction of ik chain')
    parser.add_option('--eelink', action='store', type='int', dest='eelink',
                      help='end effector link index to end extraction of ik chain')
    parser.add_option('--freeindex', action='append', type='int', dest='freeindices',default=[],
                      help='Optional joint index specifying a free parameter of the manipulator. If not specified, assumes all joints not solving for are free parameters. Can be specified multiple times for multiple free parameters.')
    parser.add_option('--iktype', action='store', dest='iktype',default='transform6d',
                      help='The iktype to generate the ik for. Possible values are: %s'%(', '.join(name for name,fn in IKFastSolver.GetSolvers().iteritems())))
    parser.add_option('--maxcasedepth', action='store', type='int', dest='maxcasedepth',default=3,
                      help='The max depth to go into degenerate cases. If ikfast file is too big, try reducing this, (default=%default).')
    parser.add_option('--lang', action='store',type='string',dest='lang',default='cpp',
                      help='The language to generate the code in (default=%default), available=('+','.join(name for name,value in CodeGenerators.iteritems())+')')
    parser.add_option('--debug','-d', action='store', type='int',dest='debug',default=logging.INFO,
                      help='Debug level for python nose (smaller values allow more text).')
    
    (options, args) = parser.parse_args()
    if options.robot is None or options.baselink is None or options.eelink is None:
        print('Error: Not all arguments specified')
        sys.exit(1)

    format = logging.Formatter('%(levelname)s: %(message)s')
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(format)
    log.addHandler(handler)
    log.setLevel(options.debug)

    solvefn=IKFastSolver.GetSolvers()[options.iktype.lower()]
    if options.robot is not None:
        try:
            env=openravepy.Environment()
            kinbody=env.ReadRobotXMLFile(options.robot)
            env.Add(kinbody)
            solver = IKFastSolver(kinbody,kinbody)
            solver.maxcasedepth = options.maxcasedepth
            chaintree = solver.generateIkSolver(options.baselink,options.eelink,options.freeindices,solvefn=solvefn)
            code=solver.writeIkSolver(chaintree,lang=options.lang)
        finally:
            openravepy.RaveDestroy()

    if len(code) > 0:
        with open(options.savefile,'w') as f:
            f.write(code)
