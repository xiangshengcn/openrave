from sympy import *

class AST:
    """Abstarct Syntax Tree class definitions specific for evaluating complex math equations.
    """
    
    class SolverBase(object):
        def GetChildrenOfType(self, classinstance):
            return []
        def GetZeroThreshold(self):
            """returns the threshold to use to check for zeros
            """
            return None
        
    class SolverSolution(SolverBase):
        """Contains equations for evaluating one unknown variable. The variable can have multiple solutions, and the solution is only valid if every equation in checkforzeros is non-zero
        """
        jointname = None
        jointeval = None
        jointevalcos = None
        jointevalsin = None
        AddPiIfNegativeEq = None
        isHinge = True
        checkforzeros = None
        thresh = None
        AddHalfTanValue = False
        dictequations = None
        presetcheckforzeros = None
        equationsused = None
        """Meaning of FeasibleIsZeros:
        If set to false, then solution is feasible only if all of these equations evalute to non-zero.
        If set to true, solution is feasible only if all these equations evaluate to zero.
        """
        FeasibleIsZeros = False
        score = None
        def __init__(self, jointname, jointeval=None,jointevalcos=None,jointevalsin=None,AddPiIfNegativeEq=None,isHinge=True,thresh=0.000001):

            self.jointname = jointname
            self.jointeval = jointeval
            self.jointevalcos = jointevalcos
            self.jointevalsin = jointevalsin
            self.AddPiIfNegativeEq = AddPiIfNegativeEq
            self.isHinge=isHinge
            self.thresh = thresh
            self.presetcheckforzeros = []
            self.dictequations = []
            self.equationsused = []
            assert(self.checkValidSolution())
        def subs(self,solsubs):
            if self.jointeval is not None:
                self.jointeval = [e.subs(solsubs) for e in self.jointeval]
            if self.jointevalcos is not None:
                self.jointevalcos = [e.subs(solsubs) for e in self.jointevalcos]
            if self.jointevalsin is not None:
                self.jointevalsin = [e.subs(solsubs) for e in self.jointevalsin]
            if self.checkforzeros is not None:
                self.checkforzeros = [e.subs(solsubs) for e in self.checkforzeros]
            self.dictequations = [(s,v.subs(solsubs)) for s,v in self.dictequations]
            self.presetcheckforzeros = [e.subs(solsubs) for e in self.presetcheckforzeros]
            self.equationsused = [e.subs(solsubs) for e in self.equationsused]
            if not self.checkValidSolution():
                raise IKFastSolver.CannotSolveError('substitution produced invalid results')
            return self
        def generate(self, generator):
            assert(self.checkValidSolution())
            return generator.generateSolution(self)
        def end(self, generator):
            return generator.endSolution(self)
        def numsolutions(self):
            n=0
            if self.jointeval is not None:
                n += len(self.jointeval)
            if self.jointevalcos is not None:
                n += 2*len(self.jointevalcos)
            if self.jointevalsin is not None:
                n += 2*len(self.jointevalsin)
            return n
        def checkValidSolution(self):
            from ikfast_IKFastSolver import IKFastSolver
            valid=True
#            exec(ipython_str)
            if self.jointeval is not None:
                valid &= all([IKFastSolver.isValidSolution(e) for e in self.jointeval])
            if self.jointevalsin is not None:
                valid &= all([IKFastSolver.isValidSolution(e) for e in self.jointevalsin])
            if self.jointevalcos is not None:
                valid &= all([IKFastSolver.isValidSolution(e) for e in self.jointevalcos])
            return valid
        def getPresetCheckForZeros(self):
            return self.presetcheckforzeros
        def getEquationsUsed(self):
            return self.equationsused
        def GetZeroThreshold(self):
            return self.thresh
        
    class SolverPolynomialRoots(SolverBase):
        """find all roots of the polynomial and plug it into jointeval. poly should be Poly
        """
        jointname = None
        poly = None
        polybackup = None # if poly does not yield and results, this polynomial will be solved instead
        jointeval = None
        jointevalcos = None # not used
        jointevalsin = None # not used
        checkforzeros = None
        postcheckforzeros = None # fail if any zero
        postcheckfornonzeros = None # fail if any nonzero
        postcheckforNumDenom = None # list of (A,B) pairs where Ax=B was used. Fail if A==0&&B!=0
        postcheckforrange = None # checks that value is within [-1,1]
        dictequations = None
        postcheckforzerosThresh = 1e-8 # threshold for checking postcheckforzeros. if abs(val) <= postcheckforzerosThresh: skip
        postcheckfornonzerosThresh = 1e-8 # threshold for checking postcheckfornonzeros. if abs(val) > postcheckfornonzerosThresh: skip
        postcheckforrangeThresh = 1e-8 # threshold for checking postcheckforrange. if  val <= -1-postcheckforrangeThresh || val > 1+postcheckforrangeThresh: skip
        postcheckforNumDenomThresh = 1e-8 # threshold for checking postcheckforNumDenom: if abs(val[0]) <= postcheckforNumDenomThresh && abs(val[1]) > postcheckforNumDenomThresh: skip
        isHinge = True
        FeasibleIsZeros = False
        AddHalfTanValue = False
        score = None
        equationsused = None
        def __init__(self, jointname, poly=None, jointeval=None,isHinge=True):
            self.poly = poly
            assert(self.poly.degree(0)>0)
            self.jointname=jointname
            self.jointeval = jointeval
            self.isHinge = isHinge
            self.dictequations = []
            self.equationsused = []
        def numsolutions(self):
            if self.polybackup is None:
                return self.poly.degree(0)
            else:
                return max(self.poly.degree(0), self.polybackup.degree(0))
            
        def subs(self,solsubs):
            if self.jointeval is not None:
                self.jointeval = [e.subs(solsubs) for e in self.jointeval]
            if self.checkforzeros is not None:
                self.checkforzeros = [e.subs(solsubs) for e in self.checkforzeros]
            if self.postcheckforzeros is not None:
                self.postcheckforzeros = [e.subs(solsubs) for e in self.postcheckforzeros]
            if self.postcheckfornonzeros is not None:
                self.postcheckfornonzeros = [e.subs(solsubs) for e in self.postcheckfornonzeros]
            if self.postcheckforrange is not None:
                self.postcheckforrange = [e.subs(solsubs) for e in self.postcheckforrange]
            if self.postcheckforNumDenom is not None:
                self.postcheckforNumDenom = [e.subs(solsubs) for e in self.postcheckforNumDenom]
            self.dictequations = [(s,v.subs(solsubs)) for s,v in self.dictequations]
            self.equationsused = [e.subs(solsubs) for e in self.equationsused]
            if self.poly is not None:
                self.poly = Poly(self.poly.subs(solsubs),*self.poly.gens)
            if self.polybackup is not None:
                self.polybackup = Poly(self.polybackup.subs(solsubs),*self.polybackup.gens)
            assert(self.checkValidSolution())
            return self
        def generate(self, generator):
            return generator.generatePolynomialRoots(self)
        def end(self, generator):
            return generator.endPolynomialRoots(self)
        def checkValidSolution(self):
            from ikfast_IKFastSolver import IKFastSolver
            valid = True
            if self.poly is not None:
                valid &= IKFastSolver.isValidSolution(self.poly.as_expr())
            if self.polybackup is not None:
                valid &= IKFastSolver.isValidSolution(self.polybackup.as_expr())
            if self.jointeval is not None:
                valid &= all([IKFastSolver.isValidSolution(e) for e in self.jointeval])
            return valid
        def getPresetCheckForZeros(self):
            # make sure that all the coefficients containing higher-order variables are not 0
            zeroeq = S.Zero
            for monom, coeff in self.poly.terms():
                if monom[0] > 0:
                    if len(self.dictequations) > 0: # bug with sympy?
                        zeroeq += abs(coeff.subs(self.dictequations))
                    else:
                        zeroeq += abs(coeff)
            if self.polybackup is not None:
                for monom, coeff in self.polybackup.terms():
                    if monom[0] > 0:
                        if len(self.dictequations) > 0: # bug with sympy?
                            zeroeq += abs(coeff.subs(self.dictequations))
                        else:
                            zeroeq += abs(coeff)
            return [zeroeq]#self.poly.LC()]
        def getEquationsUsed(self):
            return self.equationsused
        def GetZeroThreshold(self):
            return self.postcheckforzerosThresh # not really sure...
        
    class SolverCoeffFunction(SolverBase):
        """Evaluate a set of coefficients and pass them to a custom function which will then return all possible values of the specified variables in jointnames.
        """
        jointnames = None
        jointeval = None
        isHinges = True
        exportvar = None
        exportcoeffeqs = None
        rootmaxdim = None
        exportfnname = None
        jointevalcos = None # used for half angles
        jointevalsin = None # used for half angles
        checkforzeros = None
        FeasibleIsZeros = False
        score = None
        presetcheckforzeros = None
        dictequations = None
        equationsused = None
        def __init__(self, jointnames, jointeval=None, exportvar=None, exportcoeffeqs=None,exportfnname=None,isHinges=None,rootmaxdim=16,jointevalcos=None,jointevalsin=None):
            self.jointnames=jointnames
            self.jointeval = jointeval
            self.isHinges = isHinges
            self.exportvar=exportvar
            self.exportcoeffeqs=exportcoeffeqs
            self.exportfnname=exportfnname
            self.rootmaxdim=rootmaxdim
            self.jointevalsin=jointevalsin
            self.jointevalcos=jointevalcos
            self.presetcheckforzeros = []
            self.dictequations = []
            self.equationsused = []
        def numsolutions(self):
            return self.rootmaxdim
        def subs(self,solsubs):
            if self.jointeval is not None:
                self.jointeval = [e.subs(solsubs) for e in self.jointeval]
            if self.jointevalcos is not None:
                self.jointevalcos = [e.subs(solsubs) for e in self.jointevalcos]
            if self.jointevalsin is not None:
                self.jointevalsin = [e.subs(solsubs) for e in self.jointevalsin]
            if self.checkforzeros is not None:
                self.checkforzeros = [e.subs(solsubs) for e in self.checkforzeros]
            self.dictequations = [(s,v.subs(solsubs)) for s,v in self.dictequations]
            self.presetcheckforzeros = [e.subs(solsubs) for e in self.presetcheckforzeros]
            self.equationsused = [e.subs(solsubs) for e in self.equationsused]
            #if self.poly is not None:
            #    self.poly = Poly(self.poly.subs(solsubs)...)
            assert(self.checkValidSolution())
            return self
        def generate(self, generator):
            return generator.generateCoeffFunction(self)
        def end(self, generator):
            return generator.endCoeffFunction(self)
        def checkValidSolution(self):
            #if self.poly is not None:
            #    valid = IKFastSolver.isValidSolution(self.poly.as_expr())
            if self.jointeval is not None:
                valid &= all([IKFastSolver.isValidSolution(e) for e in self.jointeval])
            if self.jointevalcos is not None:
                valid &= all([IKFastSolver.isValidSolution(e) for e in self.jointevalcos])
            if self.jointevalsin is not None:
                valid &= all([IKFastSolver.isValidSolution(e) for e in self.jointevalsin])
            return valid
        def getPresetCheckForZeros(self):
            return self.presetcheckforzeros
        def getEquationsUsed(self):
            return self.equationsused

    class SolverMatrixInverse(SolverBase):
        """Take the inverse of a large matirx and set the coefficients of the inverse to the symbols in Asymbols.
        """
        A = None
        Asymbols = None # has to be same size as B
        checkforzeros = None
        def __init__(self, A, Asymbols):
            self.A = A
            self.Asymbols = Asymbols
        def subs(self,solsubs):
            return self
        def generate(self, generator):
            return generator.generateMatrixInverse(self)
        def end(self, generator):
            return generator.endMatrixInverse(self)
        def checkValidSolution(self):
            return True
        def getsubs(self,psubs):
            Asub = self.A.subs(psubs)
            d = Asub.det()
            if d == S.Zero:
                raise IKFastSolver.CannotSolveError('determinant for matrix is zero')
            
            Anew = Asub.inv()
            subs = []
            for i in range(self.A.shape[0]):
                for j in range(self.A.shape[1]):
                    if self.Asymbols[i][j] is not None:
                        subs.append((self.Asymbols[i][j],Anew[i,j]))
            return subs

    class SolverConditionedSolution(SolverBase):
        """set solutions based on evaluating equations
        """
        dictequations = None
        solversolutions = None # a list of solutions. If the solution's checkforzeros evaluates to all zeros, then that solution us used
        thresh=0.000001
        def __init__(self, solversolutions):
            self.solversolutions = solversolutions
            self.dictequations = []
        def subs(self,solsubs):
            for s in self.solversolutions:
                s.subs(solsubs)
            return self
        def generate(self, generator):
            return generator.generateConditionedSolution(self)
        def end(self, generator):
            return generator.endConditionedSolution(self)
        def GetChildrenOfType(self, classinstance):
            nodes = []
            for childnode in self.solversolutions:
                if isinstance(childnode, classinstance):
                    nodes.append(childnode)
                nodes += childnode.GetChildrenOfType(classinstance)
            return nodes
        def GetZeroThreshold(self):
            return self.thresh

    class SolverBranchConds(SolverBase):
        """
        take certain branches depending if a set of equations evaluate to zero.
        Each branch can also have dictequations
        """
        jointbranches = None # list of (checkzeroequations, branch, dictequations)
        thresh = 0.000005 # because it is && comparison, have to relax the threshold otherwise SolverCheckZeros checks that check with same thresh can fail on either equation not being 0, and this only passes when all equations are 0. so have to set threshold higher than what is used for SolverCheckZeros
        def __init__(self, jointbranches):
            self.jointbranches = jointbranches
        def generate(self, generator):
            return generator.generateBranchConds(self)
        def end(self, generator):
            return generator.endBranchConds(self)
        def GetChildrenOfType(self, classinstance):
            nodes = []
            for checkzeroequations, branch, extradictequations in self.jointbranches:
                for childnode in branch:
                    if isinstance(childnode, classinstance):
                        nodes.append(childnode)
                    nodes += childnode.GetChildrenOfType(classinstance)
            return nodes
        def GetZeroThreshold(self):
            return self.thresh

    class SolverCheckZeros(SolverBase):
        jointname = None
        jointcheckeqs = None # only used for evaluation
        zerobranch = None
        nonzerobranch = None
        anycondition=None
        dictequations=None
        thresh=None # a threshold of 1e-6 breaks hiro ik
        equationsused = None
        def __init__(self, jointname, jointcheckeqs, zerobranch, nonzerobranch,thresh=None,anycondition=True):
            self.jointname = jointname
            self.jointcheckeqs = jointcheckeqs
            self.zerobranch = zerobranch
            self.nonzerobranch = nonzerobranch
            if thresh is None:
                self.thresh = 0.000001
            else:
                self.thresh = thresh
            self.anycondition = anycondition
            self.dictequations = []
        def generate(self, generator):
            return generator.generateCheckZeros(self)
        def end(self, generator):
            return generator.endCheckZeros(self)
        def getPresetCheckForZeros(self):
            return []
        def checkValidSolution(self):
            for branch in self.nonzerobranch:
                if not branch.checkValidSolution():
                    return False
            for branch in self.zerobranch:
                if not branch.checkValidSolution():
                    return False
            return True
        def numsolutions(self):
            return 1
        def subs(self,solsubs):
            for branch in self.nonzerobranch:
                if hasattr(branch,'subs'):
                    branch.subs(solsubs)
            for branch in self.zerobranch:
                if hasattr(branch,'subs'):
                    branch.subs(solsubs)
            return self
        def getEquationsUsed(self):
            return self.equationsused
        def GetChildrenOfType(self, classinstance):
            nodes = []
            for childnode in self.nonzerobranch + self.zerobranch:
                if isinstance(childnode, classinstance):
                    nodes.append(childnode)
                nodes += childnode.GetChildrenOfType(classinstance)
            return nodes
        def GetZeroThreshold(self):
            return self.thresh
    class SolverFreeParameter(SolverBase):
        jointname = None
        jointtree = None
        def __init__(self, jointname, jointtree):
            self.jointname = jointname
            self.jointtree = jointtree
        def generate(self, generator):
            return generator.generateFreeParameter(self)
        def end(self, generator):
            return generator.endFreeParameter(self)
        def GetChildrenOfType(self, classinstance):
            nodes = []
            for childnode in self.jointtree:
                if isinstance(childnode, classinstance):
                    nodes.append(childnode)
                nodes += childnode.GetChildrenOfType(classinstance)
            return nodes
        
    class SolverRotation(SolverBase):
        T = None
        jointtree = None
        functionid=0
        def __init__(self, T, jointtree):
            self.T = T
            self.jointtree = jointtree
            self.dictequations = []
        def generate(self, generator):
            return generator.generateRotation(self)
        def end(self, generator):
            return generator.endRotation(self)

    class SolverFunction(SolverBase):
        jointtree = None
        name='innerfn'
        def __init__(self, name, jointtree):
            self.name = name
            self.jointtree = jointtree
            self.dictequations = []
        def generate(self, generator):
            return generator.generateFunction(self)
        def end(self, generator):
            return generator.endFunction(self)
        def GetChildrenOfType(self, classinstance):
            nodes = []
            for childnode in self.jointtree:
                if isinstance(childnode, classinstance):
                    nodes.append(childnode)
                nodes += childnode.GetChildrenOfType(classinstance)
            return nodes
        
    class SolverStoreSolution(SolverBase):
        """Called when all the unknowns have been solved to add a solution.
        """
        alljointvars = None
        checkgreaterzero = None # used for final sanity checks to ensure IK solution is consistent
        thresh = 0
        offsetvalues = None
        isHinge = None
        def __init__(self, alljointvars,checkgreaterzero=None,isHinge=None):
            self.alljointvars = alljointvars
            self.checkgreaterzero = checkgreaterzero
            self.isHinge=isHinge
            if isHinge is None:
                log.warn('SolverStoreSolution.isHinge is not initialized')
                self.isHinge = [True]*len(self.alljointvars)
        def generate(self, generator):
            return generator.generateStoreSolution(self)
        def end(self, generator):
            return generator.endStoreSolution(self)
        
    class SolverSequence(SolverBase):
        jointtrees = None
        def __init__(self, jointtrees):
            self.jointtrees = jointtrees
        def generate(self, generator):
            return generator.generateSequence(self)
        def end(self, generator):
            return generator.endSequence(self)
        def GetChildrenOfType(self, classinstance):
            nodes = []
            for tree in self.jointtrees:
                for childnode in tree:
                    if isinstance(childnode, classinstance):
                        nodes.append(childnode)
                    nodes += childnode.GetChildrenOfType(classinstance)
            return nodes
        
    class SolverBreak(SolverBase):
        """Terminates this scope"""
        comment = None # a comment for the reason of the break
        varsubs = None # variable substitutions that were valid at the break
        othersolvedvars = None # the solved variables already
        solsubs = None # the substitutions of the solved variables
        endbranchtree = None # a node that points to the end of the tree
        def __init__(self, comment, varsubs=list(), othersolvedvars=list(), solsubs=list(), globalsymbols=list(), endbranchtree=None):
            self.comment = comment
            self.varsubs = list(varsubs)
            self.othersolvedvars = list(othersolvedvars)
            self.solsubs = list(solsubs)
            self.endbranchtree = endbranchtree
        def generate(self,generator):
            return generator.generateBreak(self)
        def end(self,generator):
            return generator.endBreak(self)
        def checkValidSolution(self):
            return True
        
    class SolverIKChainTransform6D(SolverBase):
        solvejointvars = None
        freejointvars = None
        jointtree = None
        Tfk = None
        Tee = None
        dictequations = None
        def __init__(self, solvejointvars, freejointvars, Tee, jointtree,Tfk=None):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Tee = Tee
            self.jointtree = jointtree
            self.Tfk = Tfk
            self.dictequations = []
        def generate(self, generator):
            return generator.generateChain(self)
        def end(self, generator):
            return generator.endChain(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Tfk = Tleft*self.Tfk
            self.Tee = Tleftinv*self.Tee

    class SolverIKChainRotation3D(SolverBase):
        solvejointvars = None
        freejointvars = None
        Rfk = None
        Ree = None
        jointtree = None
        dictequations = None
        def __init__(self, solvejointvars, freejointvars, Ree, jointtree,Rfk=None):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Ree = Ree
            self.Rfk=Rfk
            self.jointtree = jointtree
            self.dictequations = []
        def generate(self, generator):
            return generator.generateIKChainRotation3D(self)
        def end(self, generator):
            return generator.endIKChainRotation3D(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Rfk = Tleft[0:3,0:3]*self.Rfk
            self.Ree = Tleftinv[0:3,0:3]*self.Ree

    class SolverIKChainTranslation3D(SolverBase):
        solvejointvars = None
        freejointvars = None
        jointtree = None
        Pfk = None
        Pee = None
        dictequations = None
        uselocaltrans = False
        def __init__(self, solvejointvars, freejointvars, Pee, jointtree,Pfk=None):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Pee = Pee
            self.jointtree = jointtree
            self.Pfk=Pfk
            self.dictequations = []
        def generate(self, generator):
            return generator.generateIKChainTranslation3D(self)
        def end(self, generator):
            return generator.endIKChainTranslation3D(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Pfk = Tleft[0:3,0:3]*self.Pfk+Tleft[0:3,3]
            self.Pee = Tleftinv[0:3,0:3]*self.Pee+Tleftinv[0:3,3]

    class SolverIKChainTranslationXY2D(SolverBase):
        solvejointvars = None
        freejointvars = None
        jointtree = None
        Pfk = None
        Pee = None
        dictequations = None
        def __init__(self, solvejointvars, freejointvars, Pee, jointtree,Pfk=None):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Pee = Pee
            self.jointtree = jointtree
            self.Pfk=Pfk
            self.dictequations = []
        def generate(self, generator):
            return generator.generateIKChainTranslationXY2D(self)
        def end(self, generator):
            return generator.endIKChainTranslationXY2D(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Pfk = Tleft[0:2,0:2]*self.Pfk+Tleft[0:2,3]
            self.Pee = Tleftinv[0:2,0:2]*self.Pee+Tleftinv[0:2,3]
            
    class SolverIKChainDirection3D(SolverBase):
        solvejointvars = None
        freejointvars = None
        jointtree = None
        Dfk = None
        Dee = None
        dictequations = None
        def __init__(self, solvejointvars, freejointvars, Dee, jointtree,Dfk=None):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Dee = Dee
            self.jointtree = jointtree
            self.Dfk=Dfk
            self.dictequations = []
        def generate(self, generator):
            return generator.generateIKChainDirection3D(self)
        def end(self, generator):
            return generator.endIKChainDirection3D(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Dfk = Tleft[0:3,0:3]*self.Dfk
            self.Dee = Tleftinv[0:3,0:3]*self.Dee

    class SolverIKChainRay(SolverBase):
        solvejointvars = None
        freejointvars = None
        jointtree = None
        Pfk = None
        Dfk = None
        Pee = None
        Dee = None
        dictequations = None
        is5dray = False # if True, then full 3D position becomes important and things shouldn't be normalized
        def __init__(self, solvejointvars, freejointvars, Pee, Dee, jointtree,Pfk=None,Dfk=None,is5dray=False):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Pee = Pee
            self.Dee = Dee
            self.jointtree = jointtree
            self.Pfk = Pfk
            self.Dfk = Dfk
            self.dictequations = []
            self.is5dray=is5dray
        def generate(self, generator):
            return generator.generateIKChainRay(self)
        def end(self, generator):
            return generator.endIKChainRay(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Pfk = Tleft[0:3,0:3]*self.Pfk+Tleft[0:3,3]
            self.Dfk = Tleft[0:3,0:3]*self.Dfk
            self.Pee = Tleftinv[0:3,0:3]*self.Pee+Tleftinv[0:3,3]
            self.Dee = Tleftinv[0:3,0:3]*self.Dee

    class SolverIKChainLookat3D(SolverBase):
        solvejointvars = None
        freejointvars = None
        jointtree = None
        Pfk = None
        Dfk = None
        Pee = None
        dictequations = None
        def __init__(self, solvejointvars, freejointvars, Pee, jointtree,Pfk=None,Dfk=None):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Pee = Pee
            self.jointtree = jointtree
            self.Pfk=Pfk
            self.Dfk=Dfk
            self.dictequations = []
        def generate(self, generator):
            return generator.generateIKChainLookat3D(self)
        def end(self, generator):
            return generator.endIKChainLookat3D(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Pfk = Tleft[0:3,0:3]*self.Pfk+Tleft[0:3,3]
            self.Dfk = Tleft[0:3,0:3]*self.Dfk
            self.Pee = Tleftinv[0:3,0:3]*self.Pee+Tleftinv[0:3,3]
            
    class SolverIKChainAxisAngle(SolverBase):
        solvejointvars = None
        freejointvars = None
        jointtree = None
        Pfk = None
        Pee = None
        dictequations = None
        angleee=None
        anglefk=None
        iktype=None
        def __init__(self, solvejointvars, freejointvars, Pee, angleee,jointtree,Pfk=None,anglefk=None,iktype=None):
            self.solvejointvars = solvejointvars
            self.freejointvars = freejointvars
            self.Pee = Pee
            self.anglefk=anglefk
            self.jointtree = jointtree
            self.Pfk=Pfk
            self.angleee=angleee
            self.dictequations = []
            self.iktype=iktype
        def generate(self, generator):
            return generator.generateSolverIKChainAxisAngle(self)
        def end(self, generator):
            return generator.endSolverIKChainAxisAngle(self)
        def leftmultiply(self,Tleft,Tleftinv):
            self.Pfk = Tleft[0:2,0:2]*self.Pfk+Tleft[0:2,3]
            self.Pee = Tleftinv[0:2,0:2]*self.Pee+Tleftinv[0:2,3]
            assert(0) # need to change angle
