#!/usr/bin/env python3
#
#  The collection of Local Lorentz Frame functions, working with Justin Stevens.
#  20190129 weg  Seems to be working, the LLF() and LorShape() classes with the 
#                helper functions.  Used first with BugAndRivet.ipynb notebook.
#


import numpy as np

# Rotation in [x, y, z]^T .
def rotmat( eulerAngles ):
    """Rotation matrix with three Euler angles as defined in physics, i.e. Goldstein, Fig. 5.7 p 209.
    Order is rotate around z, then x, then z (new z) with angles phi, theta, and psi, resp.
    Assumes 3D.  Returns matrix.
    """
    (phi, th, psi) = eulerAngles
    cc = np.cos(phi)
    ss = np.sin(phi)
    mat1 = np.array( [ [cc, -ss, 0],
                        [ss, cc, 0],
                         [0, 0, 1] ] )
    cc = np.cos(th)
    ss = np.sin(th)
    mat2 = np.array( [])
    cc = np.cos(ph)
    ss = np.sin(phi)

# Adds the time coordinate in front of the ndim array.    
# def prependTime(mypoints, tvalue):
#     """Give a single time value, and it prepends that value to the x, y points as a numpy array.
#     """
#     xx = mypoints.transpose()[0]
#     yy = mypoints.transpose()[1]
#     tt = np.full(xx.shape, tvalue)  # Should create a numpy vector/array of the right size and shape.
    
#     final= np.array([tt, xx, yy]).transpose()
#     return(final)

def prependTime(mypoints, tvalue):
    """Give a single time value, and it prepends that value to the x, y, z, ... points as a numpy array.
    Useful to create a shape in x, y, z, ... and the create the four vectors for those points.
    """
    return( np.insert(mypoints, 0, tvalue, axis=1) )  
    # Use the fancy np.insert function.
    # See https://stackoverflow.com/questions/36998260/prepend-element-to-numpy-array
    
def translateShape(mypoints, tvector):
    """Translate each point by the vector amount tvector.
    """
    if len(mypoints[0]) != len(tvector):
           print('***Error, LLFTransforms, translateShape: tvector wrong length.')
           print('        len(mypoints[0]) is {} and len(tvector) is {}.'.format( len(mypoints[0]), len(tvector) ) )
           return(mypoints)
           
    return( np.add(mypoints, tvector) ) 
              
    
    
    
# Class for the Local Lorentz Frame, especially translating S' coords to S (lab) ones.
class LLF():
    """Gives some methods and data for a local observer, a Local Lorentz Frame,
    with a position and velocity in the original 'lab' frame.  
    For now, velocity is in the x-direction and y and z are perpendiculat to the motion.
    """
    def __init__(self, beta, ndim):
        #self.pos = np.array(pos)
        self.beta = beta # np.array(beta)
        self.ndim = ndim

        # Check dimensions are same.
        #self.checkDim()
        self.gamma = self.mygamma()
        # Calculate the t-x part of the matrix, mult (t,x)^T to get (t', x')^T .
        self.mat_tx, self.inv_mat_tx = self.txmatrix()
        #
        self.mat_all, self.inv_mat_all = self.allmatrix()
            
#     def getpos(self):
#         return( self.pos )
#     def setpos(self, pos):
#         self.pos = pos
    def getbeta(self):
        return( self.beta )
    def setbeta(self, beta): # Recalculate all the matrices.
        self.beta = beta
        self.gamma = self.mygamma()
        self.mat_tx, self.inv_mat_tx = self.txmatrix()
        self.mat_all, self.inv_mat_all = self.allmatrix()
        return(True)
    def getndim(self):
        return( self.ndim )
    def setndim(self, ndim):
        if ndim<=0:
            print('***Error, setndim: ndim is 0 or negative, {}, using 2'.format(ndim) )
            self.ndim = 2
            return(self.ndim)
        self.ndim = ndim
        return( True )
    
    def getmat_tx(self):
        return(self.mat_tx)

    def getmat_all(self):
        return(self.mat_all)
    
    def getinv_mat_tx(self):
        return(self.inv_mat_tx)
    
    def getinv_mat_all(self):
        return(self.inv_mat_all)
    
    
    def checkDim(self, pos):
        # Check dimensions are same.
        if len(pos) != self.ndim:
            print('***ERROR in class LLF, position has dims {} and ndim is {}'.format(len(pos), self.ndim ) )
            return(False)
        return(True)
            
    def mygamma(self):
        return( np.sqrt(1/(1-self.beta*self.beta) ) )
    
            
    def txmatrix(self):
        """Calculate just the t-x part of the matrix.  Recall that transverse dimensions are not changed.
        """
#         if self.beta == 0:
#             amat = np.identity(2)
#             self.mat_tx = amat
#             return(amat)
        
        gg = self.gamma
        bb = self.beta
        amat = np.array( [ [gg, -gg*bb], [-gg*bb, gg] ] )
        amatinv = np.linalg.inv(amat)
        return( (amat, amatinv) )
    
    def allmatrix(self):
        """Calculate the full t, x, y, z, ... matrix.
        """
        txmat = self.mat_tx
        # Other components are 1's along the diagonal.
        amat = np.identity(self.ndim)
        amat[0,0] = txmat[0,0]
        amat[0,1] = txmat[0,1]
        amat[1,0] = txmat[1,0]
        amat[1,1] = txmat[1,1]
        
        amatinv = np.linalg.inv(amat)
        return( (amat, amatinv) )
    
    def transTXY(self, pos):
        """Transform the S point [t, x, y]^T into the "moving" frame S'.
        """
        if not self.checkDim(pos):
            return( np.zeros(self.ndim) )
        ###print('{}\n{}'.format( self.mat_all, pos ) )
        return( np.dot( self.mat_all, pos) )
    
    def transTPXPYP(self, posprime):
        """Transform from S' point [t', x', y']^T into the "lab" frame S.
        Use the inverse matrix / transform."""
        if not self.checkDim(posprime):
            return( np.zeros(self.ndim) )
        return( np.dot( self.inv_mat_all, posprime) )
    
    def transTXYmany(self, arr):
        """Takes arr, a "vector" of points (t, x, y) and transforms them to (t', x', y').
        """
        aa = []
        for apos in arr:
            aa.append( self.transTXY(apos) )
        return( np.array(aa) )
        
    def transTPXPYPmany(self, arr):
        """Takes arr, a "vector" of points (t', x', y') and transforms them to (t, x, y).
        """
        aa = []
        for apos in arr:
            aa.append( self.transTPXPYP(apos) )
        return( np.array(aa) )

    

    
# Try out a "shape" class that inherits(??) a Lorentz transformation, a Lorentz frame from above.
class LorShape(LLF):
    """Assuming a shape in the x'-y' plane of S' that is fixed and just moving vertically in S'
    """
    def __init__(self, beta, ndim, xp0, mp):
        """beta is the speed, ndim is number of dimensions to work in, usually 2 or 3,
        xp0 is the numpy array of [t', x', y', ...] at start, u=0, and mp is the slope of these
        lines in the S' frame, typicall [1, 0, 0)], as the shape at rest in S' just moves vertically.
        numpoints is the number of points in the shape.  These will be treated like lines in S' and S.
        """
        LLF.__init__(self, beta, ndim)  # Initialize the Local Lorentz frame variables and methods.
        self.ndim = ndim # ???
        self.xp0 = xp0
        self.mp = mp
        self.numpoints = len(xp0)
        
        self.checkSize()
        
        self.x0 = self.calcX0()
        self.m = self.calcM()
        
    def checkSize(self):
        if len(self.xp0) != len(self.mp):
            print('***Error, LorShape, checkSize: xp0 and mp different sizes, {}, {}, resp.'.
                  format(len(self.xp0), len(self.mp) ) )
            return(False)
        if len(self.xp0[0]) != self.ndim:
            print('***Error, LorShape, checkSize: len(xp0[0]) and ndim different, {}, {}, resp.'.
                  format(len(self.xp0[0]), self.ndim) )
            return(False)
        return(True)
    
    def calcX0(self):
        return( self.transTPXPYPmany(self.xp0) )
    
    def calcM(self):
        return( self.transTPXPYPmany(self.mp) )
        
    def ufromt(self, ptindex, tt):
        """Given a point index on the set of shape points, ptindex, calculate the u value for 
        the given t (in S) value.  Returns u or 0.0 if there is a divide by 0 problem."""
        xx0 = self.x0[ptindex]
        mm = self.m[ptindex]
        if mm[0] != 0:
            uu = (tt - xx0[0])/mm[0]
            return(uu)
        else:
            print('***Error, LorShape, ufromt: mm[0] is zero! Divide by zero.')
            print('***                       : (t-x0[0] is {}, m[0] is {})'.format((tt-xx0[0]), mm[0]) )
            return(0.0)
        
    def ufromtp(self, ptindex, ttp):
        """Given a point index on the set of shape points, ptindex, calculate the u value for 
        the given t (in S) value.  Returns u or 0.0 if there is a divide by 0 problem."""
        xxp0 = self.xp0[ptindex]
        mmp = self.mp[ptindex]
        if mmp[0] != 0:
            uu = (ttp - xxp0[0])/mmp[0]
            return(uu)
        else:
            print('***Error, LorShape, ufromt: mm[0] is zero! Divide by zero.')
            print('***                       : (t-x0[0] is {}, m[0] is {})'.format((tt-xx0[0]), mm[0]) )
            return(0.0)
        
    def shapeXAtT(self, tt):
        """Return the shape array in S frame at the given time t (in S).
        """
        # Run over the points.
        xnew = np.zeros( (len(self.x0), self.ndim) )
        for ptindex, (aa, bb) in enumerate( zip(self.x0, self.m) ):
            uu = self.ufromt(ptindex, tt)
            xnew[ptindex] = aa + bb*uu
        return(xnew)
        
    def shapeXPAtTP(self, ttp):
        """Return the shape array in S frame at the given time t (in S).
        """
        # Run over the points.
        xpnew = np.zeros( (len(self.xp0), self.ndim) )
        for ptindex, (aa, bb) in enumerate( zip(self.xp0, self.mp) ):
            uu = self.ufromtp(ptindex, ttp)
            xpnew[ptindex] = aa + bb*uu
        return(xpnew)

    
    
def doTest():
    """Do some simple test items with LLF and LorShape.
    """
    
    aLLF = LLF(0.9, 3)
    print('aLLF = LLF(0.9, 3)')
    print('aLLF.beta is {} and aLLF.gamma is {}'.format( aLLF.beta, aLLF.gamma ) )
    print('aLLF.mat_all and aLLF.inv_mat_all are\n{}\n\n{}'.format( aLLF.mat_all, aLLF.inv_mat_all ) )
    print('...and their product is\n {}'.format( np.dot(aLLF.mat_all, aLLF.inv_mat_all) ) )

    print('\n')
    # Draw a closed polygon as a starting shape at t'=0, and will move by m' .
    # Make it 1 deep and +/- 1/4 for y.  All in the x-y plane and time is t' = 0.0 say.
    xdeep = 1.0
    ywidth = 0.25
    ysize = 2 # +/- the units in extent
    xback = 2.0
    # Append t later.
    aShape = np.array([ [0.0, ysize], [0.0, ywidth], [-xdeep, ywidth], [-xdeep, -ywidth], [0.0, -ywidth], 
                       [.0, -ysize], [-xback, -ysize], [-xback, ysize], [0.0, ysize]  ])

    aShapeStart = prependTime(aShape, 0.0)
    aSlopeStart = np.array( [ np.array([1.0, 0, 0]) for uu in aShapeStart ] )
    
    aLShape = LorShape(0.8, 3, aShapeStart, aSlopeStart)
    print('aLShape = LorShape(0.8, 3, aShapeStart, aSlopeStart)')
    
    print('aLShape.mat_all is (multiplies [t, x, y]^T to give the coordinates in S\'-frame)\n{}'.format( aLShape.mat_all ) ) 
    aa = aLShape.shapeXAtT(12.0)  # Find the new shape, given in the S'-frame in the S-frame, at time 12.0.
    print('aa = aLShape.shapeXAtT(12.0)')
    
    print('aa is shape in S-frame at t=12.0\n{}'.format(aa) )
    print('shape in S\'-frame at t\'=0 is\n{}'.format( aLShape.xp0 ) ) 
          

if __name__ == '__main__':
    
    doTest()
