#!/usr/bin/env python3
#
# SpecialRelativity - Handle the coordinate systems and transformations of Special Relativity.  Working toward laying out
# Local Lorentz Frames for the twin (and triplet) problem, even if the LLF are not global.
#
# Bill Gabella 20181113
#

import numpy as np
import matplotlib.pyplot as plt

def rot(theta):
    """Rotate around the z-axis by theta radians.  Matrix traditionally multiplies [x,y]^T vector.
    """
    cos = np.cos(theta)
    sin = np.sin(theta)
    return( np.array( [[cos, sin], [-sin, cos]] ) )
           
def boost(beta):
    """Return the matrix for the boost from lab (S) frame to co-moving along x-axis frame (S') moving to the
    right (positive x) with velocity beta = v/c.
    """
    if abs(beta) >= 1:
        print('***ERROR in SpecialRelativity.boost, beta is {:.3f} .'.format(beta) )
        return( np.array( [ [1.0, 0], [0, 1.0] ] ) )
    gamma = 1.0/np.sqrt( 1.0-beta*beta )
    ch = gamma
    sh = gamma*beta
    return( np.array( [ [ch, -sh], [-sh, ch] ] ) )    
       

class txgrid():
    """Build the Special Relativity coord system, usually called S' with t' (vertical) and x' (horizontal).
    Positive beta means that S' moves to the right (pos x) in S.  Negative means is moves to the left (neg x).
    """
    def __init__(self, beta):
        if np.abs(beta) < 1:
            self.beta = beta
        else:
            self.beta = 0.0
            print('txgrid:  **Error**: beta not between -1 and +1, is {}.'.format(beta) )
            print('txgrid:             Setting beta to 0.0.')
        self.eta = np.arctanh(beta)  # Rapidity eta where beta = tanh eta.
        self.psi = np.arctan(1/beta) # rads, Angle between x and x' axis.  Positive is x' above x, negative x' below x.
        self.gamma = 1.0/np.sqrt(1.0-self.beta*self.beta)
        #self.axgrid = self.grid()
        
    def setbeta(self, beta):
        self.beta=beta
    def getbeta(self, beta):
        return( self.beta )
    
    def boost(self):
        """Multiplies vector (t,x) to give (t', x').  Use inverse to go the other way.
        """
        ch = self.gamma
        sh = self.gamma*self.beta
        return( np.array( [ [ch, -sh], [-sh, ch] ] ) )

    # Maybe https://matplotlib.org/devel/index.html is better way to go.
    def grid(self):
        """Create the Matplotlib Pyplot graph for the S' grid, assumed on the rectangulare (t,x) grid.
        """
        def startEndPoints(start, end, num):
            """Generate num points between a start and stop point.  Assumed to be (t,x) points.
            """
            ll = np.linspace(0,1,num)
            xxs = start[0]*(1-ll)+end[0]*ll
            tts = start[1]*(1-ll)+end[1]*ll
            return( np.array([xxs, tts]) )
        
        fig, ax = plt.subplots( figsize=(10,8) )
        # The slope is beta.
        # x' axis, input points in (x,t) S points
        boost = self.boost()
        invboost = np.linalg.inv( boost )
        astart = np.dot(invboost, [0, 0] )
        aend = np.dot(invboost, [0,10] )
        aaxis = startEndPoints( astart, aend, 40 )

        ax.plot(aaxis[0], aaxis[1], 'b-')
        
        # t' axis
        astart = np.dot(invboost, [0, 0])
        aend = np.dot(invboost, [10, 0])
        aaxis = startEndPoints( astart, aend, 40 )

        ax.plot(aaxis[0], aaxis[1], 'b-')
        
        # light lines on all 45 degrees
        aaxis = startEndPoints( [-10, -10], [10, 10], 40 )
        ax.plot(aaxis[0], aaxis[1], 'k--')
        aaxis = startEndPoints( [-10, 10], [10, -10], 40 )
        ax.plot(aaxis[0], aaxis[1], 'k--')
        
        # t-hat and x-hat axes
        aaxis = startEndPoints( [0, 0], [10, 0], 40 )
        ax.plot(aaxis[0], aaxis[1], 'k-')
        aaxis = startEndPoints( [0, 0], [0, 10], 40 )
        ax.plot(aaxis[0], aaxis[1], 'k-')
        ax.grid(True)
        
        return( (fig, ax) )

    

    
    
def main():
    print('This is SpecialRelativity.py module.')
    
if __name__ == '__main__':
    main()
    
