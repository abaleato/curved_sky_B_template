import numpy as np
import shts
import math

class First_order_template(object):
    """ First-order lensing template inspired by the implementation of the lensing quadratic estimators
        in appendix (A.2) of https://arxiv.org/pdf/1502.01591.pdf
        
    The estimator, q^{XY}(L), can be run on fields \bar{X} and \bar{Y} as

    q^{XY}(L) = 1/2 \int{d^2 l_X} \int{d^2 l_Y}
                    W^{XY}(l_X, l_Y, L) \bar{X}(l_X) \bar{Y}(l_Y)

    with l_X + l_Y = L.

    the weight function W^{s, XY} must be separable. For
    full-sky calculations it is encoded as

    W^{XY} = \sum_{i=0}^{N_i} \int{d^2 n}
                    {}_s^{i,X}Y_{l_X m_X}(n) W^{i,X}(l_X) *
                     {}_s^{i,Y}Y_{l_Y m_Y}(n) W^{i,Y}(l_Y) *
                      {}_s^{i,L}Y_{  L M  }(n) W_^{i,L}( L ).

    the spins s^{i,n} are stored in an array self.s[i][n] and the
    weights w^{i,n}(l) are encapsulated as functions w[i][n](l). some
    flat-sky-specific weight functions may also take lx, ly as arguments.
    """
    def __init__(self):
        pass

    def eval( self, barX, barY=None, **kwargs ):
        if barY is None:
            barY = barX

        if False:
            pass
        else:
            return self.eval_fullsky( barX, barY, **kwargs )

    def eval_fullsky( self, barX, barY ):
        """
        inputs:
             * barX            = input field \bar{X}. harmonic modes for the X field
                                 (in the 'vlm' indexing scheme, see shts.util).
             * barY            = input field \bar{Y}. harmonic modes for the X field
                                 (in the 'vlm' indexing scheme, see shts.util).
        """
        lmax   = self.lmax
        lmax_X = shts.util.nlm2lmax( len(barX) )
        lmax_Y = shts.util.nlm2lmax( len(barY) )

        nphi   = lmax_X+lmax_Y+lmax+1
        glq    = math.wignerd.gauss_legendre_quadrature( (lmax_X + lmax_Y + lmax)/2 + 1 )
        tht    = np.arccos(glq.zvec)
        phi    = np.linspace(0., 2.*np.pi, nphi, endpoint=False)

        ret = np.zeros( (lmax+1)**2, dtype=np.complex )
        for i in xrange(0, self.ntrm):
            # l_X term
            vlx = shts.util.alm2vlm( barX )
            for l in xrange(0, lmax_X+1):
                vlx[l**2:(l+1)**2] *= self.get_wlX(i,l)
            vmx = shts.vlm2map( self.get_slX(i), tht, phi, vlx )
            del vlx

            # l_Y term
            vly = shts.util.alm2vlm( barY )
            for l in xrange(0, lmax_X+1):
                vly[l**2:(l+1)**2] *= self.get_wlY(i,l)
            vmy = shts.vlm2map( self.get_slY(i), tht, phi, vly )
            del vly

            # multiply in position space
            vmm  = vmx * vmy
            del vmx, vmy

            # apply weights for harmonic integration
            for j, w in enumerate(glq.wvec):
                vmm[j,:] *= w * (2.*np.pi / nphi)

            # perform integration
            vlm = shts.map2vlm(lmax, self.get_slL(i), tht, phi, vmm)
            del vmm

            for l in xrange(0, lmax+1):
                vlm[l**2:(l+1)**2] *= 0.5*self.get_wlL(i,l)

            ret += vlm

        return ret

    def get_slX(self, i):
        return self.sl[i][0]

    def get_slY(self, i):
        return self.sl[i][1]

    def get_slL(self, i):
        return self.sl[i][2]

    def get_wlX(self, i, l, **kwargs):
        return self.wl[i][0](l, **kwargs)

    def get_wlY(self, i, l, **kwargs):
        return self.wl[i][1](l, **kwargs)

    def get_wlL(self, i, l, **kwargs):
        return self.wl[i][2](l, **kwargs)
