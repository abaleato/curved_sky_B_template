import numpy as np
import integration_functions

class B_template_weights(integration_functions.First_order_template):
    """ By A.Baleato. Template for lensing B modes. All weights include  an additional - sign to adhere to Quicklens conventions. """
    def __init__(self, cls):
        """ initialize the B-template estimator.
            * cls = spectrum from which to recover lmax.
            """
        self.cls = cls
        self.lmax = len(cls)-1
        self.ntrm = 6
        
        self.wl = { i : {} for i in xrange(0, self.ntrm) }
        self.sl = { i : {} for i in xrange(0, self.ntrm) }
        
        self.wl[0][0] = self.wo_ml; self.sl[0][0] = 0 # This is w_l1 and s_1
        self.wl[0][1] = self.onehalf ; self.sl[0][1] = 2 # This is w_l2 and s_2
        self.wl[0][2] = self.plus_i; self.sl[0][2] = 2 # This is w_L and s
        
        self.wl[1][0] = self.onehalf; self.sl[1][0] = 0
        self.wl[1][1] = self.wo_ml; self.sl[1][1] = 2
        self.wl[1][2] = self.plus_i; self.sl[1][2] = 2
        
        self.wl[2][0] = self.minus_onehalf; self.sl[2][0] = 0
        self.wl[2][1] = self.plus_i; self.sl[2][1] = 2
        self.wl[2][2] = self.wo_ml; self.sl[2][2] = 2
        
        self.wl[3][0] = self.wo_ml; self.sl[3][0] = 0
        self.wl[3][1] = self.minus_onehalf ; self.sl[3][1] = -2
        self.wl[3][2] = self.plus_i; self.sl[3][2] = -2
        
        self.wl[4][0] = self.minus_onehalf; self.sl[4][0] = 0
        self.wl[4][1] = self.wo_ml; self.sl[4][1] = -2
        self.wl[4][2] = self.plus_i; self.sl[4][2] = -2
        
        self.wl[5][0] = self.onehalf; self.sl[5][0] = 0
        self.wl[5][1] = self.plus_i; self.sl[5][1] = -2
        self.wl[5][2] = self.wo_ml; self.sl[5][2] = -2
    
    def wo_ml(self, l, lx=None, ly=None):
        return l*(l+1.)
    def onehalf(self, l, lx=None, ly=None):
        return 0.5
    def minus_onehalf(self, l, lx=None, ly=None):
        return -0.5
    def plus_i(self, l, lx=None, ly=None):
        return 1.j


class T_template_weights(integration_functions.First_order_template):
    """ By A.Baleato. Template for 1st order lensing correction to temperature.
        All weights include  an additional - sign to adhere to Quicklens conventions. """

    def __init__(self, cls):
        """ initialize the T-template estimator.
            * cls = spectrum from which to recover lmax.
            """
        self.cls = cls
        self.lmax = len(cls) - 1
        self.ntrm = 3

        self.wl = {i: {} for i in xrange(0, self.ntrm)}
        self.sl = {i: {} for i in xrange(0, self.ntrm)}

        self.wl[0][0] = self.wo_ml;
        self.sl[0][0] = 0  # This is w_l1 and s_1
        self.wl[0][1] = self.minus_one;
        self.sl[0][1] = 0  # This is w_l2 and s_2
        self.wl[0][2] = self.one;
        self.sl[0][2] = 0  # This is w_L and s

        self.wl[1][0] = self.minus_one;
        self.sl[1][0] = 0
        self.wl[1][1] = self.wo_ml;
        self.sl[1][1] = 0
        self.wl[1][2] = self.one;
        self.sl[1][2] = 0

        self.wl[2][0] = self.one;
        self.sl[2][0] = 0
        self.wl[2][1] = self.one;
        self.sl[2][1] = 0
        self.wl[2][2] = self.wo_ml;
        self.sl[2][2] = 0

    def wo_ml(self, l, lx=None, ly=None):
        return l * (l + 1.)

    def one(self, l, lx=None, ly=None):
        return 1.

    def minus_one(self, l, lx=None, ly=None):
        return -1.
