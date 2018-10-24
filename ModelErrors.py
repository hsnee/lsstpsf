import numpy as np
import matplotlib.pyplot as plt
from angles import r2d, r2arcs, d2arcs
import seaborn as sns;sns.set_style('darkgrid')
import lsst.sims.maf.stackers as stackers
import treecorr as tr, healpy as hp
%config InlineBackend.figure_format = 'retina'
matplotlib.rcParams['figure.figsize'] = 18,12
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.db as db
from lsst.sims.maf.plots import PlotHandler
from collections import defaultdict
from scipy.stats import moment

class ModelErrors():

    def __init__(self,ModelType):
        self.fwhm = 0.7 #arcsec
        self.alpha = None
        self.ModelType = ModelType
        self.survey_area = (np.radians(200),np.radians(50)) # in radians
        stars_X = np.random.rand(int(1E4))*self.survey_area[0]
        stars_Y = np.random.rand(int(1E4))*self.survey_area[1] - self.survey_area[1]
        stars = np.array((stars_X,stars_Y))
        self.stars = stars.swapaxes(1,0)
        self.e = defaultdict(list)

    def PSF(self):
        pass

    def STAR(self):
        pass

    def getModels(self):
        '''calls one of the model methods to create an overly simplified model,
        or (N/A yet) imports one'''
        # construct a circle around the central position.

        for i in range(len(self.positions)):
            if self.ModelType == 'radial':
                self.radial_pattern(position_num=i)
            # elif self.ModelType == 'horizontal':
            #         self.de1[i], self.de2[i] = horizontal_pattern(position)

    def getPositions(self,sqlWhere,database='/Users/myhome/Downloads/minion_1016_sqlite.db'):
        '''uses OpSim to find all dithered positions given some constraint sqlWhere'''
        opsdb = db.OpsimDatabase(database)
        pos = opsdb.fetchMetricData(('ditheredRA', 'ditheredDec'), sqlconstraint=sqlWhere)
        pos = np.array(pos)
        if any(row[1] == x for row in pos):
            pos2 = zip(*pos)[1]
        if any(row[2] == x for row in pos):
            pos3 = zip(*pos)[2]
        pos = np.array((pos2,pos3))
        self.positions = pos.swapaxes(1,0)

    def avgOver(self):
        ''' should rewrite '''
        for i in self.e.keys(): self.e[i] = (np.mean(self.e[i][0]),np.mean(self.e[i][1]))

    def process(self,sqlWhere):
        ''' runs all analysis methods '''
        self.getPositions(sqlWhere)
        self.getModels()
        self.avgOver()
        self.getTraces()
        self.getRhos()
        self.rhos2errors()


    def getRhos(self):
        ''' method to get the rho statistics, needs a model, and traces for rhos 2 through 5.'''
        e1 = np.zeros(len(self.e.keys()))
        e2 = np.zeros(len(self.e.keys()))
        for j,i in enumerate(self.e.keys()):
            e1[j] = self.e[i][0]
            e2[j] = self.e[i][1]
        de = np.sqrt(e1**2 + e2**2)
        cat = tr.Catalog(k=de, ra=np.array(self.e.keys())[:,0], dec=np.array(self.e.keys())[:,1],ra_units='radians',dec_units='radians')
        dede = tr.KKCorrelation(min_sep=0.01, max_sep=50, nbins=100, sep_units='degrees')
        dede.process(cat)
        self.r = r = np.exp(dede.meanlogr)
        self.rho1 = dede.xi

        self.rho2 = None
        self.rho3 = None
        self.rho4 = None
        self.rho5 = None

    def getTraces(self):
        '''gets traces from sizes etc to be used in getting rho statistics'''
        self.T_PSF = None
        self.T_gal = None # this could be moved to init, but keeping it here for consistency

    def rhos2errors(self):
        '''propagates rho statistics into shear errors'''
        self.delta_xip = 2 * (self.delta_T_psf / self.T_gal) * self.xip \
                    + (self.T_PSF/self.T_gal)**2 * rho1                 \
                    - self.alpha * (self.T_PSF/self.T_gal) * self.rho2  \
                    + (self.T_PSF/self.T_gal)**2 * self.rho3            \
                    + (self.T_PSF/self.T_gal)**2 * self.rho4            \
                    - self.alpha * (self.T_PSF/self.T_gal) * self.rho5

    def radial_pattern(self,position_num):
        ''' method to create a radial pattern (one of the overly simplified models),
            e = 0.05* distance from origin.
            Should probably combine with horizontal pattern.'''
        from angles import r2d, r2arcs, d2arcs
        center = self.positions[position_num]
        stars = self.stars[np.where(((self.stars[:,0]-center[0])**2+(self.stars[:,1]-center[1])**2)<np.radians(1.2))[0]]
        try:
            rel_X, rel_Y = stars[:,0] - center[0], stars[:,1] - center[1]
        except:
            return
        r = np.sqrt(rel_X**2+rel_Y**2)
        theta = np.arctan((rel_Y)/(rel_X))
        e1 = r*np.cos(2*theta)/20
        e2 = r*np.sin(2*theta)/20
        for i in range(len(e1)):
            self.e[(stars[i][0],stars[i][1])].append((e1[i],e2[i]))

    def horizontal_pattern(self): ## untested for new algorithm
        ''' method to create a horizontal pattern (one of the overly simplified models)
            e = 0.1'''
        from angles import r2d, r2arcs, d2arcs
        e1 = np.ones(1000)*0.1
        e2 = np.zeros(len(self.stars_X))
        return e1,e2
