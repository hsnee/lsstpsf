from collections import defaultdict
import numpy as np
from angles import arcs2r, arcs2d, d2arcs, d2r, r2d
import treecorr as tr
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.db as db
import random
from lsst.sims.utils import angularSeparation
from itertools import groupby


def arcm2r(theta):
    """adding this for consistency

    Args:
        theta (float): angle in arc minutes

    Returns:
        float: angle in radians
    """
    return arcs2r(theta*60)


class ModelErrors():

    """Summary
    This no longer supports OpsimV3 runs (e.g. the minion runs)

    Attributes:
        alpha (float): Description
        counter (defaultdict): Description
        delta_xip (TYPE): Description
        delta_xip_sigma (): Description
        DitherPattern (str): 'field' (none), or 'hexagonal', 'radial', 'spiral'
        DitherPatterns (list): all possible diher patterns
        FOVradius (float): field of view of LSST camera, 1.75 degrees
        fwhm (float): full width at half max
        ksstatistic (list): D-statistic for ks test
        Maker (str): person/organisation that made the database
        ModelType (str): 'radial' or 'horizontal'
        OpsimRun (str): name of opsim run database, without extension
        positions (ndarray): positions of objects to observe
        pvalues (list): Description
        r (ndarray): log(theta) for angles to calculate correlations at.
        rho1 (ndarray): rho statistic
        rho1_im (ndarray): rho statistic imaginary part
        rho1_sigma (ndarray): rho statistic std
        rho2 (ndarray): rho statistic
        rho2_im (ndarray): rho statistic
        rho2_sigma (ndarray): rho statistic std
        rho3 (ndarray): rho statistic
        rho3_im (ndarray): rho statistic
        rho3_sigma (ndarray): rho statistic std
        rho4 (ndarray): rho statistic
        rho4_im (ndarray): rho statistic
        rho4_sigma (ndarray): rho statistic std
        rho5 (ndarray): rho statistic
        rho5_im (ndarray): rho statistic
        rho5_sigma (ndarray): rho statistic std
        rotDitherPattern (TYPE): Description
        rotTelPos (ndarray): telescope angles to use for rotational dithering.
        runName (str): name of opsim run including directory, without ext
        savedStarsAngles (dict): keys representing star positions and values
                                 representing the angle distribution
                                 to FOV centers
        sigma (TYPE): Description
        size_error_ratio (float): Description
        Stacker (method): Description
        Stackers (dict): Description
        star_num (int): Description
        stars (TYPE): Description
        trace_ratio (float): Description
        TrM (float): 2nd ellipticity moment trace
        xip (ndarray): xi+
        year (str): year of survey
    """

    def __init__(self, ModelType, DitherPattern,
                 rotDithers, OpsimRun, objects_base,
                 year, random_seed=1000):
        """Summary

        Args:
            ModelType (str): residual model: 'radial' or 'horizontal'
            DitherPattern (str): usually 'random_visit', look for
                                 self.Stackers for a full list.
            rotDithers (bool): Whether to use rotational dithers or not.
            OpsimRun (str): the name of the OpSim run, without path or .db.
            objects_base (str): depth cut made based on 'Y10' or 'actual' year.
            year (int): considers the survey up to this year,
                        must be between [1, 10].

        Raises:
            ValueError: if objects base not in {'Y10', 'actual'}
            FileNotFoundError: if file containing star positions is not found
        """
        self.runName = '/global/cscratch1/sd/husni/OpsimRuns/'+OpsimRun
        self.OpsimRun = OpsimRun
        self.random_seed = random.seed(random_seed)
        self.xipList = []
        self.fwhm = 0.7  # arcsec
        self.sigma = self.fwhm/(2*np.sqrt(2*np.log(2)))
        self.TrM = 2*self.sigma**2
        self.trace_ratio = 2.8933775060156068  # calculated in GalSize.ipynb
        self.PSF.TrM = self.TrM
        self.STAR.TrM = self.TrM
        self.FOVradius = 1.75  # degrees
        self.alpha = 0.01
        self.ModelType = ModelType
        self.star_num = 50000
        self.bootstrap_iterations = 300
        self.size_error_ratio = 0.001

        if DitherPattern is 'alt_sched' \
                or DitherPattern is 'alt_sched_rolling' \
                or DitherPattern is 'altsched_riswap' \
                or DitherPattern is 'altsched_rolling_riswap':
            self.Maker = 'Daniel'
        elif DitherPattern is 'rolling_10yrs_opsim' \
                or DitherPattern is 'rolling_mis10yrs_opsim':
            self.Maker = 'Peter'
        else:
            self.Maker = 'OpSim'

        try:
            if objects_base == 'actual':
                stars_pos = np.load('newcutnpys/'+self.OpsimRun+str(
                                     objects_base)+'.npy')
            elif objects_base == 'Y10':
                stars_pos = np.load('newcutnpys/'+self.OpsimRun+'10.npy')
            else:
                raise ValueError('Cannot understand the base year, \
                                 choose either "Y10" or "actual"')
        except ValueError:
            try:
                stars_pos = np.load('objs/'+self.OpsimRun+'.npy')

            except FileNotFoundError:
                stars_pos = np.load('/global/cscratch1/sd/husni/OpsimRuns/'
                                    + self.OpsimRun + '.npy')

        self.stars = np.array(random.choices(list(stars_pos),
                                             self.random_seed,
                                             k=self.star_num))
        self.DELTA.e = defaultdict(np.array)
        self.DELTA.M = defaultdict(np.array)
        self.STAR.M = defaultdict(np.array)
        self.STAR.e = defaultdict(np.array)
        self.PSF.M = defaultdict(np.array)
        self.PSF.e = defaultdict(np.array)
        self.counter = defaultdict(int)
        self.ksstatistic = []
        self.pvalues = []
        self.DitherPatterns = {
            'random_night': 'randomDitherFieldPerNight',
            'random_visit': 'randomDitherFieldPerVisit',
            'spiral_night': 'spiralDitherFieldPerNight',
            'spiral_visit': 'spiralDitherFieldPerVisit',
            'hex_night': 'hexDitherFieldPerNight',
            'hex_visit': 'hexDitherFieldPerVisit',
            'field': 'field'
        }

        fieldidcolumn = 'fieldId'
        deg = True

        self.Stackers = {
         'random_night': stackers.RandomDitherFieldPerNightStacker(
                                        fieldIdCol=fieldidcolumn, degrees=deg),
         'random_visit': stackers.RandomDitherFieldPerVisitStacker(
                                        degrees=deg),
         'spiral_night': stackers.SpiralDitherFieldPerNightStacker(
                                        fieldIdCol=fieldidcolumn, degrees=deg),
         'spiral_visit': stackers.SpiralDitherFieldPerVisitStacker(
                                        fieldIdCol=fieldidcolumn, degrees=deg),
         'hex_night': stackers.HexDitherFieldPerNightStacker(
                                        fieldIdCol=fieldidcolumn, degrees=deg),
         'hex_visit': stackers.HexDitherFieldPerVisitStacker(
                                        fieldIdCol=fieldidcolumn, degrees=deg),
         'field': None
         }

        if DitherPattern not in self.DitherPatterns.keys():
            raise ValueError(
               'Dither Pattern not supported, please choose from: field, \
               random_night, random_visit, spiral_night, spiral_visit, \
               hex_night, hex_visit')
        else:
            self.DitherPattern = self.DitherPatterns[DitherPattern]
            self.Stacker = [self.Stackers[DitherPattern]]

        self.rotDitherPattern = rotDithers
        # if rotDithers is True:
        #     self.Stacker.append(
        #         stackers.RandomRotDitherPerFilterChangeStacker()
        #         )

        self.savedStarsAngles = {tuple(k): [] for k in self.stars}
        self.year = year

    class PSF:
        '''empty namespace to organise results into
        '''
        pass

    class STAR:
        '''empty namespace to organise results into
        '''
        pass

    class DELTA:
        '''empty namespace to organise results into
        '''
        pass

    def getPositions(self, sqlWhere):
        '''uses OpSim to find all dithered positions
        given some constraint sqlWhere

        Args:
            sqlWhere (str): sql query

        Returns:
            None: records results to self.positions
        '''
        print('getting the dither positions from the database: ',
              self.DitherPattern)
        print('using stackers:', self.Stacker)

        if self.Maker is not 'OpSim':

            opsdb = db.OpsimDatabase(self.runName+'.db')

            data = opsdb.fetchMetricData(['fieldRA', 'fieldDec',
                                          'filter', 'night'])

            posRA = data['fieldRA']
            posDec = data['fieldDec']
            posRA = np.array([pos[0]*np.radians(1) for pos in posRA])
            posDec = np.array([pos[0]*np.radians(1) for pos in posDec])
            filters = data['filter']
            nights = data['night']

            if self.Maker is 'Peter':
                notes = opsdb.fetchMetricData(['note'])
                ddf_cond = ['DD' not in (
                    str(note)
                    ).split(':')[0] for note in notes]
                cond = np.logical_and.reduce((ddf_cond,
                                              filters == 'i',
                                              nights <= int(self.year*365)))
            else:
                cond = np.logical_and(filters == 'i',
                                      nights <= int(self.year*365))
            posRA = posRA[cond]
            posDec = posDec[cond]

            pos = np.array((posRA, posDec))
            self.positions = pos.swapaxes(1, 0)

            return

        if self.Stacker[0] is None:
            self.Stacker = self.Stacker[1:]

        if self.Maker is 'OpSim':
            outDir = 'temp'
            myBundles = {}
            nside = 128
            resultsDb = db.ResultsDb(outDir=outDir)

            opsdb = db.OpsimDatabase(self.runName+'.db')
            slicer = slicers.HealpixSlicer(
                lonCol='fieldRA', latCol='fieldDec', nside=nside
                )

            metric = metrics.PassMetric()
            print('stacker is ', self.Stacker)
            if self.Stacker == [None]:
                myBundles['metric bundle'] = metricBundles.MetricBundle(
                    metric,
                    slicer,
                    constraint=sqlWhere,
                    runName=self.runName,
                    metadata='running metric')
            else:
                myBundles['metric bundle'] = metricBundles.MetricBundle(
                    metric,
                    slicer,
                    constraint=sqlWhere,
                    stackerList=self.Stacker,
                    runName=self.runName,
                    metadata='running metric')
            bgroup = metricBundles.MetricBundleGroup(myBundles, opsdb,
                                                     outDir=outDir,
                                                     resultsDb=resultsDb)
            bgroup.runAll()

            if self.DitherPattern == 'field':
                posRA = bgroup.simData['fieldRA']
                posDec = bgroup.simData['fieldDec']
            else:
                posRA = bgroup.simData[self.DitherPattern+'Ra']
                posDec = bgroup.simData[self.DitherPattern+'Dec']
            pos = np.array((posRA, posDec))

            pos *= np.radians(1)
            self.positions = pos.swapaxes(1, 0)

        # ROTATIONAL DITHERING
        if self.rotDitherPattern is True:
            print('getting random rotational dithers')
            data = opsdb.fetchMetricData(['filter'])
            filters = data['filter']
            groupedFilters = [
                    (
                      filt, len(list(occurence))
                    ) for filt, occurence in groupby(filters)
                ]

            d = []
            cutoff = 1440
            for j in groupedFilters:
                if j[1] <= cutoff:
                    d.append(j)
                else:
                    d = d + [(j[0], cutoff)]*(j[1]//cutoff)
                    d.append((j[0], j[1] % cutoff))

            self.rotTelPos = []
            for j in d:
                self.rotTelPos = self.rotTelPos + \
                    [np.random.uniform(-90, 90)]*j[1]
                self.rotTelPos = np.array(self.rotTelPos)
        else:
            print('not using rotational dithering')
            self.rotTelPos = None
        print('Summary: we are using the OpsimRun at {}, \
               with the translational dither pattern: {}. \
               Rotation: {}'.format(self.runName,
                                    self.DitherPattern,
                                    self.rotDitherPattern))

    def M2e(self):
        '''go back from 2nd moment space to elipticities
        '''
        print('moving back from moment space to elipticities')
        for pos in self.PSF.M.keys():
            Mxx, Mxy, Myy = self.DELTA.M[pos]
            self.DELTA.e[pos] = np.array([(Mxx-Myy)/self.PSF.TrM,
                                         Mxy*2/self.PSF.TrM])
            Mxx, Mxy, Myy = self.STAR.M[pos]
            self.STAR.e[pos] = np.array([(Mxx-Myy)/self.STAR.TrM,
                                        Mxy*2/self.STAR.TrM])
            Mxx, Mxy, Myy = self.PSF.M[pos]
            self.PSF.e[pos] = np.array([(Mxx-Myy)/self.PSF.TrM,
                                        Mxy*2/self.PSF.TrM])

    def process(self, sqlWhere):
        '''runs all analysis methods

        Args:
            sqlWhere (str): SQL query for the database.
        '''
        self.getPositions(sqlWhere)
        print('creating the models at every dither, this will take a while')
        for i in range(len(self.stars)):
            self.getModel(position_num=i)
        self.M2e()
        print('finding rhos/errors '+str(self.bootstrap_iterations)+' times')
        for bootstrap_iteration in range(self.bootstrap_iterations):
            self.getRhos()
            self.rhos2errors()

    def getRhos(self):
        '''method to get the rho statistics, needs a model,
        and traces for rhos 2 through 5.
        '''
        # print('finding rhos')
        seed = random.seed()
        itms = random.choices(
            list(self.DELTA.e.items()), seed, k=self.star_num
            )
        X, Y = [it[0][0] for it in itms], [it[0][1] for it in itms]
        de1, de2 = [it[1][0] for it in itms], [it[1][1] for it in itms]

        psfe1, psfe2 = np.array(
            [self.PSF.e[(x, y)] for x, y in zip(X, Y)]
            ).swapaxes(1, 0)
        stare1, stare2 = np.array(
            [self.STAR.e[(x, y)] for x, y in zip(X, Y)]
            ).swapaxes(1, 0)

        decat = tr.Catalog(g1=de1, g2=de2, ra=X, dec=Y,
                           ra_units='radians', dec_units='radians')
        psfcat = tr.Catalog(g1=psfe1, g2=psfe2, ra=X, dec=Y,
                            ra_units='radians', dec_units='radians')
        starcat = tr.Catalog(g1=stare1, g2=stare2, ra=X, dec=Y,
                             ra_units='radians', dec_units='radians')

        min_sep = 0.01  # in degrees
        max_sep = 10  # in degrees
        nbins = 24  # number of bins

        dedecorr = tr.GGCorrelation(min_sep=min_sep, max_sep=max_sep,
                                    nbins=nbins, sep_units='degrees')
        dedecorr.process(decat)
        dede_xip = dedecorr.xip
        dede_xim = dedecorr.xim
        self.r = np.exp(dedecorr.meanlogr)

        self.rho1 = dede_xip

        edecorr = tr.GGCorrelation(min_sep=min_sep, max_sep=max_sep,
                                   nbins=nbins, sep_units='degrees')
        edecorr.process(psfcat, decat)
        ede_xip = edecorr.xip
        ede_xim = edecorr.xim

        self.rho2 = ede_xip

        edtedtcorr = tr.GGCorrelation(min_sep=min_sep, max_sep=max_sep,
                                      nbins=nbins, sep_units='degrees')
        edtedtcorr.process(psfcat, psfcat)
        edtedt_xip = edtedtcorr.xip * self.size_error_ratio**2
        edtedt_xim = edtedtcorr.xim * self.size_error_ratio**2

        self.rho3 = edtedt_xip

        deedtcorr = tr.GGCorrelation(min_sep=min_sep, max_sep=max_sep,
                                     nbins=nbins, sep_units='degrees')
        deedtcorr.process(decat, psfcat)
        deedt_xip = deedtcorr.xip * self.size_error_ratio
        deedt_xim = deedtcorr.xim * self.size_error_ratio

        self.rho4 = deedt_xip

        eedtcorr = tr.GGCorrelation(min_sep=min_sep, max_sep=max_sep,
                                    nbins=nbins, sep_units='degrees')
        eedtcorr.process(psfcat, psfcat)
        eedt_xip = eedtcorr.xip * self.size_error_ratio
        eedt_xim = eedtcorr.xim * self.size_error_ratio

        self.rho5 = eedt_xip

        # star shape correlation:
        starcorr = tr.GGCorrelation(min_sep=min_sep, max_sep=max_sep,
                                    nbins=nbins, sep_units='degrees')
        starcorr.process(starcat)
        self.xip = starcorr.xip

    def rhos2errors(self):
        '''propagates rho statistics into shear errors
        '''

        delta_xip = 2 * self.size_error_ratio*self.trace_ratio * self.xip\
            + (self.trace_ratio)**2 * self.rho1\
            - self.alpha * (self.trace_ratio) * self.rho2\
            + (self.trace_ratio)**2 * self.rho3\
            + (self.trace_ratio)**2 * self.rho4\
            - self.alpha * (self.trace_ratio) * self.rho5

        self.xipList.append(delta_xip)

    def getModel(self, position_num):
        '''method to create a radial pattern (one of the simplified models),
        e = 0.05* distance from origin.
        New algorithm loops over each star, finding the dither positions
        at which it would be visible -- then

        Args:
            position_num (int): the dither position to consider

        Returns:
            None: records results in self.counter,
                  self.STAR.M, self.PSF.M, self.DELTA.M.
        '''
        star_pos = self.stars[position_num]
        cond = angularSeparation(
                            self.positions[:, 0]*np.degrees(1),
                            self.positions[:, 1]*np.degrees(1),
                            star_pos[0]*np.degrees(1),
                            star_pos[1]*np.degrees(1)
                            ) < self.FOVradius

        innerDithers = self.positions[cond]
        if len(innerDithers) < 2:
            return

        if self.rotDitherPattern is True:
            rotDithers = self.rotTelPos[cond]
        else:
            pass

        if self.ModelType == 'radial':
            stare1, stare2, psfe1, psfe2 = self.RadialModel(
                star_pos=star_pos,
                innerDithers=innerDithers
                )

        elif self.ModelType == 'horizontal':
            stare1, stare2, psfe1, psfe2 = self.HorizontalModel(
                star_pos=star_pos,
                innerDithers=innerDithers,
                rotDithers=rotDithers
                )

        star_pos = tuple(star_pos)
        starMxx = 0.5*self.STAR.TrM*(stare1+1)
        starMxy = 0.5*self.STAR.TrM*stare2
        starMyy = 0.5*self.STAR.TrM*(-stare1+1)
        psfMxx = 0.5*self.PSF.TrM*(psfe1+1)
        psfMxy = 0.5*self.PSF.TrM*psfe2
        psfMyy = 0.5*self.PSF.TrM*(-psfe1+1)
        angles = np.array(
            [r2d(0.5*np.arctan2(e2, e1)) for e2, e1 in zip(stare2, stare1)]
            )
        self.savedStarsAngles[star_pos] = angles

        deltaMxx = starMxx - psfMxx
        deltaMxy = starMxy - psfMxy
        deltaMyy = starMyy - psfMyy

        innerDithers = [tuple(i) for i in list(innerDithers)]
        self.counter[star_pos] = len(set(innerDithers))
        innerDithers = np.array(list(innerDithers))
        self.STAR.M[star_pos] = np.array([np.mean(starMxx),
                                          np.mean(starMxy),
                                          np.mean(starMyy)])
        self.PSF.M[star_pos] = np.array([np.mean(psfMxx),
                                         np.mean(psfMxy),
                                         np.mean(psfMyy)])
        self.DELTA.M[star_pos] = np.array([np.mean(deltaMxx),
                                           np.mean(deltaMxy),
                                           np.mean(deltaMyy)])

    def RadialModel(self, star_pos, innerDithers):
        """Summary

        Args:
            star_pos (tuple): (ra, dec) positions for stars to consider
            innerDithers (tuple): (ra, dec) positions for dithers
                                  that can see the star

        Returns:
            stare1: e1 ellipticity for truth
            stare2: e2 ellipticity for truth
            psfe1: e1 ellipticity for the PSF
            psfe2: e2 ellipticity for the PSF
        """
        r = angularSeparation(
            star_pos[0]*np.degrees(1), star_pos[1]*np.degrees(1),
            innerDithers[:, 0]*np.degrees(1), innerDithers[:, 1]*np.degrees(1))
        r *= np.radians(1)
        r[r < 0.0237] = 0  # 0.0237 is 80% of the LSST FOV radius in radians.
        r[r >= 0.0237] = 0.08
        rel_X = star_pos[0] - innerDithers[:, 0]
        rel_Y = star_pos[1] - innerDithers[:, 1]
        theta = np.arctan((rel_Y)/(rel_X))
        stare1 = r*np.cos(2*theta)
        stare2 = r*np.sin(2*theta)
        psfe1 = stare1/1.06
        psfe2 = stare2/1.06
        return stare1, stare2, psfe1, psfe2

    def HorizontalModel(self, star_pos, innerDithers, rotDithers):
        '''Summary
        this method creates a horizontal model model of e_psf = 0.06
        and the residual is 3%.

        Args:
            star_pos (tuple): (ra, dec) positions for stars to consider
            innerDithers (list of tuples): (ra, dec) positions for dithers
                                           that can see the star
            rotDithers (list of tuples): angles for dither positions
                                         that can see the star

        Returns:
            stare1: e1 ellipticity for truth
            stare2: e2 ellipticity for truth
            psfe1: e1 ellipticity for the PSF
            psfe2: e2 ellipticity for the PSF
        '''
        stare1 = np.cos(2*rotDithers)/5
        stare2 = np.sin(2*rotDithers)/5
        psfe1 = stare1/1.03
        psfe2 = stare2/1.03
        return stare1, stare2, psfe1, psfe2

    def getRequirements(self):
        '''Getting requirements on rhos and xi_+ from HSC data.

        Returns:
            reqs_r (ndarray): theta seperations where
                              requirements are computed.
            rho134_reqs (ndarray): requirements on rho1, 3 and 4.
            rho25_reqs (ndarray): requirements on rho 2 and 5.
        '''

        HSCCosmicShear = np.loadtxt(
            'HSCS16A_combinedarea_1000rea_full_xi_p.mean_sqrtvar')
        reqs_r = arcm2r(1)*HSCCosmicShear[:, 0]
        HSCCosmicShear_xip = HSCCosmicShear[:, 1]/2.6
        rho_reqs = np.loadtxt('rho_requirements.txt')
        lsst_area = 18000.  # in degrees
        hsc_area = 136.  # in degrees
        rho25_reqs = 0.02/self.alpha * np.sqrt(hsc_area/lsst_area) * \
            rho_reqs[:, 2]/(self.trace_ratio)**2
        rho134_reqs = np.sqrt(hsc_area/lsst_area) * \
            rho_reqs[:, 1]/(self.trace_ratio)
        return reqs_r, rho134_reqs


def getCounterAndDeltaXips(model,
                           year,
                           DitherPattern,
                           OpsimRun,
                           rotDithers,
                           objects_base='Y10'):
    """
    Args:
        model (str): 'radial' or 'horizontal'
        year (int): year of survey
        DitherPattern (str): 'field', 'spiral', 'hexagonal' or 'random'
        OpsimRun (str): name of the Opsim run's file without extension
        rotDithers (bool): Whether to use rotational dithering or not
        objects_base (str): year for making cut on objects 'Y10' or 'actual'
        proposal_format (str): 'default uses a proposalDict with values
                               representing'

    Returns:
        dict: the distribution of visits per object

    Raises:
        ValueError: if string args are not understood
    """
    proposalDict = {'baseline2018a': 3, 'colossus_2664': 2, 'colossus_2665': 1,
                    'colossus_2667': 1, 'kraken_2026': 3, 'kraken_2035': 3,
                    'kraken_2036': 3, 'kraken_2042': 2, 'kraken_2044': 1,
                    'mothra_2045': 1, 'nexus_2097': 1, 'pontus_2002': 1,
                    'pontus_2489': 3, 'pontus_2502': 2, 'pontus_2579': 3}
    countersDict = {}
    nightsNum = year*365
    if OpsimRun not in list(proposalDict.keys()):
        proposal_format = 'none'
    elif OpsimRun == 'pontus_2502':
        proposal_format = 'pontus_2502'
    else:
        proposal_format = 'default'
    if proposal_format == 'default':
        sqlWhere = 'night < '+str(nightsNum)+' and \
        filter = "i" and proposalId = '+str(proposalDict[OpsimRun])
    elif proposal_format == 'pontus_2502':
        sqlWhere = 'night < '+str(nightsNum)+' and \
        filter = "i" and proposalId != '+str(proposalDict[OpsimRun])
    elif proposal_format == 'none':
        sqlWhere = 'night < '+str(nightsNum)+' and \
        filter = "i"'
    else:
        raise ValueError('Cannot understand proposal_format')
    directory = 'newcutnpys/'
    outName = directory+OpsimRun+DitherPattern+str(rotDithers)+str(year)+'.npy'

    print('analysing'+OpsimRun)
    errors_object = ModelErrors(ModelType=model,
                                DitherPattern=DitherPattern,
                                OpsimRun=OpsimRun,
                                rotDithers=rotDithers,
                                year=year,
                                objects_base=objects_base)
    errors_object.process(sqlWhere)
    countersDict[OpsimRun] = errors_object.counter
    np.save(outName, errors_object.xipList)
    print('there are now {} runs total for this strategy'.format(
                                                len(np.load(outName))))
    return countersDict, errors_object


def getPositions(runName, year):
    proposalDict = {'baseline2018a': 3, 'colossus_2664': 2, 'colossus_2665': 1,
                    'colossus_2667': 1, 'kraken_2026': 3, 'kraken_2035': 3,
                    'kraken_2036': 3, 'kraken_2042': 2, 'kraken_2044': 1,
                    'mothra_2045': 1, 'nexus_2097': 1, 'pontus_2002': 1,
                    'pontus_2489': 3, 'pontus_2502': 2, 'pontus_2579': 3}
    directory = '/global/cscratch1/sd/husni/OpsimRuns/'
    nights = year*365 + 1
    opsdb = db.OpsimDatabase(directory+runName+'.db')
    outDir = 'temp'
    resultsDb = db.ResultsDb(outDir=outDir)
    nside = 256
    myBundles = {}
    if runName in list(proposalDict.keys()):
        sqlconstraint = 'filter = "i" and night < ' + str(nights) + \
            ' and proposalId = ' + str(proposalDict[runName])
    else:
        sqlconstraint = 'filter = "i" and night < '+str(nights)
    metric = metrics.ExgalM5(lsstFilter='i')
    dustMap = maps.DustMap(interp=False, nside=nside)
    stackerList = []
    slicer = slicers.HealpixSlicer(nside=nside)
    myBundles['field dither'] = metricBundles.MetricBundle(
        metric=metric,
        slicer=slicer,
        constraint=sqlconstraint,
        stackerList=stackerList,
        runName=runName,
        metadata='field dither',
        mapsList=[dustMap])
    bgroup = metricBundles.MetricBundleGroup(myBundles,
                                             opsdb,
                                             outDir=outDir,
                                             resultsDb=resultsDb)
    bgroup.runAll()
    bundle = myBundles['field dither']
    bundle.metricValues
    vminDict = {1: 24.5, 3: 25, 6: 25.5, 10: 26}
    vmin = vminDict[year]
    cond = np.logical_and.reduce((bundle.slicer.getSlicePoints()['ebv'] < 0.2,
                                  bundle.metricValues.mask is False,
                                  bundle.metricValues.data > vmin,
                                  bundle.metricValues.data < 28))
    condx = (bundle.slicer.getSlicePoints()['ra'])[cond]
    condy = (bundle.slicer.getSlicePoints()['dec'])[cond]
    conds = [np.array([cxi, cyi]) for cxi, cyi in zip(condx, condy)]
    a = np.array(random.sample(conds*100, 500000)) + \
        np.random.normal(0, 0.009, (500000, 2))
    np.save('newcutnpys/'+runName+str(year)+'.npy', a)
