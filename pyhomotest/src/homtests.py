import numpy as np
import lmoments3 as lm
from lmoments3 import distr
import math
import pandas as pd

'''
Eq. 4.3
'''
def getWeightedLMoments(regionalDurationStatisticsDF):
    l2r=0.
    tr=0.
    t3r=0.
    t4r=0.
    sumN = 0.
    sumNT = 0.
    
    for row in regionalDurationStatisticsDF.itertuples():
        l2r += row.nyears*row.l2/row.l1
        tr += row.nyears*row.t
        t3r += row.nyears*row.t3
        t4r += row.nyears*row.t4
        sumN += row.nyears

    l2r = l2r/sumN
    tr = tr/sumN
    t3r = t3r/sumN
    t4r = t4r/sumN
    return {'l2r':l2r, 'tr':tr, 't3r':t3r, 't4r':t4r}


def calculateWeightedLMoments(nYears=None, l2=None, t=None, t3=None, t4=None ):
    l2r = 0.
    tr = 0.
    t3r = 0.
    t4r = 0.
    sumN = 0.

    for i in len(nYears):
        nyrs = nYears[i]
        l2r += nyrs * l2[i]
        tr += nyrs * t[i]
        t3r += nyrs * t3[i]
        t4r += nyrs * t4[i]
        sumN += nyrs

    l2r = l2r / sumN
    tr = tr / sumN
    t3r = t3r / sumN
    t4r = t4r / sumN
    return {'l2r': l2r, 'tr': tr, 't3r': t3r, 't4r': t4r}

'''
Extracted from the distr.py - line 462 

Description
This function estimates the parameters of the Kappa distribution given the L-moments of the data in an L-moment 
object such as that returned by lmoms. The relations between distribution parameters and L-moments are seen under 
lmomkap, but of relevance to this documentation, the upper bounds of L-kurtosis (τ_4) and a function of L-skew (τ_3) is given by

τ_4 < \frac{5τ_3^2+1}{6}

This bounds is equal to the Generalized Logistic distribution (parglo) and failure occurs if this upper bounds is exceeded. 
However, the argument snap.tau4, if set, will set τ_4 equal to the upper bounds of τ_4 of the distribution to the relation above. 
This value of τ_4 should be close enough numerically The argument nudge.tau4 is provided to offset τ_4 downward just a little. 
This keeps the relation operator as “<” in the bounds above to match Hosking's tradition as his sources declare “≥” as above the GLO. 
The nudge here hence is not zero, which is a little different compared to the conceptually similar snapping in paraep4.
'''
def clean_lmom_ratio(lmom):
    lmom_ratios = [1, lmom['tr'], lmom['t3r'], lmom['t4r']]
    T3 = lmom_ratios[2]
    T4 = lmom_ratios[3]
    
    if lmom_ratios[1] <= 0: 
        raise ValueError("lmom_ratios[1] <= 0")
    elif abs(T3) >= 1:
        raise ValueError('abs(T3) >= 1')
    elif abs(T4) >= 1:
        raise ValueError('abs(T4) >= 1')
    elif T4 <= (5 * T3 * T3 - 1) / 4:
        raise ValueError('T4 <= (5 * T3 * T3 - 1) / 4')
    elif T4 >= (5 * T3 * T3 + 1) / 6:
        new_t4 = (5 * T3 * T3 + 1) / 6. - 0.000001
        lmom['t4r'] = new_t4
        # TAU3 and TAU4 are above Generalized Logistic line.
        # https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/fExtremes/demo/funLmomco.R?revision=27&root=rmetrics&pathrev=1281
        # https://rdrr.io/cran/lmomco/man/parkap.html
        #raise ValueError('T4 >= (5 * T3 * T3 + 1) / 6 - t3={0} t4={1}'.format(T3, T4))
        
    return lmom


# https://rdrr.io/cran/homtest/src/R/HOMTESTS.R
class HomogeneityTest(object):
    def __init__(self, regionalDurationStatisticsDF, rlmoments):
        self.__v_array = []
        self.__regionalDurationStatisticsDF = regionalDurationStatisticsDF
        self.__rlmoments = rlmoments
        self.__vObs = self.get_V1(regionalDurationStatisticsDF,rlmoments)
        
        
        kapDistParams = distr.kap.lmom_fit(lmom_ratios=[1, rlmoments['tr'], rlmoments['t3r'], rlmoments['t4r']])
        #kapDistParams = distr.kap.lmom_fit(lmom_ratios=rlmoments)
        #print(kapDistParams)
        self.__kapDist = distr.kap(**kapDistParams)
    
    
    
    def get_V1(self, regionalDurationStatisticsDF, rlmoments):
        v=0.
        sumN = 0.
        for row in regionalDurationStatisticsDF.itertuples():
            sumN += row.nyears
            v += row.nyears * (row.t - rlmoments['tr'])**2
        v = math.sqrt(v/sumN)
        return v
        
    def __run_a_simulation(self):
        simulatedDurationStatisticsDF = pd.DataFrame()
        for row in self.__regionalDurationStatisticsDF.itertuples():
            years = row.nyears
            # Generate synethic AMS
            ams = self.__kapDist.rvs(years)
            ratios = lm.lmom_ratios(ams, nmom=5)
            data = {'nyears':years,
                    'l1':ratios[0], 
                    'l2':ratios[1], 
                    't':ratios[1]/ratios[0],
                    't3':ratios[2], 
                    't4':ratios[3], 
                    't5':ratios[4]}
            simulatedDurationStatisticsDF = simulatedDurationStatisticsDF.append(data, ignore_index=True)
        # Calculate the regionalized moments
        simulatedRlmoments = getWeightedLMoments(simulatedDurationStatisticsDF)

        v = self.get_V1(simulatedDurationStatisticsDF, simulatedRlmoments)
        self.__v_array.append(v)
    
    def calculate_heterogeneity(self, Nsims):
        # Run NSims
        for i in range(Nsims):
            self.__run_a_simulation()
        vSims = np.array(self.__v_array)
        muV =  vSims.mean()
        sigmaV =  vSims.std()
        
        H = (self.__vObs - muV)/sigmaV
        
        return H
        
        
