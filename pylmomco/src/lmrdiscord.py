
from scipy.stats import f
import numpy as np
from numpy.linalg import inv


# https://rdrr.io/cran/lmomco/src/R/lmrdiscord.R
# https://rdrr.io/cran/lmomco/man/lmrdiscord.html
def lmrdiscord(sites=None, t2=None, t3=None, t4=None, alpha1=0.10, alpha2=0.01):

    if alpha1 < alpha2:
        alpha = alpha1
        alpha1 = alpha2
        alpha2 = alpha
        print('Warning - inverted significance levels, flipping for you')
    
    nsites =len(sites)
    
    if nsites < 5:
        raise ValueError("Too few sites (<5), although >10 is preferred")
    
    # Hosking and Wallis (1997, table 3.1)
    Dcr_table_31 = np.array([1.3333, 1.6481, 1.9166, 2.1401, 2.3287,
                       2.4906, 2.6321, 2.7573, 2.8694, 2.9709, 3.])
    
    # Setup the matrix of t, t3 and t4
    lmrs = np.column_stack((t2, t3, t4))

    lmrmeans = lmrs.mean(0)
    tmp = np.subtract(lmrs,lmrmeans)

    A = np.zeros((3,3))
    for j in range(0,3):
        for k in range(0,3):
            for site in range(0,nsites):
                A[j,k] = A[j,k] + tmp[site,j] * tmp[site,k]

    Ainv = inv(A)

    # Hosking and Wallis (1997, eq. 3.3)
    D = np.zeros(nsites)
    for site in range(0,nsites):
        for j in range(0,3):
            for k in range(0,3):
                D[site] = D[site] + tmp[site,j] * tmp[site,k] * Ainv[j,k]
    
    # Finally calculate the discordancy measure for each station.
    D = D*nsites/3.

    # Determine whether D value exceeds critical value
    Dcrit = Dcr_table_31.max()

    if nsites < 15:
        Dcrit = Dcr_table_31[nsites-5]

    isD = np.full(nsites, False, dtype=bool)
    for site in range(nsites):
        if D[site] > Dcrit:
            isD[site] = True

    # Statistical tests to determine significance
    V = alpha1/nsites
    Z = f.isf(V, 3, nsites-4.)
    D1_by_Fdist = (nsites - 1)*Z/(nsites - 4 + 3*Z)


    V = alpha2/nsites
    Z = f.isf(V, 3, nsites-4.)
    D2_by_Fdist = (nsites - 1)*Z/(nsites - 4 + 3*Z)


    signif = []
    for site in range(nsites):
        if D[site] > D2_by_Fdist:
            signif.append("**")
        elif D[site] > D1_by_Fdist:
            signif.append("*")
        else:
            signif.append("-")
    
    Dmax = np.ones(nsites) * (nsites - 1.)/3.
    Dalpha1 = np.ones(nsites) * D1_by_Fdist
    Dalpha2 = np.ones(nsites) * D2_by_Fdist
    Dcrit = np.ones(nsites) * Dcrit
    
    retvalDf = pd.DataFrame({'site':sites, 't2':t2, 't3':t3, 't4':t4,'D':D, 'Dmax':Dmax, 'Dalpha1':Dalpha1, 'Dalpha2':Dalpha2, 'Dcrit':Dcrit,'signif':signif})
    
    return retvalDf