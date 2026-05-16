# ocean_functions.py
# Version 1.0.0
# Last modified: May 14th, 2026
# Authors: Matthew J. Graham, Adrien Hélias

import irlb
import numpy as np
from scipy.stats import skew, kurtosis
from spectrum import dpss
import sys
import warnings
warnings.simplefilter("ignore", RuntimeWarning) 


def power_iteration(input_matrix: np.ndarray,
                    vector: np.ndarray,
                    error_tol: float = 1e-16,
                    max_iterations: int = 1000
                    ) -> tuple[float, np.ndarray]:
    """
    Finds the largest eigenvalue and corresponding eigenvector
    of matrix input_matrix given a random vector in the same space.
    Will work so long as vector has component of largest eigenvector.
    input_matrix must be either real or Hermitian.

    Inputs
    ----------
    input_matrix : 2D Numpy array
        M x M input matrix whose largest eigenvalue we will find.
    vector : 1D Numpy array
        Random initial vector of size M in same space as the matrix.
    error_tol : float, optional
        Error tolerance. The default is 1e-16.
    max_iterations : int, optional
        Maximum number of iterations. The default is 1000.

    Outputs
    ----------
    largest_eigenvalue : float
        Largest eigenvalue of the matrix input_matrix.
    largest_eigenvector : 1D Numpy array
        Eigenvector of size M corresponding to the largest eigenvalue.
    """
    
    # Ensure matrix is square.
    assert np.shape(input_matrix)[0] == np.shape(input_matrix)[1]
    # Ensure proper dimensionality.
    assert np.shape(input_matrix)[0] == np.shape(vector)[0]
    # Ensure inputs are either both complex or both real.
    assert np.iscomplexobj(input_matrix) == np.iscomplexobj(vector)
    is_complex = np.iscomplexobj(input_matrix)
    if is_complex:
        # Ensure complex input_matrix is Hermitian.
        assert np.array_equal(input_matrix, input_matrix.conj().T)

    # Set convergence to False. Will define convergence when we exceed max_iterations
    # or when we have small changes from one iteration to next.

    convergence = False
    lambda_previous = 0
    iterations = 0
    error = 1e12

    while not convergence:
        # Multiple matrix by the vector.
        w = np.dot(input_matrix, vector)
        # Normalize the resulting output vector.
        vector = w / np.linalg.norm(w)
        # Find Rayleigh quotient.
        # (faster than usual b/c we know vector is normalized already)
        vector_h = vector.conj().T if is_complex else vector.T
        lambda_ = np.dot(vector_h, np.dot(input_matrix, vector))

        # Check convergence.
        error = np.abs(lambda_ - lambda_previous) / lambda_
        iterations += 1

        if error <= error_tol or iterations >= max_iterations:
            convergence = True

        lambda_previous = lambda_

    if is_complex:
        lambda_ = np.real(lambda_)

    return lambda_, vector


def qmatrix(c: float,
            jscale: int,
            mu: float,
            mesh: np.ndarray
            ) -> np.ndarray:
    """
    Constructs the concentration matrix Q used in Slepian wavelet calculations.
    This function builds a sinc-based matrix that represents the
    prolate spheroidal wave equation in discretized form.
    
    Inputs
    ----------
    c : float
        c is a constant such that 2c is an integer. Usually set to 1.
    jscale : int
        The scale index (octave level) for the wavelet. Used in exponent
        operations like 2**(jscale + 1) and 2**jscale to set frequency bounds.
    mu : float
        The mean sampling interval.
    mesh : 1D Numpy array
        The time points over which the wavelet is defined. Should have exactly
        M = c*2**jscale elements. In practice, this is a subset of the
        observation times extracted from the full time series.
    
    Output
    ----------
    slep : 2D Numpy array
        M x M matrix where M = c*2**jscale. A symmetric, real-valued matrix
        representing the concentration matrix for the prolate spheroidal
        wave equation.
        This matrix is constructed to have maximum energy concentration in
        the frequency band [fl, fr], and its largest eigenvector becomes the
        Slepian wavelet used for analyzing the time series.
    """
    
    fl = 1. / 2**(jscale + 1) / mu
    fr = 1. / 2**jscale / mu
    mm = (2 ** jscale) * c
    mat = np.repeat(mesh, mm).reshape(mm, mm) - np.repeat(mesh, mm).reshape(mm, mm).T
    slep = (np.sin(2. * np.pi * mat * fr) - np.sin(2. * np.pi * mat * fl)) / (np.pi * mat)
    np.fill_diagonal(slep, 2. * (fr - fl))
    iden = np.identity(mm) - 1./mm
    slep = np.dot(iden, np.dot(slep, iden))
    return slep


def eigenvector(c: float,
                jscale: int,
                mu: float,
                mesh: np.ndarray
                ) -> np.ndarray:
    """
    Extracts the largest eigenvector from the concentration matrix Q using
    power iteration. This computes which Slepian wavelet direction has the
    most energy.
    The function calls power_iteration() with a random starting vector to
    find the largest eigenvector of the concentration matrix produced by
    qmatrix().
    
    Inputs
    ----------
    c : float
        c is a constant such that 2c is an integer. Usually set to 1.
    jscale : int
        The scale index (octave level) for the wavelet. Used in exponent
        operations like 2**(jscale + 1) and 2**jscale to set frequency bounds.
    mu : float
        The mean sampling interval.
    mesh : 1D Numpy array
        The time points over which the wavelet is defined. Should have exactly
        M = c*2**jscale elements. In practice, this is a subset of the
        observation times extracted from the full time series.

    Output
    ----------
    vec : 1D Numpy array
        The largest eigenvector of the concentration matrix computed via
        power iteration. The vector is of size M.
    """
    slep = qmatrix(c, jscale, mu, mesh)
    test = np.random.rand(slep.shape[0])
    power, vec = power_iteration(slep, test)
    return vec
  

def eigenvector2(c: float,
                 jscale: int,
                 mu: float,
                 mesh: np.ndarray
                 ) -> np.ndarray:
    """
    Computes the largest eigenvector for large matrices using an optimized
    algorithm (IRLB). This is an alternative to eigenvector() for efficiency.
    Instead of power iteration, it uses the irlb (Implicitly Restarted Lanczos
    Bidiagonalization) package, which is faster for larger matrices.
    It extracts the first singular vector from the result.

    Inputs
    ----------
    c : float
        c is a constant such that 2c is an integer. Usually set to 1.
    jscale : int
        The scale index (octave level) for the wavelet. Used in exponent
        operations like 2**(jscale + 1) and 2**jscale to set frequency bounds.
    mu : float
        The mean sampling interval.
    mesh : 1D Numpy array
        The time points over which the wavelet is defined. Should have exactly
        M = c*2**jscale elements. In practice, this is a subset of the
        observation times extracted from the full time series.

    Output
    ----------
    1D Numpy array
        The largest transposed singular vector from the IRLB algorithm.
        The vector is of size M.
    """
    # Use irlb to calculate largest eigenvector for large matrix
    slep = qmatrix(c, jscale, mu, mesh)
    S = irlb.irlb(slep, 1, tol = 1.e-5)
    return S[0].T


def slepian(c: float,
            jscale: int,
            mu: float,
            mesh: np.ndarray
            ) -> np.ndarray:
    """
    Calculates a normalized Slepian wavelet by retrieving its eigenvector
    and scaling it.
    We normalize it by dividing by sqrt(mu*2**jscale) to produce the final 
    wavelet.
    
    Inputs
    ----------
    c : float
        c is a constant such that 2c is an integer. Usually set to 1.
    jscale : int
        The scale index (octave level) for the wavelet. Used in exponent
        operations like 2**(jscale + 1) and 2**jscale to set frequency bounds.
    mu : float
        The mean sampling interval.
    mesh : 1D Numpy array
        The time points over which the wavelet is defined. Should have exactly
        M = c*2**jscale elements. In practice, this is a subset of the
        observation times extracted from the full time series.

    Output
    ----------
    evs : 1D Numpy array
        The normalized, largest eigenvector of the concentration matrix computed via
        power iteration.
    """
    # Calculate Slepian wavelets directly
    evs = eigenvector(c, jscale, mu, mesh) 
    return evs / np.sqrt(mu*2.**jscale)


def slepian2(c: float,
             jscale: int,
             mu: float,
             mesh: np.ndarray
             ) -> np.ndarray:
    """
    Alternative Slepian wavelet calculation using the IRLB method for
    potentially faster computation on large matrices.
    Identical to slepian(), but calls eigenvector2() instead, making it more
    efficient for larger datasets.
    
    Inputs
    ----------
    c : float
        c is a constant such that 2c is an integer. Usually set to 1.
    jscale : int
        The scale index (octave level) for the wavelet. Used in exponent
        operations like 2**(jscale + 1) and 2**jscale to set frequency bounds.
    mu : float
        The mean sampling interval.
    mesh : 1D Numpy array
        The time points over which the wavelet is defined. Should have exactly
        M = c*2**jscale elements. In practice, this is a subset of the
        observation times extracted from the full time series.

    Output
    ----------
    evs : 1D Numpy array
        The normalized, largest transposed singular vector from the IRLB algorithm.
    """
    # Calculate Slepian wavelets with the alternate IRLB method
    evs = eigenvector2(c, jscale, mu, mesh) 
    return evs[0] / np.sqrt(2. ** jscale * mu)


def slepwav(times: np.ndarray,
            obsx: np.ndarray,
            obsxerr: np.ndarray
            ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Performs the complete Slepian Wavelet Variance analysis on time series
    data. This is the main wrapper function.
    Applied to quasar light curves, this produces the characteristic wavelet
    variance scaling relationship used to characterize quasar variability
    across different timescales.

    Inputs
    ----------
    times : 1D Numpy array
        The time (X) values of the time series.
    obsx : 1D Numpy array
        The flux/magnitude (Y) values of the time series.
    obsxerr : 1D Numpy array
        The errors on the flux/magnitude (Y) values.

    Outputs
    -------
    t2 : 1D Numpy array
        The log2 values of the timescales of the wavelet coefficients.
    v2 : 1D Numpy array
        The log2 values of the variance of the wavelet coefficients.
    verr2 : 1D Numpy array
        The error on the log2 values of the variance of the wavelet coefficients.
    skw : 1D Numpy array
        The skewness of the wavelet coefficients.
    kur : 1D Numpy array
        The kurtosis of the wavelet coefficients.
    """
    
    n = len(obsx) - 2
    times = times - times[0]
    meandelta = times[n + 1] / (n + 1)  # Mean dt
    cons = 1
    maxscale = int(np.floor(np.log2(n + 2))) + 3
    sLj = [2**(j + 1) * cons for j in range(maxscale)]
    sMj = [n - sLj[j] + 1 for j in range(maxscale)]
    mu = meandelta
    dps = []
    lamplus = []
    interplim = 11  # Different eigensolver used for jscale + 1 > interplim
    
    # Calculate Slepian wavelets
    sw_var = np.zeros(maxscale)
    sw_skew = np.zeros(maxscale)
    sw_kur = np.zeros(maxscale)
    tau = np.zeros_like(sw_var)
    sw_dis = np.zeros_like(sw_var)
    wc = np.empty((maxscale, n))
    for jscale in range(maxscale):
        if sMj[jscale] / 2. <= 3.5: continue # NW < N / 2
        a, b = dpss(sMj[jscale], 3.5, 5) # N, NW, k
        dps.append(a)
        lamplus.append(np.repeat(1, sMj[jscale]).dot(dps[jscale]))
        for tpt in range(sMj[jscale]):
            index = range(tpt + 1, tpt + sLj[jscale] + 1)
            if jscale + 1 > interplim:
                wf = slepian2(cons, jscale + 1, mu, times[index])
            else:  
                wf = slepian(cons, jscale + 1, mu, times[index])
            try:
                wc[jscale, tpt] = np.dot(obsx[index], wf)
            except ValueError:
                print(jscale, tpt, index, obsx[index], wf)
                sys.exit(-1)
        ser = wc[jscale, 0:sMj[jscale]] ** 2.
        sw_var[jscale] = np.mean(ser)
        sw_skew[jscale] = skew(wc[jscale, 0:sMj[jscale]])
        sw_kur[jscale] = kurtosis(wc[jscale, 0:sMj[jscale]])
        tau[jscale] = 2. ** (jscale) * meandelta
        Jser = ser.dot(dps[jscale])
        u0 = Jser.dot(lamplus[jscale].conj().T) / np.sum(lamplus[jscale] ** 2.)
        sw_dis[jscale] = np.mean((Jser - np.array(u0) * lamplus[jscale]) ** 2.) / sMj[jscale]
    idx = np.where(((sw_var != np.inf) & (sw_var != np.nan) & (sw_var > 0)))[0]
    t2 = np.log2(tau[idx])
    v2 = np.log2(sw_var[idx])
    verr2 = np.sqrt(sw_dis[idx]) / (np.log(2) * sw_var[idx])
    skw = sw_skew[idx]
    kur = sw_kur[idx]
    return t2, v2, verr2, skw, kur


def run_slepian_wavelet_variance(X: np.ndarray,
                                 Y: np.ndarray,
                                 Yerr: np.ndarray,
                                 three_sigma_filter: bool = False,
                                 bin_light_curve: bool = False,
                                 bin_value: float = 1
                                 ) -> np.ndarray:
    """
    Runs Slepian Wavelet Variance (SWV) analysis on a time series.

    Inputs
    ----------
    X : np.ndarray
        The time values of the time series.
    Y : np.ndarray
        The flux/magnitude values of the time series.
    Yerr : np.ndarray
        The errors on the flux/magnitude values.
    three_sigma_filter : bool, optional
        An option to apply a 3-sigma filter on the time series before SWV
        processing. The default is False.
    bin_light_curve : bool, optional
        An option to bin the time series by a given time interval bin_value,
        before SWV processing. The default is False.
    bin_value : float, optional
        The time interval to bin the time series with, in the units of X. The default is 1.

    Output
    -------
    ans : 1D Numpy array
        An array containing all the parameters obtained from the SWV analysis.
        See the documentation of the function slepwav for more details.
    """
    
    print("Running Slepian Wavelet Variance analysis...")
    if three_sigma_filter == True:
        cleaned_Y = abs(Y - np.mean(Y)) < 3*np.std(Y)
        Y = cleaned_Y
    if bin_light_curve == True:
        # Binning ZTF light curve, weighting by inverse variance
        t_min = X.min()
        t_max = X.max()
        n_bins = int(np.ceil((t_max - t_min)/bin_value))
        edges = t_min + np.arange(n_bins + 1)*bin_value
        inds = np.digitize(X, edges) - 1
        
        # Prepare outputs
        bin_center = (edges[:-1] + edges[1:])/2.0
        weighted_mean = np.full(n_bins, np.nan)
        weighted_err = np.full(n_bins, np.nan)
        count = np.zeros(n_bins, dtype=int)
        w = 1/(Yerr**2)
        
        # Aggregate per bin
        for b in range(n_bins):
            sel = inds == b
            if not np.any(sel):
                continue
            sw = w[sel].sum()
            if sw == 0:
                continue
            weighted_mean[b] = (w[sel]*Y[sel]).sum()/sw
            weighted_err[b] = np.sqrt(1/sw)
            count[b] = sel.sum()
            
        X_binned = bin_center
        Y_binned = weighted_mean
        Yerr_binned = weighted_err
        mask = ~np.isnan(Y_binned)  # True where Y is not NaN
        X_binned_nonan = X_binned[mask]
        Y_binned_nonan = Y_binned[mask]
        Yerr_binned_nonan = Yerr_binned[mask]
        X = X_binned_nonan
        Y = Y_binned_nonan
        Yerr = Yerr_binned_nonan
    
    t, v, verr, sk, k = slepwav(X, Y, Yerr)
    ans = np.array([t, v, verr, sk, k])
    # print(ans.T)
    # dt = np.diff(X)
    # print("Mean dt: " + str(np.mean(dt)) + " days")
    print("Done!")
    return ans


