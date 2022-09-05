"""
Calculate Signal-to-Clutter Ratio (SCR) from a set of interferograms
"""
from psps import sario

import glob
import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from scipy.ndimage import uniform_filter
from scipy.special import erfc
PSDIR = 'ps_data'

def amp_dispersion(amplist,nrow,ncol):
    """
    Calculate amplitude dispersion

    Parameters
    ----------
    amplist : list of string
        list of amplitude filenames, assuming the amplitude files are stored
        in float32 format
    nrow : int
        number of rows of each amplitude image
    ncol : int
        number of columns of each amplitude image

    Returns
    -------
    amp_disp : 2D array (nrow x ncol)
        amplitude dispersion estimated from the amplitude image list
    """
    n = len(amplist)
    npts = nrow*ncol
    amp_sum = np.zeros(nrow*ncol)
    amp2_sum = np.zeros(nrow*ncol)
    for ampfile in amplist:
        amp = np.fromfile(ampfile,dtype='float32')
        amp = amp/np.mean(amp)
        amp_sum = amp_sum + amp
        amp2_sum = amp2_sum + amp**2
    amp_std = np.sqrt((amp2_sum - amp_sum**2/n)/n)
    amp_mean = amp_sum/n
    amp_disp = amp_std/(sys.float_info.epsilon+amp_mean)
    amp_disp = np.reshape(amp_disp,(nrow,ncol))
    amp_disp = amp_disp.astype('float32')
    return amp_disp

def conv_kernel(n,m,L):
    '''
    Generate a 2d box-car filter kernel associated with a 2D image of size
    n x m

    Parameters
    ----------
    n : int
        number of rows of the image
    m : int
        number of column of the image
    L : int
        window size (square box)

    Returns
    -------
    k : 2D array (nfft x mfft)
        2D box-car filter kernel: nfft is the smallest power of 2 greater or
        equal to n and mfft is the smallest power of 2 greater or equal to m
    '''
    nfft = int(2**(np.ceil(np.log2(n+L-1))))
    mfft = int(2**(np.ceil(np.log2(m+L-1))))
    k = np.zeros((nfft,mfft))
    halfL = int(L/2)
    if halfL > 0:
        k[0:halfL,0:halfL] = 1./L**2
        k[-halfL:,0:halfL] = 1./L**2
        k[0:halfL,-halfL:] = 1./L**2
        k[-halfL:,-halfL:] = 1./L**2
    else:
        k[0,0] = 1.
    return k

def filter2(ifg,windowsize,k_fft=None):
    """
    Apply a 2D box-car filter on an interferogram

    Parameters
    ----------
    ifg : 2D complex array
        an interferogram
    windowsize : int
        size of the box-car filter
    k_fft : 2D array (optional)
        kernel of the box-car filter. If k_fft is None, `conv_kernel` will be
        called to generate a kernel

    Returns
    -------
    ifg_filtered : 2D complex array
        filtered interferogram
    """
    if k_fft is None:
        k = conv_kernel(ifg.shape[0],ifg.shape[1],windowsize)
        k_fft = np.fft.fft2(k)
    temp = np.zeros(k_fft.shape,dtype=np.complex64)
    temp[0:ifg.shape[0],0:ifg.shape[1]] = ifg
    ifg_fft = np.fft.fft2(temp)
    ifg_filtered = np.fft.ifft2(ifg_fft*np.conj(k_fft))
    ifg_filtered = ifg_filtered[0:ifg.shape[0],0:ifg.shape[1]]
    return ifg_filtered
    
def phase_residue(ifglist,nrow,ncol,windowsize=11,output_dir=PSDIR):
    """
    Calculate phase residue of interferograms in `ifglist`. Each interferogram
    is filtered using a box-car filter. The moving-average phase is subtracted
    from the original phase to get the phase residue, which can be used for
    SCR calculation.

    Parameters
    ----------
    ifglist : list of string
        list of interferogram files
    nrow : int
        number of rows of each amplitude image
    ncol : int
        number of columns of each amplitude image
    windowsize : int (optional)
        box-car filter size. windowsize is 11 by default
    output_dir : string
        directory to save phase residue files
    """
    if not os.path.exists(output_dir):
        try:
            os.system('mkdir '+output_dir)
        except Exception as e:
            print(e)
            exit()
    k = conv_kernel(nrow,ncol,windowsize)
    k_fft = np.fft.fft2(k)
    for ifgfile in ifglist:
        print(ifgfile)
        respfile = os.path.split(ifgfile)[-1]+'.resphase'
        ifg = sario.read_ifg(ifgfile,nrow,ncol)
        ifg_filtered = filter2(ifg,windowsize,k_fft)
        ifg_residue = ifg*np.conj(ifg_filtered)/(1e-5+np.abs(ifg_filtered))
        sario.save_ifg(ifg_residue,os.path.join(output_dir,respfile))

def pdf_gaussian_signal(gamma,phi):
    """
    pdf of interferometric phase assuming both signal and noise follow complex
    Gaussian distribution

    Parameters
    ----------
    gamma : float
        signal-to-clutter ratio (SCR)
    phi : float or float array
        interferometric phase

    Returns
    -------
    Probability density estimated at phi given SCR = gamma
    """
    rho = gamma/(1+gamma)
    beta = rho*np.cos(phi)
    f = (1-rho*rho)/(2*np.pi*(1-beta**2))* \
            (1+beta*np.arccos(-beta)/np.sqrt(1-beta**2))
    return f

def phase_distr_const(gamma,theta):
    """
    pdf of SAR phase assuming a constant signal and a noise term following
    the complex Gaussian distribution

    Parameters
    ----------
    gamma : float
        signal-to-clutter ratio (SCR)
    theta : float or float array
        SAR phase

    Returns
    -------
    Probability density estimated at theta given SCR = gamma
    """
    p = 1/(2*np.pi)*np.exp(-np.sin(theta)**2*gamma)* \
            (np.exp(-np.cos(theta)**2*gamma)+ \
             np.sqrt(np.pi*gamma)*np.cos(theta)*erfc(-np.sqrt(gamma)*np.cos(theta)))
    return p

def joint_phase_distr_const(gamma,phi,phi_sum):
    p = 0.5*(phase_distr_const(gamma,(phi+phi_sum)/2)* \
            phase_distr_const(gamma,(phi_sum-phi)/2)+ \
            phase_distr_const(gamma,(phi+phi_sum)/2+np.pi)* \
            phase_distr_const(gamma,(phi_sum-phi)/2+np.pi))
    return p

def int_phase_distr_const(gamma,phi,n):
    f = 0.
    for i in range(n):
        phi_sum = 2*i*np.pi/n-np.pi
        f = f + 2*np.pi/n*joint_phase_distr_const(gamma,phi,phi_sum)
    return f

def pdf_constant_signal(gamma,phi):
    """
    pdf of interferometric phase assuming a constant signal and a noise term
    following the complex Gaussian distribution

    Parameters
    ----------
    gamma : float
        signal-to-clutter ratio (SCR)
    phi : float or float array
        interferometric phase

    Returns
    -------
    Probability density estimated at phi given SCR = gamma
    """
    n = 100
    phi_test = np.arange(n+1)*2*np.pi/n-np.pi
    lookup_tbl = int_phase_distr_const(gamma,phi_test,n)
    idx = np.round((phi+np.pi)/(2*np.pi/n)).astype(int)
    idx = np.minimum(idx,len(lookup_tbl))
    idx = np.maximum(idx,0)
    f = lookup_tbl[idx] 
    return f

def _cal_scr(phi,model='constant'):
    rho = np.linspace(0.0,0.99,49)
    scrs = rho/(1-rho)
    if len(phi.shape) == 1:
        phi = np.expand_dims(phi,axis=1)
    loglikelihood = np.zeros((len(scrs),phi.shape[1]))
    for i in range(len(scrs)):
        if model == 'gaussian':
            p = pdf_gaussian_signal(scrs[i],phi)
        elif model == 'constant':
            p = pdf_constant_signal(scrs[i],phi)
        # SCR calculation will be the same as Piyush's mlestack code if the
        # phase distribution model is 'sarpixel'
        elif model == 'sarpixel':
            p = phase_distr_const(scrs[i],phi)
        loglikelihood[i,:] = np.sum(np.log(p),axis=0)
    optscr = scrs[np.argmax(loglikelihood,axis=0)]
    return optscr

def cal_scr(ifglist, nrow, ncol, ntilerow=1, ntilecol=1, windowsize = 11,
            candidate_mask=None, model='constant',output_dir=PSDIR):
    """
    Calculate signal to clutter ratio based on a list of interferograms over the
    same area. In case of large interferograms, set ntilerow and ntilecol
    greater than 1. The interferogram will be partitioned into ntilerow x
    ntilecol tiles and processed them separately.

    Parameters
    ----------
    ifglist: list of str
        list of interferograms
    nrow: int
        number of rows of each interferogram
    ncol: int
        number of columns of each interferogram
    ntilerow: int (optional, default:1)
        number of rows of tiles
    ntilecol: int (optional, default:1)
        number of columns of tiles 
    windowsize: int (optional, default: 11)
        size of the box-car window to estimate the local average phase
    candidate_mask: 2d boolean array (optional)
        candidate_mask[i,j] = True if the pixel (i,j) needs to be processed
    model: str (optional)
        model of interferometric phase used for SCR calculation ('gaussian' or
        'constant')
    output_dir: string (optional)
        path of phase residue files

    Returns
    -------
    scr : 2D array (nrow x ncol)
        singal-to-clutter ratio (SCR)
    """
    phase_residue(ifglist,nrow,ncol,windowsize=windowsize,output_dir=output_dir)
    resphase_list = [os.path.join(output_dir,os.path.split(s)[-1]+'.resphase') \
                     for s in ifglist]
    n = len(ifglist)
    n1 = int(np.ceil(nrow/ntilerow))
    n2 = int(np.ceil(ncol/ntilecol))
    idx_list1 = np.arange(ntilerow+1)*n1
    idx_list1[-1] = nrow
    idx_list2 = np.arange(ntilecol+1)*n2
    idx_list2[-1] = ncol
    scr = np.zeros((nrow,ncol),dtype=np.float32)
    if candidate_mask is None:
        candidate_mask = np.ones((nrow,ncol),dtype=bool)
    for i in range(ntilerow):
        for j in range(ntilecol):
            print('processing row {}/{}, column {}/{}'.format(i+1,ntilerow,j+1,ntilecol))
            idx11 = idx_list1[i]
            idx12 = idx_list1[i+1]
            idx21 = idx_list2[j]
            idx22 = idx_list2[j+1]
            phi = sario.read_tile(resphase_list,nrow,ncol,idx11,idx12,idx21,idx22)
            cmean = np.mean(phi,axis=0)
            phi = np.angle(phi*np.conj(cmean))

            mask_tile = candidate_mask[idx11:idx12,idx21:idx22].flatten()
            phi = np.reshape(phi,(phi.shape[0],phi.shape[1]*phi.shape[2]))
            phi = phi[:,np.where(mask_tile)[0]]
            scr_tile = _cal_scr(phi,model) 
            tmp = np.zeros((idx12-idx11)*(idx22-idx21))
            tmp[mask_tile] = scr_tile
            scr[idx11:idx12,idx21:idx22] = np.reshape(tmp,(idx12-idx11,idx22-idx21))
    os.system('rm '+os.path.join(output_dir,'*.resphase'))
    return scr

