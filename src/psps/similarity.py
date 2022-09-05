import numpy as np
from psps import sario,cppsimilarity

def median_similarity(ifglist,nrow,ncol,ps0,ntilerow=1,ntilecol=1,N=20,
                      rdmin=3,rdmax=50):
    """
    For each PS pixel, median_similarity compute its phase similarity with
    nearby PS and return the median of similarity measurements

    Parameters
    ----------
    ifglist: list of str
        list of interferograms
    nrow: int
        number of rows of each interferogram
    ncol: int
        number of columns of each interferogram
    ps0 : 2D boolean array (nrow x ncol)
        PS mask: ps0[i,j] = True if pixel (i,j) is a PS pixel
    ntilerow: int (optional, default:1)
        number of rows of tiles
    ntilecol: int (optional, default:1)
        number of columns of tiles 
    N : int
        number of neighbor PS candidates used for median phase similarity
        calculation
    rdmin : int (optional, default:3)
        minimum search radius: PS pixels with distance <= rdmin to the center
        pixel are excluded for similarity calculation
    rdmax : int (optional, default:50)
        maximum search radius: PS pixels with distance > rdmax to the center
        pixel are excluded for similarity calculation

    Returns
    -------
    sim_median: 2D float array (nrow x ncol)
        median of similarity (0 at non-PS pixels)
    """
    n = len(ifglist)
    n1 = int(np.ceil(nrow/ntilerow))
    n2 = int(np.ceil(ncol/ntilecol))
    idx_list1 = np.arange(ntilerow+1)*n1
    idx_list1[-1] = nrow
    idx_list2 = np.arange(ntilecol+1)*n2
    idx_list2[-1] = ncol
    sim_median  = np.zeros((nrow,ncol),dtype=np.float32)
    for i in range(ntilerow):
        for j in range(ntilecol):
            print('processing row {}/{}, column {}/{}'.format(i+1,ntilerow,j+1,ntilecol))
            idx11 = idx_list1[i]
            idx12 = idx_list1[i+1]
            idx21 = idx_list2[j]
            idx22 = idx_list2[j+1]
            stack = sario.read_tile(ifglist,nrow,ncol,idx11,idx12,idx21,idx22)
            stack = np.angle(stack)
            sim_median_tile = cppsimilarity.median_similarity(stack,ps0,N,
                                                              rdmin,rdmax)
            sim_median[idx11:idx12,idx21:idx22] = np.array(sim_median_tile)
    return sim_median

def max_similarity(ifglist,nrow,ncol,ps0,ntilerow=1,ntilecol=1,
                   threshold=0.5,N=20,rdmin=3,rdmax=50,maxiter=20):
    """
    Find all pixels with similar phase time series to previously selected PS
    pixels based on a similarity threshold `threshold`. The function runs
    iteratively. At each iteration, non-PS pixels within a radius of rdmax
    around PS pixels identified in the previous iteration are explored. 

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
    ps0 : 2D boolean array (nrow x ncol)
        PS mask: ps0[i,j] = True if pixel (i,j) is a PS pixel
    threshold : float (optional, default:0.5)
        similarity threshold
    N : int
        maximum number of neighbor PS pixels for maximum phase similarity
        calculation
    rdmin : int
        minimum search radius: PS pixels with distance <= rdmin to the center
        pixel are excluded for similarity calculation
    rdmax : int
        maximum search radius: PS pixels with distance > rdmax to the center
        pixel are excluded for similarity calculation
    maxiter : int
        maximum iteration before PS search terminates

    Returns
    -------
    sim_max: 2D float array (nrow x ncol)
        maximum similarity to nearby PS pixels
    """
    n = len(ifglist)
    n1 = int(np.ceil(nrow/ntilerow))
    n2 = int(np.ceil(ncol/ntilecol))
    idx_list1 = np.arange(ntilerow+1)*n1
    idx_list1[-1] = nrow
    idx_list2 = np.arange(ntilecol+1)*n2
    idx_list2[-1] = ncol
    sim_max  = np.zeros((nrow,ncol),dtype=np.float32)
    for i in range(ntilerow):
        for j in range(ntilecol):
            print('processing row {}/{}, column {}/{}'.format(i+1,ntilerow,j+1,ntilecol))
            idx11 = idx_list1[i]
            idx12 = idx_list1[i+1]
            idx21 = idx_list2[j]
            idx22 = idx_list2[j+1]
            stack = sario.read_tile(ifglist,nrow,ncol,idx11,idx12,idx21,idx22)
            stack = np.angle(stack)
            sim_max_tile = cppsimilarity.max_similarity(stack,ps0,threshold,N,
                                                        rdmin,rdmax,maxiter,False)
            sim_max[idx11:idx12,idx21:idx22] = np.array(sim_max_tile)
    return sim_max

def nonps_similarity(ifglist,nrow,ncol,nonps,ntilerow=1,ntilecol=1,
                     rdmin=3,rdmax=50):
    """
    Calculates the maximum phase similarity threshold distribution for non-PS
    pixels

    Parameters
    ----------
    ifglist: list of str
        list of interferograms
    nrow: int
        number of rows of each interferogram
    ncol: int
        number of columns of each interferogram
    nonps : 2D boolean array (nrow x ncol)
        mask of calibration pixels: nonps[i,j] = True if pixel (i,j) is a
        calibration non-PS pixel
    ntilerow: int (optional, default:1)
        number of rows of tiles
    ntilecol: int (optional, default:1)
        number of columns of tiles 
    rdmin : int
        minimum search radius: PS pixels with distance <= rdmin to the center
        pixel are excluded for similarity calculation
    rdmax : int
        maximum search radius: PS pixels with distance > rdmax to the center
        pixel are excluded for similarity calculation

    Returns
    -------
    sim_max: 2D float array (nrow x ncol)
        maximum phase similarity evaluated at calibration pixels
    """
    n = len(ifglist)
    n1 = int(np.ceil(nrow/ntilerow))
    n2 = int(np.ceil(ncol/ntilecol))
    idx_list1 = np.arange(ntilerow+1)*n1
    idx_list1[-1] = nrow
    idx_list2 = np.arange(ntilecol+1)*n2
    idx_list2[-1] = ncol
    sim_max  = np.zeros((nrow,ncol),dtype=np.float32)
    for i in range(ntilerow):
        for j in range(ntilecol):
            print('processing row {}/{}, column {}/{}'.format(i+1,ntilerow,j+1,ntilecol))
            idx11 = idx_list1[i]
            idx12 = idx_list1[i+1]
            idx21 = idx_list2[j]
            idx22 = idx_list2[j+1]
            stack = sario.read_tile(ifglist,nrow,ncol,idx11,idx12,idx21,idx22)
            stack = np.angle(stack)
            sim_max_tile = cppsimilarity.max_similarity(stack,nonps,1,(rdmax*2)**2,
                                                        rdmin,rdmax,1,True)
            sim_max[idx11:idx12,idx21:idx22] = np.array(sim_max_tile)
    return sim_max

