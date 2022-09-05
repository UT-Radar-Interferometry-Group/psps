import glob
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

def np2rgb(m,cmap='jet',vmin=None,vmax=None):
    """
    convert a 2D numpy array to a rgba matrix, pixels with nan value are
    considered transparent
    
    Parameters
    ----------
    m: 2D array
        The image data.
    cmap: str or Colormap, optional (default: 'jet')
        A Colormap instance or registered colormap
    vmin,vmax: float, optional
        vmin and vmax set the color scaling by fixing the values that map
        to the colormap color limits
    """

    if vmin is None:
        vmin = np.percentile(m,1)
    if vmax is None:
        vmax = np.percentile(m,99)
    norm = plt.Normalize(vmin=vmin,vmax=vmax)
    temp = np.copy(m)
    temp[np.isnan(m)] = vmin
    if isinstance(cmap,str):
        colors = getattr(plt.cm,cmap)(norm(temp))
    else:
        colors = cmap(norm(temp))
    colors[:,:,3] = ~np.isnan(m)
    return colors

def savem_to_image(m,fname,cmap='jet',vmin=None,vmax=None):
    """
    save a 2D numpy array as an image file
    
    Parameters
    ----------
    m: 2D array
        The image data.
    fname: str
        The name of the output image file
    cmap: str or Colormap, optional (default: 'jet')
        A Colormap instance or registered colormap
    vmin,vmax: float, optional
        vmin and vmax set the color scaling by fixing the values that map
        to the colormap color limits
    """

    colors = np2rgb(m,cmap,vmin,vmax)
    plt.imsave(fname,colors)

def read_ifg(fname,nrow,ncol):
    ifg = np.fromfile(fname,dtype='float32')
    ifg = np.reshape(ifg,(nrow,2*ncol))
    ifg = ifg[:,0:2*ncol:2] + 1j * ifg[:,1:2*ncol:2]
    return ifg

def save_ifg(ifg,fname):
    naz,nrg = ifg.shape
    tosave = np.zeros((naz,nrg*2),dtype='float32')
    tosave[:,0:2*nrg:2] = np.real(ifg)
    tosave[:,1:2*nrg:2] = np.imag(ifg)
    f = open(fname,'w')
    tosave.tofile(f)
    f.close()

def read_tile(ifglist,nrow,ncol,idx11,idx12,idx21,idx22):
    n = len(ifglist)
    s = np.zeros((n,idx12-idx11,idx22-idx21),dtype=np.complex64)
    for i,ifgfile in enumerate(ifglist):
        res = read_ifg(ifgfile,nrow,ncol)
        s[i,:,:] = res[idx11:idx12,idx21:idx22]
    return s

def cpxlooks(imgbg,nlookaz,nlookrg):
    naz,nrg = imgbg.shape
    newnaz = np.floor(naz/nlookaz).astype(int)
    newnrg = np.floor(nrg/nlookrg).astype(int)
    imgaz = np.zeros((newnaz,nrg),dtype=imgbg.dtype)
    imgsm = np.zeros((newnaz,newnrg),dtype=imgbg.dtype)
    if nlookaz>1:
        for i in np.arange(0,newnaz):
            imgaz[i,:] = np.sum(imgbg[i*nlookaz:(i+1)*nlookaz,:],axis=0)
    else:
        imgaz = imgbg
    if nlookrg>1:
        for i in np.arange(0,newnrg):
            imgsm[:,i] = np.sum(imgaz[:,i*nlookrg:(i+1)*nlookrg],axis=1)
    else:
        imgsm = imgaz
    return imgsm/nlookrg/nlookaz

def powlooks(imgbg,nlookaz,nlookrg):
    return np.abs(cpxlooks(np.abs(imgbg)**2,nlookrg,nlookaz))

def multilooks(imgfile,outfile,nr,nc,nrlook,nclook,chuncksize=1e9):
    """
    Take multilook of a large SLC image (>10 G) to get the multilooked
    amplitude

    Parameters
    ----------
    imgfile : string
        slcfile to read
    outfile : string
        output amplitude file to write
    nr : int
        number of rows of the image
    nc : int
        number of columns of the image
    nrlook : int
        number of looks to take in the row direction
    nclook : int
        number of looks to take in the column direction
    chunckisize : int, optional (default 1e9 [bit])
        number of bytes to read each time
    """

    nrpatch = int(int(chuncksize/(nc*4))//nrlook*nrlook)
    nr_read = 0
    with open(outfile,'wb') as fout:
        with open(imgfile,'r') as f:
            while nr_read < nr-nrlook+1:
                rcount = int(np.minimum(nrpatch,nr-nr_read))
                nr_read = nr_read + rcount
                patch = np.fromfile(f, dtype=np.float32, count=rcount*nc*2)
                patch = np.reshape(patch,(rcount,2*nc))
                patch = patch[:,0:2*nc:2] + 1j * patch[:,1:2*nc:2]
                patch = np.sqrt(powlooks(patch,nrlook,nclook)).astype(np.float32)
                fout.write(patch.tobytes())

def ifg_cmap(n=256):
    """
    A colorbar used for displaying interferograms

    Parameters
    ----------
    n: int, optional (default: 256)
        number of discrete colors

    Returns
    -------
    cmap: Colormap instance
        A colormap

    """
    J = [[1.0,0.7,0.7],
	 [1.0,1.0,0.4],
	 [0.4,1.0,1.0],
	 [1.0,0.4,1.0],
	 [1.0,0.7,0.7]]
    J = np.array(J)
    c = np.ones((n,4))
    for i in range(3):
        c[:,i] = np.interp(np.linspace(0,1,n),np.linspace(0,1,5),J[:,i])
    cmap = ListedColormap(c)
    return cmap

ifgcmap = ifg_cmap(256)
