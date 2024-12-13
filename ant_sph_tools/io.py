"""ant_sph_tools.io -- Read data in TIRCA .sph format"""
import numpy as np  


def _generate_j(mmax: int, nmax: int) -> np.ndarray:
    """ Generate array of compressed mode coefficients j.

    Args:
        mmax (int): Maximum azimuthal mode order.
        nmax (int) Maximum elevation mode order.

    Returns:
        j (np.ndarray): 1D array of compressed mode coeffs.
                        len: 2*(nmax * (nmax+1) + mmax-1) + 2
                        dtype: int32
    """
    m = 0
    s = 2
    
    j = np.zeros(2 * (nmax * (nmax + 1) + mmax - 1) + 2, dtype='int32')
    j3D = np.zeros_like(j, dtype=object)
    
    cnt = 0
    for n in range(1, nmax+1):
        for s in (1, 2):
            j[cnt] =  2 *(n*(n+1)+ m-1) + s 
            j3D[cnt] =  f"{s} {m} {n}"
            cnt += 1
    
    for mmode in range(1, mmax + 1):  
        for n in range(abs(mmode), nmax + 1):  
            for m in range(-mmode, mmode + 1, 2 * mmode):  
               for s in range(1, 3):
                    j[cnt] =  2 *(n*(n+1)+ m-1) + s
                    j3D[cnt] =  f"{s} {m} {n}"  
                    cnt += 1
    return j, j3D


def read_sph(fn: str) -> dict:  
    """  Reads in spherical mode data in TIRCA .sph format. 
    
    Reads spherical mode Q coefficients (See GRASP documentation).  
    
   Args:  
       fn (str): File name of the .sph file  
    
   Returns:
       sph_data (dict): Dictionary with following entries:
           Q (numpy array): Spherical mode coefficients  
           j (numpy array): Compressed mode indices   
           j3D (numpy array): 3D mode indices   
           mmax (int): Maximum azimuthal mode number  
           nmax (int): Maximum elevation mode number  
    """  
  
    # Read the file, skipping the first two lines  
    with open(fn, 'r') as f:  
        lines = f.readlines()[2:]  
    
    # Parse the header lines  
    header = lines[:6]
    data = lines[6:]
    
    nmax = int(header[0].split()[2])  
    mmax = int(header[0].split()[3])  
    
    # Initialize arrays  
    Q = np.zeros((2 * (nmax * (nmax + 1) + mmax - 1) + 2), dtype=np.complex128)  
    j, j3D = _generate_j(mmax, nmax)
    
    idx = 0
    for line in data:
        n_l = len(line.split())
        if n_l == 4:
            Q[2*idx:2*idx+2] = np.fromstring(line, sep=' ', dtype='float64').view('complex128')
            idx += 1

    # Remap Q by compressed mode coefficient order
    Q[j-1] = np.copy(Q)
    
    return {'Q': Q, 'j': j, 'j3D': j3D, 'mmax': mmax, 'nmax': nmax}