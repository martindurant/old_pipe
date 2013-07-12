""" 
Dynamic library interface to the old_pipe motion correction/registration code 
based on intensity-variation optic flow equations (See Martel, Froh et al.).

Provides the following improvements over the original pure C code:

- instead of accessing analyze files on-disc, pass in-memory arrays directly,
enhancing efficiency 

- ability to get/set the mapping, so that one regisration solution can be 
applied to other data, e.g., to a whole dynamic sequence

- call-backs to assess progress

- remove all argument parsing code and unnecessary code (e.g., file output)
to simplify and make more readable.

This work required de-constructing the C-code package, which is copied into 
this directory and modified in-place. Should produce identical results to
prog.registration.old_pipe for the same parameters
"""

import os
import numpy as np
import ctypes as ct
here = os.path.dirname(os.path.abspath(__file__)) + os.sep

c_fl_p = ct.POINTER(ct.c_float)
c_i_p = ct.POINTER(ct.c_int)
nd = np.ctypeslib.ndpointer(np.float32,3,None,'C')
_lib = ct.cdll.LoadLibrary(here+'_old_pipe.so')

_lib.make_fixed.argtypes = [nd,ct.c_int,ct.c_int,ct.c_int,
                            ct.c_float,ct.c_float,ct.c_float]
_lib.make_moving.argtypes = [nd,ct.c_int,ct.c_int,ct.c_int,
                            ct.c_float,ct.c_float,ct.c_float]
_lib.make_mask.argtypes = [nd,ct.c_int,ct.c_int,ct.c_int,
                            ct.c_float,ct.c_float,ct.c_float]

"""for each Image structure, the first three elements are ints giving
the array size, and the next is a pointer to the data; the ints are
int32, and the pointer to the data is at fixed+16 bytes, also an int32
"""

"""Suspect that breakage is in the removal of if clauses dealing with
ROI, which was referenced in setting up the number of nodes and the
extremes of the ranges to use.
Should repair them, and simply provide an ROI that is all 1.0"""

def put_data(fixed,moving,vox):
    """Make the global variables fixed,moving and mask point to the numpy array
    data given, which must be float32, C-shape and 3d.
    vox are the voxel dimentsions (mm floats), and mask may be disbanded, since
    it will always be all 1 anyway.
    """
    assert moving.shape == fixed.shape
    z,y,x = fixed.shape
    k,j,i = vox
    _lib.make_fixed(fixed,x,y,z,i,j,k)
    _lib.make_moving(moving,x,y,z,i,j,k)
    mask = np.ones_like(fixed)
    _lib.make_mask(mask,x,y,z,i,j,k)
    return mask

def get_map():
    """Return all the information in the current map, stored in global variables
    in _lib: 
    
    - a set of integers (including nnodes, the number of nodes)
    
    - a set of binary strings of length nnodes float32 values (nnodes*4 bytes)
    representing the arrays.
    
    The values of these arrays can be changed in place.
    """
    nnodes, XP, YP, ZP, SIDE, MCOLS, MROWS, MSLICES = [ct.c_int.in_dll(_lib,u).value for u in
        ('NUMBER_OF_NODES','XP', 'YP', 'ZP', 'SIDE', 'MCOLS', 'MROWS', 'MSLICES')]
    ints = [nnodes, XP, YP, ZP, SIDE, MCOLS, MROWS, MSLICES]
    strs = [ct.string_at(c_fl_p.in_dll(_lib,u)) for u  in 
        ('XC', 'YC', 'ZC', 'MAP', 'UT', 'VT', 'WT', 'AT')]
    return ints,strs

def put_map(ints,strs):
    """Return all the information in the current map, stored in global variables
    in _lib: 
    
    - a set of integers (including nnodes, the number of nodes)
    
    - a set of float arrays of length nnodes float32 values (nnodes*4 bytes)
    representing the arrays.
    """
    for (label,value) in zip(('NUMBER_OF_NODES','XP', 'YP', 'ZP', 'SIDE', 'MCOLS', 'MROWS', 'MSLICES'),ints):
        ct.c_int.in_dll(_lib,label).value = value
    for (label,value) in zip(('XC', 'YC', 'ZC', 'MAP', 'UT', 'VT', 'WT', 'AT'),strs):
        st = ct.create_string_buffer(value)
        ct.c_int.in_dll(_lib,label).value = ct.byteref(st)

def test():
    from prog import volumes
    v = volumes.load('/home/durant/501/12778',[3])
    v = v.select([0,1]).segment(left=1)
    print "selected"
    mask = put_data(v[0].astype(np.float32),v[1].astype(np.float32),v.vox)
    print "putted"
    _lib.go()
    return mask,v