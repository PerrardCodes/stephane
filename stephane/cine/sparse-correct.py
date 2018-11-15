#!/usr/bin/env python
#from pylab import *

from numpy import *
from sparse4d import Sparse4D, make_block
from scipy import ndimage
import time
import sys, os
import argparse, cine
import multipool

def if_int(str):
    if not str: return None
    else: return int(str)

def eval_slice(s, N): 
    return range(*slice(*[if_int(x) for x in s.split(':')]).indices(N))
    
def work_func(args):
    frame_num, dat, args = args
    
    #dat = (asarray(dat, 'f') * (1/255.)) ** args.gamma
    
    dat = (asarray(dat, 'f') - args.min_val)/ (2.**args.real_bpp - args.min_val)
    
    nz, ny, nx = dat.shape

    z, y, x  = mgrid[:nz, :ny, :nx]
    x = x / float(nx)
    y = y / float(ny)
    z = z / float(nz)
    zp = (z - 0.5) * nz / nx
    
    x += zp * args.shear
    x = (x - 0.5) * (1 + zp * args.perspective) + 0.5
    y = (y - 0.5) * (1 + zp * args.perspective) + 0.5
        
    dat2 = clip(ndimage.map_coordinates(dat, [z*nz, y*ny, x*nx], order=args.order), 0, 1)

    if args.background_amplitude != 1.0:
        z = (arange(dat.shape[1])[::-1] * args.scale)[newaxis, :, newaxis]
        z -= args.z0
        ba = args.background_amplitude
        
        Ie =  ba + (1-ba) * exp(-z**2 / (2. * args.sigma)**2)
        #print Ie
        
        dat2 /= Ie
        
        #I_func = lambda x: (args.background_amplitude + (1.-args.background_amplitude * exp(-(x - 72)**2/(2. * 38**2)))/69.0
        
    
    dat2u8 = (clip((dat2 * args.brighten) ** (1./args.gamma), 0, 1) * 255).astype('u1')
    
    
    return (frame_num, make_block(dat2), make_block(dat2u8))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a CINE to corrected S4Ds. Produces a "-float" and "-uint8" version.')
    parser.add_argument('input', metavar='input', type=str, nargs='+',help='input files')
    parser.add_argument('-r', dest='range', type=str, default=":", help='range of frames to convert, in python slice format [:]')
    parser.add_argument('-s', dest='shear', type=float, default=0.15, help='XZ shear correction (camera angle) [0.15]')
    parser.add_argument('-p', dest='perspective', type=float, default=0.08, help='Z scale correction (camera FOV) [0.08]')
    parser.add_argument('-o', dest='order', type=int, default=1, help='interpolation spline order [1]')
    parser.add_argument('-g', dest='gamma', type=float, default=2.0, help='gamma correction, only applies to uint8 version [2.0]')
    parser.add_argument('-t', dest='threads', type=int, default=1, help='number of threads (recommended <= 1 thread/2GB RAM!) [1]')
    parser.add_argument('-d', dest='depth', type=int, default=300, help='stacks per volume')
#    parser.add_argument('-M', dest='top_clip', type=int, default=1024, help='max val of rescale, only applies to uint8 version (top clip) [1024]')
    parser.add_argument('-b', dest='brighten', type=float, default=1., help='amount to brighten image before gamma/uint8 conversion (> 1 results in clipping) [1.0]')
    parser.add_argument('-m', dest='min_val', type=int, default=40, help='min val of rescale (bottom clip) [40]')
    parser.add_argument('-D', dest='displayframes', type=str, default="14:270", help='range of z frames to save in volume [:]')
#    parser.add_argument('-h', dest='histogram', type=bool, default=False, help='display a histogram and exit')
#    parser.add_argument('-o', dest='output', type=str, default="%s.s4d", help='output filename')
    parser.add_argument('--skip', dest='skip', type=int, default=0, help='skip this many frames in the file; used to fix sync offsets [0]')
    parser.add_argument('-S', dest='scale', type=float, default=0.6, help='mm per voxel at center [0.6]')
    parser.add_argument('-z', dest='z0', type=float, default=60., help='z offset of center of laser brightness in mm [60]')
    parser.add_argument('--sigma', dest='sigma', type=float, default=38.0, help='gaussian width (sigma) of laser brightness function in mm [38]')
    parser.add_argument('-B', dest='background_amplitude', type=float, default=0.22, help='fraction of max laser brightness that is background instead of gaussian peak; set to 1 for no brightness correction, 0 for fitting to a straight gaussian [default: 0.22 != 0 because of spherical-aberation/best fit]')
    args = parser.parse_args()
 
    for fn in args.input:
        base, ext = os.path.splitext(fn)

        if ext.lower().startswith('.tif'):
            sys.setrecursionlimit(10**5)
            source = cine.Tiff(fn)    
        else:
            source = cine.Cine(fn)
            
         
        try:
            args.real_bpp = source.real_bpp
        except:
            args.real_bpp = soure[0].dtype.itemsize * 8
        
        max_frame = (len(source) - args.skip) // args.depth
        saved_frames = eval_slice(args.range, max_frame)
        
        print '%s (%d frames):' % (fn, len(saved_frames))
 
        frame_offsets = array(eval_slice(args.displayframes, args.depth))

        
        if args.threads == 1:
            results = []
            start = time.time()
        else:
            pool = multipool.WorkerPool(work_func, args.threads)
        
        for i, frame_num in enumerate(saved_frames):
            raw_frame = array([source[j] for j in (frame_offsets + frame_num * args.depth + args.skip)])
            
            if args.threads == 1:
                results.append(work_func((i, raw_frame, args)))
                print '%4d' % i,
                sys.stdout.flush()
            else:
                pool.send_job((i, raw_frame, args))
        
        if args.threads != 1:
            results = pool.return_results()
            pool.close()
        else:
            print '-> Done in %.1f s.' % (time.time() - start)


        results.sort()
        
        header = vars(args)
        header['input'] = fn
        header['command'] = ' '.join(sys.argv)

        
        output_f = Sparse4D(base + '-float.s4d', 'w', header)
        output_u8 = Sparse4D(base + '-uint8.s4d', 'w', header)
                
        for i, block_f, block_u8 in results:
            output_f.append_block(block_f)
            output_u8.append_block(block_u8)
            
        output_f.close()
        output_u8.close()
