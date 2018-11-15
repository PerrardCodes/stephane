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
    
    dat = (asarray(dat, 'f') / 255.) ** args.gamma
    
    nz, ny, nx = dat.shape

    z, y, x  = mgrid[:nz, :ny, :nx]
    x = x / float(nx)
    y = y / float(ny)
    z = z / float(nz)
    zp = (z - 0.5) * nz / nx
    
    x += zp * args.shear
    x = (x - 0.5) * (1 + zp * args.perspective) + 0.5
    y = (y - 0.5) * (1 + zp * args.perspective) + 0.5
        
    dat2 = ndimage.map_coordinates(dat, [z*nz, y*ny, x*nx], order=args.order)
    
    return (frame_num, make_block(dat2))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Correct perspective on an S4D by reinteroplating the image.  Outputs float S4D (not for use in the viewer!).')
    parser.add_argument('input', metavar='input', type=str, nargs='+',help='input files')
    parser.add_argument('-r', dest='range', type=str, default=":", help='range of frames to convert, in python slice format [:]')
    parser.add_argument('-s', dest='shear', type=float, default=0.15, help='XZ shear correction (camera angle) [0.15]')
    parser.add_argument('-p', dest='perspective', type=float, default=0.08, help='Z scale correction (camera FOV) [0.08]')
    parser.add_argument('-o', dest='order', type=int, default=1, help='interpolation spline order [1]')
    parser.add_argument('-g', dest='gamma', type=float, default=2.0, help='gamma correction [2.0]')
    parser.add_argument('-t', dest='threads', type=int, default=1, help='number of threads (recommended <= 1 thread/2GB RAM!) [1]')
    args = parser.parse_args()

    for fn in args.input:
        base, ext = os.path.splitext(fn)
        input = Sparse4D(fn)
        frames = eval_slice(args.range, len(input))
        ofn = base + '-pcorrect' + ext
        
        header = input.header
        header['command'] = ' '.join(sys.argv)
        header['shear'] = args.shear
        header['perspective'] = args.perspective
        header['gamma'] = args.gamma
    
        print '%s -> %s (%d frames)' % (fn, ofn, len(frames))
        #print input.header
        #print output.header

        #sys.exit()

        x = None
        
        if args.threads == 1:
            results = []
            start = time.time()
        else:
            pool = multipool.WorkerPool(work_func, args.threads)
        
        for i in frames:
            if args.threads == 1:
                results.append(work_func((i, input[i], args)))
                print '%4d' % i,
                sys.stdout.flush()
            else:
                pool.send_job((i, input[i], args))
        
        if args.threads != 1:
            results = pool.return_results()
            pool.close()
        else:
            print '-> Done in %.1f s.' % (time.time() - start)


        results.sort()

        
        output = Sparse4D(ofn, 'w', header)
        
        for i, block in results:
            output.append_block(block)
            
        output.close()
        
        #    start = time.time()
        #    
        #    dat = (asarray(input[i], 'f') * (1./255)) ** args.gamma
        #    if x is None:
        #        nz, ny, nx = dat.shape
        #
        #        z, y, x  = mgrid[:nz, :ny, :nx]
        #        x = x / float(nx)
        #        y = y / float(ny)
        #        z = z / float(nz)
        #        zp = (z - 0.5) * nz / nx
        #        
        #        x += zp * args.shear
        #        x = (x - 0.5) * (1 + zp * args.perspective) + 0.5
        #        y = (y - 0.5) * (1 + zp * args.perspective) + 0.5
        #        
        #    dat2 = ndimage.map_coordinates(dat, [z*nz, y*ny, x*nx], order=args.order)
        #    
        #    output.append_array(dat2)
        #    print '%d (%.1fs)...' % (i, time.time() - start),
        #    sys.stdout.flush()
        #print 'done.'

#print datf.header
#dat = datf[50][:, 128:, :]
#dat = (asarray(dat, 'f')*(1./255)) **2
#
#
#nz, ny, nx = dat.shape
#
#TEXTURE_SHEAR = 0.15
#TEXTURE_PERSPECTIVE_SCALE = 0.08
#GAMMA = 2.
#
#
#
#zp = (z - 0.5) * nz / nx
#
#x += zp * TEXTURE_SHEAR
#x = (x - 0.5) * (1 + zp * TEXTURE_PERSPECTIVE_SCALE) + 0.5
#y = (y - 0.5) * (1 + zp * TEXTURE_PERSPECTIVE_SCALE) + 0.5
#
#dat2 = ndimage.map_coordinates(dat, [z*nz, y*ny, x*nx], order=1)
#print time.time() - start



#for i in range(3):
#    subplot(220+i)
#    
#    dats = dat2.sum(i)
#    datss = dats.sum()
#    y, x = ogrid[:dats.shape[0], :dats.shape[1]]
#    
#    yavg = (y * dats).sum() / datss
#    xavg = (x * dats).sum() / datss
#    
#    imshow(sqrt(dats))
#    plot(xavg, yavg, 'rx')
#    
#
#vert_sum = dat2.sum(0)
#
#
#z, y, x  = ogrid[:nz, :ny, :nx]
#x = x / float(nx)
#y = y / float(ny)
#z = z / float(nz)
#
#
#
#    
#show()