#!/usr/bin/env python
#from pylab import *

from numpy import *
from scipy import ndimage
import time
import sys, os
import argparse, cine
import multipool
import setupfile

from cine import Sparse4D

def if_int(str):
    if not str: return None
    else: return int(str)

def eval_slice(s, N): 
    return range(*slice(*[if_int(x) for x in s.split(':')]).indices(N))
    
def work_func(args):
    frame_num, dat, setup_list, args = args
    
    setup = setupfile.eval_str_commands(setup_list)
    
    nz, ny, nx = dat.shape

#    z, y, x  = mgrid[:nz, :ny, :nx]
    y, x  = mgrid[:ny, :nx]
    
    ys = ny / float(nx)
    zs = nz / float(nx)
    #print ys, zs
    
    x0 = (x / float(nx)) * 2 - 1
    y0 = ((y / float(ny)) * 2 - 1) * ys
    #z0 = ((z / float(nz)) * 2 - 1) * zs
    
    #x = setup['x_func'](x0, y0, z0) * 0.5 + 0.5
    #y = setup['y_func'](x0, y0, z0)/ys * 0.5 + 0.5
    #z = setup['z_func'](x0, y0, z0)/zs * 0.5 + 0.5
        
    if args.order > 1:
        dat = ndimage.spline_filter(dat, args.order)

    dat2 = zeros_like(dat)
    #This avoids having to have the whole damn thing in memory        
    for zz in arange(nz):
        z0 = ((zz / float(nz)) * 2 - 1) * zs
        
        x = setup['x_func'](x0, y0, z0) * 0.5 + 0.5
        y = setup['y_func'](x0, y0, z0)/ys * 0.5 + 0.5
        z = setup['z_func'](x0, y0, z0)/zs * 0.5 + 0.5

        dat2[zz] = ndimage.map_coordinates(dat, [z*nz, y*ny, x*nx], order=args.order, prefilter=False)
    
    block = cine.make_block(dat2)
    
    #Attempt to prevent data leaks -- shouldn't really be necessary, but...
    del x, y, z, dat, dat2
    
    return (frame_num, block)
    

def crawl(p, ext=None):
    if os.path.isdir(p):
        fns = []
        for f in os.listdir(p):
            fns += crawl(os.path.join(p, f), ext)
        return fns
    else:
        if (ext is None) or os.path.splitext(p)[1] == ext:
            return [p]
        else:
            return []
            
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Correct perspective on an S4D by reinteroplating the image.  Outputs 12 bit S4d (not for viewer).')
    parser.add_argument('input', metavar='input', type=str, nargs='+',help='input files')
    parser.add_argument('-c', dest='crawl', action='store_true', default=False, help='crawl directories for files with extension .cine')
    parser.add_argument('-s', dest='setup', type=str, required=True, help='setup file, in 3dsetup format')    
    parser.add_argument('-t', dest='threads', type=int, default=8, help='number of threads (recommended <= 1 thread/2GB RAM!) [4]')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true', default=False, help='overwrite existing cines')
    parser.add_argument('-o', dest='order', type=int, default=1, help='interpolation spline order [1]')
    args = parser.parse_args()
#    print args
    
    setup_file, setup, setup_list = setupfile.get_setup(args.setup, skip_filters=True, get_string=True)
    
    if setup_file is None: raise RuntimeError("Couldn't read setupfile!")
    
    #print setup_list
    #print setupfile.eval_str_commands(setup_list)
#    sys.exit()

    if args.crawl:
        fns = []
        #print args.input
        for fn in args.input: fns += crawl(fn, ext='.sparse')
        #print fns
        #sys.exit()
    else:
        fns = args.input

    for fn in fns:
        base, ext = os.path.splitext(fn)
        ofn = base + '-pc.s4d'
        if os.path.exists(ofn) and not args.overwrite:
            print "%s\n --> %s already exists, ignoring..." % (fn, ofn)
            continue

        input = cine.Sparse(fn)
        
        #if input[0].shape != setup['frame_shape']:
        #    print "%s does not have the correct frame shape, ignoring..." % fn
        #    continue
        
        frames = len(input) // setup['cine_depth']
        #frames = 85521 // setup['cine_depth']
        #frames = 3000 // setup['cine_depth']
        print frames, len(input)
        
        header = input.header
        header['command'] = ' '.join(sys.argv)
        header['3dsetup'] = setup_list
    
        print '%s\n --> %s %3d/%3d' % (fn, ofn, 0, frames),
        sys.stdout.flush()

        #sys.exit()
        #print input.header
        #print output.header

        #sys.exit()

        x = None
        
        
        start = time.time()
            
        if args.threads == 1:
            results = []
        else:
            pool = multipool.WorkerPool(work_func, args.threads, print_result_number=False)
        
        for i in range(frames):
            raw = array([input[j + i * setup['cine_depth']] for j in setup['display_frames']])
            job = (i, raw, setup_list, args)
            
            if args.threads == 1: results.append(work_func(job))
            else: pool.send_job(job)
        
            #Attempt to fix memory error...
            del raw
            del job
        
            print '\r --> %s %3d/%3d' % (ofn, i+1, frames),
            sys.stdout.flush()
            #time.sleep(0.01)
        
        if args.threads != 1:
            results = pool.return_results()
            pool.close()

        print '\r --> %s  -->  Done in %.1f s.' % (ofn, time.time() - start)

        results.sort()

        output = cine.Sparse4D(ofn, 'w', header)
        
        for i, block in results:
            output.append_block(block)
            
        output.close()
        