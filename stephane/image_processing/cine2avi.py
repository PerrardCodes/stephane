#!/usr/bin/env python

from numpy import *
import os, sys
import argparse
import time
import matplotlib.pyplot as plt
import glob

import image#, ImageDraw, ImageFont

import stephane.Image_processing.cine as cine
import stephane.Image_processing.mjpeg as mjpeg
import stephane.Image_processing.tiff as tiff

#import stephane.Image_processing as Im_proc

"""
Cine2avi :
Convert CINE file(s) to an AVI.  Also works on TIFFs.
Can either be used in command line (python cine2avie cines [options]) 
or from a python interface (import cine2avi, cine2avi.make(cinename))
See also gen_movies.py for multi movies method, and other compressing options (to mp4)
"""

script_dir = os.path.dirname(os.path.realpath(__file__))

def gen_parser():    
    parser = argparse.ArgumentParser(description="Convert CINE file(s) to an AVI.  Also works on TIFFs.")
    parser.add_argument('cines', metavar='cines', type=str, nargs='+', help='input cine file(s), append [start:end:skip] (python slice notation) to filename to convert only a section')
    parser.add_argument('-range', metavar='range', type=str, nargs='+', help='input cine file(s), append [start:end:skip] (python slice notation) to filename to convert only a section')
    parser.add_argument('-o', dest='output', type=str, default='%s.avi', help='output filename, may use %%s for input filename w/o extension or %%d for input file number')
    parser.add_argument('-g', dest='gamma', type=float, default=1., help='gamma of output, assumes input is gamma 1, or: I -> I**(1/gamma); use 2.2 to turn linear image to "natural" for display [default: 1]')
    parser.add_argument('-f', dest='framerate', type=int, default=30, help='frames per second [default: 30]')
    parser.add_argument('-q', dest='quality', type=int, default=75, help='JPEG quality (0-100) [default: 75]')
    parser.add_argument('-c', dest='clip', type=float, default=0, help='histogram clip, clip the specified fraction to pure black and pure white, and interpolate the rest; applied before gamma; recommended is 1E-4 - 1E-3 [default: 0]')
    parser.add_argument('-s', dest='hist_skip', type=int, default=10, help='only check every Nth frame when computing histogram clip [default: 5]')
    parser.add_argument('-r', dest='rotate', type=int, default=0, help='amount to rotate in counterclockwise direction, must be multiple on 90 [default: 0]')
    parser.add_argument('-t', dest='timestamp', default=False, action='store_true', help='write a timestamp on each frame')
    parser.add_argument('--font', dest='font', type=str, default=os.path.join(script_dir, 'Helvetica.ttf'), help='font (ttf) used for the timestamp')
    parser.add_argument('--ts', dest='ts', type=int, default=40, help='text size for timestamp in pixels [25]')
    parser.add_argument('--tx', dest='tx', type=int, default=25, help='timestamp x origin [25]')
    parser.add_argument('--ty', dest='ty', type=int, default=50, help='timestamp y origin [25]')
    parser.add_argument('--td', dest='td', type=int, default=None, help='digits in timestamp [determined from framerate]')
    parser.add_argument('--tb', dest='tb', type=int, default=255, help='Test brightness, 0-255 [255=white].')
    parser.add_argument('--type', dest='type', type=str, default='single', help='Type of movie. single : one avi for each cine file. multiple : one avi panel for all the cine files given')
    parser.add_argument('--crop', dest='rect', type=str, default=None, help='Crop of the initial movies')

#    print(parser)   
    args = parser.parse_args()
    
    print(args)
    return args

    
def fmt_time(t):
    return '%d:%02d:%02d' % (int(t/3600.), int(t/60.)%60, int(t)%60)

#------------------------------------------------------------------------------
# Class for printing pretty status counters
#------------------------------------------------------------------------------

class StatusPrinter(object):
    def __init__(self, count, msg='Calculating...', min_time=1./30):
        if hasattr(count, '__len__'):
            self.data = count
            self.count = len(count)
        elif hasattr(count, '__iter__'):
            self.data = list(count)
            self.count = len(self.data)
        else:
            self.count = count
        
        self.msg = msg
        self.current = -1
        self.start = time.time()
        self.last_time = time.time()
        self.max_len = 0
        self.min_time = min_time
        self.extra_msg = ''

        if not hasattr(self, 'data'): self.update()
    
    def message(self, msg):
        self.extra_msg = msg
    
    def print_msg(self, msg, final=False):
        if len(msg) > self.max_len: self.max_len = len(msg)
        
        msg = ('%%-%ds' % self.max_len) % msg
        
        if final: print '\r' + msg
        else: print '\r' + msg,
            
        sys.stdout.flush()        


    def update(self, extra_msg=''):
        self.current += 1
        t = time.time()
        
        if self.current < self.count and (t - self.last_time) < self.min_time:
            return None
        
        self.last_time = t
        
        percent = int(100 * self.current / self.count + 0.5)
        
        if not extra_msg: extra_msg = self.extra_msg
        if extra_msg: extra_msg = ' [%s]' % extra_msg
        
        elapsed = t - self.start

        if self.current == self.count:
            self.print_msg(self.msg + ' done in %s. ' % fmt_time(elapsed) + extra_msg, final=True)
        elif self.current <= 0:
            self.print_msg(self.msg + ' %2d%% ' % percent + extra_msg)
        elif self.current < self.count:
            est = elapsed / self.current * (self.count - self.current)
                
            self.print_msg(self.msg + ' %2d%% (%s remaining) ' % (percent, fmt_time(est)) + extra_msg)

    def __iter__(self):
        return self
    
    def next(self):
        self.update()
        
        if self.current < self.count:
            if hasattr(self, 'data'): return self.data[self.current]
            else: return self.current
        else:
            raise StopIteration

def noneint(s):
    return None if not s else int(s)

def single_avi(args):
    if args.timestamp:
        print(os.path.join(script_dir, 'Helvetica.ttf'))
        #font = ImageFont.truetype(args.font, args.ts)

    for i, fn in enumerate(args.cines):
        fn = fn.strip()
    
        frame_slice = slice(None)
        if '[' in fn:
            if fn[-1] == ']':
                fn, s = fn.split('[')
                liste = [int(i) for i in s[:-1].split(':')]
                print(liste)
                #try:
                frame_slice = slice(liste[0],liste[1],liste[2])
#                frame_slice = slice(*map(noneint, liste))
                #frame_slice = slice(liste)
                #except:
                #   raise ValueError("Couldn't convert '[%s' to slice notation" % s)

            else:
                print "Warning, found '[' in input, but it didn't end with ']', so I'll assume you didn't mean to give a frame range."
    
        base, ext = os.path.splitext(fn)
        ext = ext.lower()
    
        if not os.path.exists(fn):
            print "File %s not found, ignoring." % fn
            continue
    
        output = args.output
        if '%s' in args.output: output = output % base
        elif '%' in args.output: output = output % i
    
        bpp = None
    
        if ext in ('.cin', '.cine'):
            inp = cine.Cine(fn)
            bpp = inp.real_bpp
            if bpp < 8 or bpp > 16: bpp = None #Just in case

            td = args.td if args.td else int(ceil(log10(inp.frame_rate)))
           # frame_text = lambda i: 't: %%.%df s' % td % (i/float(input.frame_rate))
        
          #  frame_text = lambda i: 't: %%.%d s' % (inp.get_time(i))
            t0 = 0.
            frame_text = lambda i: 't: %%.%df s, ' % td % (inp.get_time(i)-t0) +  'Dt: %f ms' % (round(((inp.get_time(i)-inp.get_time(i-1))*1000)*10)/10)# 't: %f s \n Dt :  ms' % 
        
        elif ext in ('.tif', '.tiff'):
            inp = tiff.Tiff(fn)
            frame_text = lambda i: str(i)
                
        bpps = inp[0].dtype.itemsize * 8
        if bpp is None: bpp = bpps
    
        frames = range(*frame_slice.indices(len(inp)))
    
        if args.clip == 0:
            mappa = linspace(0., 2.**(bpps - bpp), 2**bpps)
        else:
            counts = 0
            bins = arange(2**bpps + 1)
        
            for i in frames[::args.hist_skip]:
                c, b = histogram(inp[i], bins)
                counts += c
        
            counts = counts.astype('d') / counts.sum()
            counts = counts.cumsum()
        
            bottom_clip = where(counts > args.clip)[0]
            if not len(bottom_clip): bottom_clip = 0
            else: bottom_clip = bottom_clip[0]

            top_clip = where(counts < (1 - args.clip))[0]
            if not len(top_clip): top_clip = 2**bpps
            else: top_clip = top_clip[-1]

            #print bottom_clip, top_clip
            #import pylab
            #pylab.plot(counts)
            #pylab.show()
            #sys.exit()

            m = 1. / (top_clip - bottom_clip)
            mappa = clip(-m * bottom_clip + m * arange(2**bpps, dtype='f'), 0, 1)
            
        mappa = mappa ** (1./args.gamma)
    
        mappa = clip(mappa * 255, 0, 255).astype('u1')

        #print '%s -> %s' % (fn, output)
    
        ofn = output
        output = mjpeg.Avi(output, framerate=args.framerate, quality=args.quality)
    
        if args.rect is not None:
            rect = [int(i) for i in args.rect[1:-1].split(':')] 
            print(rect)
        
        #print frames
        for i in StatusPrinter(frames, os.path.basename(ofn)):
            frame = inp[i]
            if args.rotate:
                frame = rot90(frame, (args.rotate%360)//90)
                        
            frame = mappa[frame]       
        #
            if args.rect==None:
                frame = asarray(frame)
            else:       
                frame = asarray(frame[rect[0]:rect[1],rect[2]:rect[3]])

            if args.timestamp:
                frame = image.fromarray(frame)
               #draw = ImageDraw.Draw(frame)
                #draw.text((args.tx, args.ty), frame_text(i), font=font, fill=args.tb)
        
        
            frame = asarray(frame)
    #        print(type(frame))
            output.add_frame(frame)
        
        output.close()
       
def multiple_avi(args): 
    """
    Generate avi from multiple cinefiles
    
    INPUT
    -----
    args : Namespace object
        contain each of the arguments used to generate an avi. See the top of this file for details, or from a terminal cine2avi --help
    OUTPUT
    -----
    None
    """
    if args.timestamp:
        print(os.path.join(script_dir, 'Helvetica.ttf'))
        #font = ImageFont.truetype(args.font, args.ts)

    files = args.cines
    #for i, fn in enumerate(files):
    fn = files[0].strip()
    frame_slice = slice(None)
    if '[' in fn:
        if fn[-1] == ']':
            fn, s = fn.split('[')
            try:
                frame_slice = slice(*map(noneint, s[:-1].split(':')))
            except:
                raise ValueError("Couldn't convert '[%s' to slice notation" % s)

        else:
            print "Warning, found '[' in input, but it didn't end with ']', so I'll assume you didn't mean to give a frame range."
    
    base, ext = os.path.splitext(fn)
    ext = ext.lower()
    
    if not os.path.exists(fn):
        print "File %s not found, ignoring." % fn
      #  continue
    
    output = args.output
    if '%s' in args.output: output = output % base
    elif '%' in args.output: output = output % i
    
    base, extout = os.path.splitext(output)
    output = output[:-6] + '_multiple' + extout
    print(output)
    
    bpp = None
    
    inp = [None for fn in files]
    
    if ext in ('.cin', '.cine'):
        inp_ref = cine.Cine(fn)
        
        for i,fn in enumerate(files):
            print(fn)
            fn = fn.strip()
            inp[i] = cine.Cine(fn)

        bpp = inp_ref.real_bpp
        if bpp < 8 or bpp > 16: bpp = None #Just in case
        td = args.td if args.td else int(ceil(log10(inp_ref.frame_rate)))
        t0 = 0.
        frame_text = lambda i: 't: %%.%df s, ' % td % (inp_ref.get_time(i)-t0) +  'Dt: %f ms' % (round(((inp_ref.get_time(i)-inp_ref.get_time(i-1))*1000)*10)/10)# 't: %f s \n Dt :  ms' % 
        
    elif ext in ('.tif', '.tiff'):
        inp = tiff.Tiff(fn)
        frame_text = lambda i: str(i)
                
    bpps = inp_ref[0].dtype.itemsize * 8
    if bpp is None: bpp = bpps
    
    lengths = [len(i) for i in inp]
    print('Movie lengths :'+str(lengths))
    Nmax = min(lengths)
    frames = range(*frame_slice.indices(Nmax))
    
    if args.clip == 0:
        mappa = linspace(0., 2.**(bpps - bpp), 2**bpps)
    else:
        counts = 0
        bins = arange(2**bpps + 1)
        
        for i in frames[::args.hist_skip]:
            c, b = histogram(inp_ref[i], bins)
            counts += c
        
        counts = counts.astype('d') / counts.sum()
        counts = counts.cumsum()
        
        bottom_clip = where(counts > args.clip)[0]
        if not len(bottom_clip): bottom_clip = 0
        else: bottom_clip = bottom_clip[0]

        top_clip = where(counts < (1 - args.clip))[0]
        if not len(top_clip): top_clip = 2**bpps
        else: top_clip = top_clip[-1]

            #print bottom_clip, top_clip
            #import pylab
            #pylab.plot(counts)
            #pylab.show()
            #sys.exit()
        m = 1. / (top_clip - bottom_clip)
        mappa = clip(-m * bottom_clip + m * arange(2**bpps, dtype='f'), 0, 1)
            
    mappa = mappa ** (1./args.gamma)
    
    mappa = clip(mappa * 255, 0, 255).astype('u1')

        #print '%s -> %s' % (fn, output)
    
    ofn = output
    output = mjpeg.Avi(output, framerate=args.framerate, quality=args.quality)
    
    if args.rect is not None:
        rect = [int(i) for i in args.rect[1:-1].split(':')] 
        print(rect)

        #print frames
        
    for i in StatusPrinter(frames, os.path.basename(ofn)):
        frame = [None for p in inp]
        for p,inpp in enumerate(inp):
            framep = inpp[i]
            if args.rotate:
                framep = rot90(framep, (args.rotate%360)//90)
            
            framep = mappa[framep]       
    #        print(type(frame))
            
        
            if args.rect==None:
                frame[p] = asarray(framep)
            else:       
                frame[p] = asarray(framep[rect[0]:rect[1],rect[2]:rect[3]])
                
        fra = concatenate(tuple(frame), axis=0)
        if args.timestamp:
            fra = image.fromarray(fra)
            #draw = ImageDraw.Draw(fra)
            #draw.text((args.tx, args.ty), frame_text(i), font=font, fill=args.tb)
            
        fra = asarray(fra)
        
    #        print(type(frame))
       # plt.imshow(fra)
    #    plt.show()
    #    input()
        output.add_frame(fra)
        
    output.close()
    
def make(cinename,**kwargs):
    """
    generate avi file for each given cine file.
    INPUT
    -----
    cinename : str
        cine file name. Can contain regular expression.
    kwargs : dict of key word arguments. 
        See on top of this file for options
        type : either multiple or single. multiple : generate one avi from all the matching cine. single : one avi for each cine file
    OUTPUT
    -----
    None
    """
    cineList = glob.glob(cinename)
    
    args = set_parser(cineList)
#    args = python_parser()
    print(args.cines)
    argList = dir(args)
    for kwarg in kwargs.keys():
        if kwarg in argList:
            print('Value found, set to the value given in kwarg')
            setattr(args,kwarg,kwargs[kwarg])    
            print(kwarg,getattr(args,kwarg))
    print(args.rect)
    main(args)

def set_parser(cineList):
    """
    Generate the default arguments for generating an avi from a cine file
    INPUT 
    -----
    cineList : str list
        list of cine path
    OUTPUT
    -----
    args : Namespace object
        contain each of the arguments used to generate an avi. See the top of this file for details, or from a terminal cine2avi --help
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    #default values of the parameters
    argList = dict(cines=cineList, clip=0, font=os.path.join(script_dir, 'Helvetica.ttf'), framerate=30, gamma=1.0, hist_skip=10, output='%s.avi', quality=75, rect=None, rotate=0, tb=255, td=None, timestamp=False, ts=40, tx=25, ty=50, type='single')
    args = type('Parser',(object,),{})
    
    for key in argList.keys():
        setattr(args,key,argList[key])

    return args

def main(args):
#    print(args.getattr)
#    print(args.__getattr__())
#    print(getattr(args))
#    for arg in kwargs:
#        if arg in args.keys():
    
    if args.type=='single':
        single_avi(args)
    if args.type=='multiple':
        multiple_avi(args)
        
if __name__ == '__main__':
    args = gen_parser()    
    
    main(args)