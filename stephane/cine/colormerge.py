#!/usr/bin/env python
from numpy import *
from sparse4d import Sparse4D
from scipy import ndimage
import argparse
import os, sys


def if_int(str):
    if not str: return None
    else: return int(str)

def eval_slice(s, N):
    if s and ':' not in s:
        return (int(s), )
    return range(*slice(*[if_int(x) for x in s.split(':')]).indices(N))
        
parser = argparse.ArgumentParser(description='Convert a 4D image to S4D')
parser.add_argument('input', metavar='input', type=str, nargs=2, help='input files (s4d); second volume is used to color first')
parser.add_argument('-r', dest='range', type=str, default=":", help='range of frames to convert, in python slice format [:]')
parser.add_argument('-b', dest='sigma', type=float, default=5.0, help='blur of coloring (second) file [5]')
parser.add_argument('-s', dest='sat', type=float, default=300., help='saturation multiplier [300]')
args = parser.parse_args()

input1 = Sparse4D(args.input[0])
input2 = Sparse4D(args.input[1])

header = input1.header
header.update({'input1':args.input[0], 'input2':args.input[1]})

def common_sub(a, b):
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]: return a[:i]
    return a[:i+1]
    
output_fn = common_sub(*args.input).rstrip('_-') + '-' + os.path.splitext(sys.argv[0])[0] + '.s4d'
#print args.input, output_fn
#sys.exit()

output = Sparse4D(output_fn, 'w', header)

for i in eval_slice(args.range, len(input1)):
    f1 = input1[i]
    f2 = input2[i].astype('f') / 255.
    f2 = clip(ndimage.gaussian_filter(f2, 5) * args.sat, 0, 1)
    
    frame = zeros(f1.shape + (4,), dtype='u1')
    
    
    frame[..., 0] = f1
    frame[..., 1] = (f1 * (1 - f2)).astype('u1')  #G
    frame[..., 2] = f1
    frame[..., 3] = f1
    
    output.append_array(frame)
    print i    
    
output.close()  