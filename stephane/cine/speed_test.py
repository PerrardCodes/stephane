from numpy import *
import time
from cine import Tiff, Cine
import sys

sys.setrecursionlimit(10**5)
#ls *.cinesource = Tiff('threefold_sub_103FPV_rm.tif')
source = Cine('kf6_60k_250fps_216perV_ratio1.cine')

def tick():
    global TICK_TIME
    TICK_TIME = time.time()
    
def tock(message=''):
    global TICK_TIME
    now = time.time()
    print message + '%.3f' % (now - TICK_TIME)
    TICK_TIME = now
    
test = source[0]
print test.shape
print test.dtype
a = zeros((256,) + test.shape, dtype=test.dtype)

tick()
for n in range(len(a)):
    a[n] = source[n+256]
tock('Load frame: ')

a = a.astype('u4')
tock('Upconvert to floats: ')

a = (a - 2)* 10
tock('Rescale: ')

a = clip(a, 0, 255)
tock('Clip: ')

a = a.astype('u1')
tock('Drop to 8 bits: ')