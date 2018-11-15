#!/usr/bin/env python
from numpy import *
import v4d_shader as v4d
from scipy.interpolate import interp1d
import sys, os

def gram_schmidt(V):
    V = list(array(v, dtype='d') for v in V)
    U = []
    
    for v in V:
        v2 = v.copy()
        
        for u in U:
            v -= dot(v2, u) * u
        
        U.append(norm(v))
        
    return array(U)
    
def normalize_basis(V):
    V = gram_schmidt(V[:2])
    
    #Ensures right handed coordinate system
    V = vstack((V, cross(V[0], V[1])))
    
    return V
    
def mag(x):
    return sqrt((x*x).sum(len(x.shape) - 1))
    
def norm(x):
    return x / mag(x)
    
def interp_basis(b1, b2, sub_divides=4, zero_tol=0.01):
    B = [gram_schmidt(b1[:2]), gram_schmidt(b2[:2])]
    
    if sqrt(((B[1] - B[0])**2).sum()) < 1E-6:
        return lambda x: normalize_basis(B[1])
    
    for i in range(sub_divides):
        B2 = []
        
        for j, b in enumerate(B):
            B2.append(b)
            
            if (j + 1) < len(B):
                b2 = B[j+1]
                
                bm = (b + b2) / 2.
                m = mag(bm)
                
                z = cross(b[0], b[1])
                if (m < zero_tol).all():
                    bm[0] += 2 * z * zero_tol
                    bm[1] += 2 * b[0] * zero_tol
                elif m[0] < zero_tol:
                    bm[0] += 2 * z * zero_tol
                elif m[1] < zero_tol:
                    bm[1] += 2 * z * zero_tol
                
                    
                B2.append(gram_schmidt(bm))
                
        B = B2
        
    B = map(normalize_basis, B)
    
    d = [0]
    
    for i in range(len(B) - 1):
        dB = B[i] - B[i+1]
        d.append(d[-1] + sqrt((dB**2).sum()))
        
    interp = interp1d(array(d), array(B), axis=0)
    return lambda x: normalize_basis(interp(x * d[-1]))
    
def strip_comments(s):
    if '#' in s:
        s = s[:s.find('#')]
    return s
    
def unpack_commands(s):
    s = s.replace('\t', ' '*8)
    lines = filter(str.strip, map(strip_comments, s.splitlines()))
    
    commands = []
    while lines:
        line = lines.pop(0).strip()
        if ':' in line:
            cmd, arg = map(str.strip, line.split(':', 1))
        else: cmd, arg = line, ''
    
        kwargs = {}    
        while lines and lines[0].startswith(' '):
            line = lines.pop(0).strip()
            if ':' in line: key, val = map(str.strip, line.split(':', 1))
            else: key, val = line, ''
            kwargs[key.lower()] = val
        
        commands.append((cmd.lower(), arg, kwargs))
    
    return commands

#Use eval_kwargs instead!
#def check_kwargs(d, special=[]):
#    for k in d.keys():
#        if k not in v4d.valid_movie_options and k not in special:
#            raise ValueError('invalid movie frame keyword "%s"' % k)
        
def eval_kwargs(d):
    for k, v in d.iteritems():
        if k in v4d.valid_movie_options:
            d[k] = v4d.valid_movie_options[k](eval(v))
        else:
            raise ValueError('invalid movie frame keyword "%s"' % k)
            
rot_d = {
    'x': v4d.rot_x,
    'y': v4d.rot_y,
    'z': v4d.rot_z
}
            
def find_previous(s, k):
    for step in reversed(s):
        if k in step: return step[k]
    else:
        raise ValueError('movie keyword "%s" must be previously defined to execute a sequence.  (e.g. it should appear in a "single" command before "steps")' % k)
    
if __name__ == '__main__':
    
    
    
    #sequence = [{'r':v4d.rot_y(a, eye(3)[:2]), 'frame':30, 'z':2, '3d':True, 'brightness':.5} for a in 2 * pi * arange(30)/ 30]
    
    
    
    #v4d.make_movie('2011_10_19_square_240fpv_250vps_B.s4d', sequence, 'test_movie')
    
    source = None
    sequence = []
    window_kwargs = {}
    
    for cmd, arg, kwargs in unpack_commands(open(sys.argv[1]).read()):
        if cmd == 'source':
            source = arg
    
        if cmd == 'size':
            w, h = map(int, arg.split(','))
            window_kwargs['width'] = w
            window_kwargs['height'] = h
    
        elif cmd == 'single':
            eval_kwargs(kwargs)
            sequence.append(kwargs)
            
        elif cmd == 'steps':
            if 'spin' in kwargs:
                axes, rotations = kwargs['spin'].split(',')
                axes = axes.lower()
                rotations = float(rotations)

                if axes in rot_d:
                    spin_func = lambda x, R: rot_d[axes](x * 2 * pi * rotations, R)
                else:
                    raise ValueError('spin option should have value: (xyz), number of spins')
                
                del kwargs['spin']
            else: spin_func = None
            
            eval_kwargs(kwargs)

            if spin_func is not None and 'r' not in kwargs:
                kwargs['r'] = find_previous(sequence, 'r')

            N = int(eval(arg))
            steps = [{} for n in range(N)]
            
            for k, v in kwargs.iteritems():
                #print k, type(v)
                if type(v) in (list, tuple): v = array(v)
                
                if k == 'r':
                    r0 = find_previous(sequence, 'r')
                    interp_func = interp_basis(r0, v)
                    for i in range(N):
                        x = (i+1) / float(N)
                        R = interp_func(x)
                        
                        if spin_func is not None:
                            R = spin_func(x, R)
                            
                        steps[i]['r'] = R
                    
                elif type(v) in (float, ndarray) or k == 'frame':
                    v0 = find_previous(sequence, k)
                    for i in range(N):
                        x = (i+1) / float(N)
                        steps[i][k] = x * v + (1-x) * v0  
                        
                else:
                    steps[0][k] = v
            
            sequence += steps
            
#    for s in sequence:
#        print s
    
    v4d.make_movie(source, sequence, os.path.splitext(sys.argv[1])[0], window_kwargs=window_kwargs)
    
    #
    #b1 = eye(3)
    #b2 = eye(3) * (-1, -1, 1)
    ##b2 = rot_y(pi/2.)
    #
    #basis_func = interp_basis(b1, b2)
    #
    #for z in linspace(0, 1, 5):
    #    print basis_func(z)