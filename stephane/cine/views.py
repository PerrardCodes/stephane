#!/usr/bin/env python
import sys, os, subprocess, shutil, glob

steps = 300

movie_file = '''
source: %s

single:
    box: True
    fov: 0
    R: %s
    brightness: 2.0
    frame: 0
    z: 3

steps: %d
    frame: %d
'''

views = {
    'front':'(( 1.000, 0.000, 0.000), ( 0.000, 1.000, 0.000))',
    'side': '(( 0.000, 0.000, 1.000), ( 0.000, 1.000, 0.000))',
    'top':  '(( 1.000, 0.000, 0.000), ( 0.000, 0.000, 1.000))',
}

for input in sys.argv[1:]:
    basename = os.path.splitext(input)[0]
    tmp_dir = 'temp-' + os.path.split(basename)[-1] + '-' + os.path.split(sys.argv[0])[-1]
    
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    
    for name, R in views.iteritems():
        mname = os.path.join(tmp_dir, name)
        f = open(mname + '.txt', 'w')
        f.write(movie_file % (input, R, steps, steps))
        f.close()
        
        subprocess.Popen(['python', 'movie.py', mname+'.txt']).wait()
        subprocess.Popen(['img2avi', '-q95', '-o', basename + '-' + name + '.avi'] + glob.glob(os.path.join(mname, '*.tga'))).wait()
        
        
    shutil.rmtree(tmp_dir)