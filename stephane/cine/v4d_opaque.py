#!/usr/bin/env python
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.arrays import vbo
from OpenGL.GL import shaders
from OpenGL.GL.ARB.framebuffer_object import *
from OpenGL.GL.ARB.multitexture import *

try: from OpenGL.GLUT.freeglut import *
except: print "No freeglut!"

from OpenGL.GLU import *
from numpy import *
import sys, os, shutil, glob
import time
#from cine import Cine, Tiff
#from sparse4d import Sparse4D
from cine import Sparse4D
import argparse
import Image
import traceback
#import setupfile


import datetime




#try:
#    import aligned
#except:
#    USE_ALIGNED = False
#    print "WARINING: failed to import aligned module, display may be slow!"
#else:
#    USE_ALIGNED = True
    
#USE_ALIGNED = False

ESCAPE = '\033'

#XR = 0
#YR = 0
#ZR = 0
    
R = eye(3)


OPACITY = 0.1
#TEXTURE_SHEAR = 0.15
#TEXTURE_PERSPECTIVE_SCALE = 0.08
FRAME_RATE = 10
DISPLAY_3D = False
EYE_SPLIT = 0.04
CURRENT_FRAME = 0.
BRIGHTNESS = 1

PLAYING = False
LAST_TIME = None
FULLSCREEN = None

MOUSE_CAPTURED = False
MOUSE_X = 0
MOUSE_Y = 0
FOV = 45
Z_SHIFT = 4.
Z_RATE = 0.25
DRAW_BOX = True

AUTOROTATE = False
AUTOROTATE_SPEED = 2 * pi / 5 #0.2 Hz

Z_STEP = 1.

ROTATE_AMOUNT = 10 * pi / 180

SCREENSHOT = None

BASENAME = os.path.splitext(sys.argv[0])[0]

SEQUENCE = None #Used for movies
FRAME_COUNT = 0

COLOR_SCALE = 0
COLOR_CYCLE = False

EDGES = array([[-1, 1.]]*3)


BACKGROUND_COLOR = 0.05
BACKGROUND_COLORS = (0.0, 0.05, 0.25, 1.0)

GREEN_MAGENTA = False

PERSPECTIVE_CORRECT = True

def display_rot(x = 0, y = 0, z = 0):
    global R
    R = rot_xyz(x, y, z, R)

def if_int(str):
    if not str: return None
    else: return int(str)
    
from scipy.interpolate import interp1d
h_rgb_i = interp1d(arange(7), (
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 1, 1),
    (0, 0, 1),
    (1, 0, 1),
    (1, 0, 0),
), axis=0)

h_rgb = lambda h: h_rgb_i(h%6)
    
    
basic_namespace = {'datetime':datetime.datetime}
for name in ('sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2', 'pi'): basic_namespace[name] = globals()[name]#getattr(N, name)

def get_line(f):
    while True:
        line = f.readline()
        if not line: return None
        line = line.strip('\n\r')
        if '#' in line: line = line[:line.index('#')]
        line = line.strip()
        if line: return line
        

def eval_str_commands(l):
    ns = basic_namespace.copy()
    dummy = eval('True', ns) #import __builtins__, etc.
    import_vars = ns.keys()   

    for k, v in l: ns[k] = eval(v, ns)
    for var in import_vars: del ns[var]
    
    return ns
    

def get_setup(filter='*.setup', search_dirs=[''], filter_vars={}, default_namespace=basic_namespace, first_line='**setup**', verbose=False, default={}, skip_filters=False, get_string=False):
    default_namespace = default_namespace.copy()
    dummy = eval('True', default_namespace) #import __builtins__, etc.
    import_vars = default_namespace.keys()
    
    f = None
    #print search_dirs
    str_commands = filter_vars.items()
    
    for dir in search_dirs:
        for fn in sorted(glob.glob(os.path.join(dir, filter))):
            if f is not None: f.close()
            f = open(fn)
            
            line = get_line(f)
            if first_line and not line.startswith(first_line): continue
            else: line = get_line(f)
            
            if verbose >= 1: print "--%s--" % fn
            
            while line.startswith('filter:'):
                line = line[len('filter:'):].strip()
                if not skip_filters:
                    try:
                        filter = eval(line, filter_vars)
                        if verbose >= 1: print "   '%s' -> %s" % (line, bool(filter))
                        if not filter:
                            break
                    except:
                        if verbose >= 1: print "   Couldn't evaluate line '%s'" % line
                        break
                    
                line = get_line(f)
                
            else: #Passed filters!
                if verbose >= 1: print "   Passed all filters, evaluating!"
                namespace = basic_namespace.copy()
                namespace.update(filter_vars)
                
                try:
                    while line is not None:
                        var, statement = map(str.strip, line.split('=', 1))
                        if verbose >= 2:
                            print 'Evaluating: "%s" -> "%s"' % (var, statement)
                        namespace[var] = eval(statement, namespace)
                        if verbose >= 1:
                            print '%s = %s' % (var, namespace[var])
                        str_commands.append((var, statement))
                        
                        line = get_line(f)
                except:
                    print "Error evaulating '%s' in setup file '%s' -- aborting setup that passed filters!" % (line, fn)
                    continue
                
                
                for var in import_vars: del namespace[var]
                
                if get_string:
                    return (fn, namespace, str_commands)
                else:
                    return (fn, namespace)
                    
    else:
        if verbose >= 1: print "No valid setup files found, returning default."
        if get_string:
            return (None, default, str_commands)
        else:
            return (None, default)

    
def make_movie(source, sequence, export_dir, s=slice(None), window_kwargs={}):
    global DATA, WINDOW, BASENAME, SEQUENCE, FRAMES, EXPORT_DIR, EDGES
    
    SEQUENCE = sequence
    
    print 'Exporting movie to %s' % export_dir
    if not os.path.exists(export_dir):
        os.mkdir(export_dir)
    else:
        print "Export directory already exists!  Overwrite? [y/N] (All files will be deleted!)"
        answer = raw_input()
        if answer.lower().startswith('y'):
            shutil.rmtree(export_dir)
            os.mkdir(export_dir)
        else:
            sys.exit()
            
    EXPORT_DIR = export_dir 

    glutInit(sys.argv)

    
    WINDOW = GLUTWindow('V4D Movie Generator', **window_kwargs)
    WINDOW.draw_func = movie_draw_scene
    WINDOW.idle_func = movie_draw_scene
    WINDOW.keyboard_func = movie_key_pressed

    DATA = Image4D(source)
#    if 'shear' in DATA.header:
#        TEXTURE_SHEAR = 0
#    if 'perspective' in DATA.header:
#        TEXTURE_PERSPECTIVE_SCALE = 0
    EDGES *= array([[DATA.shader_vars['size_x'], DATA.shader_vars['size_y'], DATA.shader_vars['size_z']]]).T

    FRAMES = range(*s.indices(len(DATA)))

    WINDOW.start()

def show(source, s=slice(None), window_name='4d viewer', mem_buffer=True, window_kwargs={}):
    global DATA, WINDOW, FRAMES, EDGES #, TEXTURE_SHEAR, TEXTURE_PERSPECTIVE_SCALE 
    
    glutInit(sys.argv)

    
    WINDOW = GLUTWindow(window_name, **window_kwargs)
    WINDOW.draw_func = draw_scene
    WINDOW.idle_func = draw_scene
    WINDOW.keyboard_func = key_pressed
    WINDOW.special_func = special_pressed
    WINDOW.mouse_func = mouse_func
    WINDOW.motion_func = motion_func
    WINDOW.mouse_wheel_func = mouse_wheel_func #Not supported by standard GLUT ):

    DATA = Image4D(source, mem_buffer=mem_buffer)
    EDGES *= array([[DATA.shader_vars['size_x'], DATA.shader_vars['size_y'], DATA.shader_vars['size_z']]]).T

    #if 'shear' in DATA.header:
    #    TEXTURE_SHEAR = 0.0
    #if 'perspective' in DATA.header:
    #    TEXTURE_PERSPECTIVE_SCALE = 0.0
    
    FRAMES = range(*s.indices(len(DATA)))
     
    print '''----Keys----
wasdqe -> Rotate volume
zx -> Zoom
i -> Reset rotation
o -> Go to first frame
arrows -> Skip forward/backward
space -> Pause/play
+- -> Adjust brightness
<> -> Adjust playback speed
ESC -> exit
3 -> Activate/deactivate stereo anaglyph
[] -> Adjust eye distance
\ -> Invert eye position
Tab -> Toggle fullscreen
Left mouse/Drag -> Rotate
r -> Autorotate
ty -> Adjust autorotate speed
jkl -> Coarsen / reset / make finer the z stack (default is "correct").
1 -> Take screenshot
p -> Toggle perspective
m -> Create / add entry to movie file
f -> Print current frame
2 -> Print frame display time
c -> Cycle color coding
v -> Toggle cyclic color
gh -> Adjust FOV
QAWSEDRFTGYH -> Adjust display extents
u -> Reset display extents
b -> Toggle outline
n -> Toggle background color
'''
#90 -> Adjust x shear
#78 -> Adjust perspective correction

    WINDOW.start()    
        
def main():
    global BASENAME

    
    parser = argparse.ArgumentParser(description='View a 4D image.')
    parser.add_argument('input', metavar='input', type=str, help='input file (s4d)')
    parser.add_argument('-r', dest='range', type=str, default=":", help='range of frames to display, in python slice format')
    parser.add_argument('-m', dest='mem_buffer', type=bool, nargs='?', default=True, const=False, help='Buffer all frames in OpenGL instead of unpacking each on demand -- not recommended unless you have an excess of video/system memory, and even then its only marginally faster.')
    
    args = parser.parse_args()
    
#    print args.mem_buffer
#    sys.exit()

    BASENAME = os.path.splitext(args.input)[0]
    
    s = slice(*[if_int(x) for x in args.range.split(':')])
    
    show(args.input, s, mem_buffer=args.mem_buffer, window_name=BASENAME)

    
    

def rot_x(a, x=eye(3)):
    x = array(x)
    rx = x.copy()
    rx[..., 1] = cos(a) * x[..., 1] - sin(a) * x[..., 2]
    rx[..., 2] = cos(a) * x[..., 2] + sin(a) * x[..., 1]
    return rx
    
def rot_y(a, x=eye(3)):
    x = array(x)
    rx = x.copy()
    rx[..., 2] = cos(a) * x[..., 2] - sin(a) * x[..., 0]
    rx[..., 0] = cos(a) * x[..., 0] + sin(a) * x[..., 2]
    return rx
 
def rot_z(a, x=eye(3)):
    x = array(x)
    rx = x.copy()
    rx[..., 0] = cos(a) * x[..., 0] - sin(a) * x[..., 1]
    rx[..., 1] = cos(a) * x[..., 1] + sin(a) * x[..., 0]
    return rx

def rot_xyz(ax, ay, az, x=eye(3)):
    return gram_schmidt(rot_z(az, rot_y(ay, rot_x(ax, x))))


#UNIT_BOX = array([(x, y, z) for x in (1, -1) for y in (1, -1) for z in (1, -1)])
#BOX_FACES = array([
#    (1, 3, 7, 5),
#    (0, 4, 6, 2),
#    (2, 6, 7, 3),
#    (0, 4, 5, 1),
#    (4, 5, 7, 6),
#    (0, 2, 3, 1)
#], dtype='u4')

UNIT_BOX = array([(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)], dtype='f')

BOX_FACES = array([
    (0, 1, 3, 2), #-X
    (4, 6, 7, 5), #+X
    (0, 4, 5, 1), #-Y
    (2, 3, 7, 6), #+Y
    (0, 2, 6, 4), #-Z
    (1, 5, 7, 3), #+Z
], dtype='u4')

BOX_EDGES = [(i, j) for i in range(8) for j in range(8) if ((i < j) and sum((UNIT_BOX[i] - UNIT_BOX[j])**2) == 1.)]

#BOX_EDGES = array([
#    (0, 1), (1, 3), (3, 2), (2, 0),
#    (4, 6), (6, 7), (7, 5), (5, 4),
#    ()
#])



glUniforms = {
    (int,   1):glUniform1i,
    (int,   2):glUniform2i,
    (int,   3):glUniform3i,
    (int,   4):glUniform4i,
    (float, 1):glUniform1f,
    (float, 2):glUniform2f,
    (float, 3):glUniform3f,
    (float, 4):glUniform4f,
}

def set_uniforms(program, **vars):
    
    for key, val in vars.iteritems():
        val = asarray(val)
        if not val.shape: val.shape = (1,)
        
        if len(val) > 4: raise ValueError('at most 4 values can be used for set_uniforms')
        if val.dtype in ('u1', 'u2', 'u4', 'u8', 'i1', 'i2', 'i4', 'i8'): dt = int
        elif val.dtype in ('f', 'd'): dt = float
        else: raise ValueError('values for set_uniforms should be ints or floats')

        #print key, dt, len(val)
        glUniforms[dt, len(val)](glGetUniformLocation(program, key), *val)

#modification of OpenGL.GL.shaders.compileProgram, which does not allow for setting variables before validation.
def compileProgram(*s, **vars):
    program = glCreateProgram()
    
    for ss in s:
        glAttachShader(program, ss)
        
    glLinkProgram(program)
        
    shaders.glUseProgram(program)
    set_uniforms(program, **vars)
    shaders.glUseProgram(0)
        
    glValidateProgram(program)
    
    validation = glGetProgramiv(program, GL_VALIDATE_STATUS)
    if validation == GL_FALSE:
        raise RuntimeError(
            """Validation failure (%s): %s"""%(
            validation,
            glGetProgramInfoLog(program),
        ))
    
    link_status = glGetProgramiv(program, GL_LINK_STATUS)
    if link_status == GL_FALSE:
        raise RuntimeError(
            """Link failure (%s): %s"""%(
            link_status,
            glGetProgramInfoLog(program),
        ))
    
    for ss in s:
        glDeleteShader(ss)
        
    return shaders.ShaderProgram(program)

TICK_TIME = None

def tick():
    global TICK_TIME
    TICK_TIME = time.time()
    
def tock(message=''):
    global TICK_TIME

    if TICK_TIME is None:
        TICK_TIME = time.time()
        return None

    now = time.time()
    print message + '%.3f' % (now - TICK_TIME)
    TICK_TIME = now
    
I3x = lambda x, y, z: x
I3y = lambda x, y, z: y
I3z = lambda x, y, z: z


def eval_vals(s):
    setup = {}
    dummy = {}
    exec("from math import *", dummy)
    exec(s, dummy, setup)
    return setup

def try_eval(s):
    while type(s) == str:
        try:
            s = ast.literal_eval(s)
        except:
            return s
    return s

class Image4D(object):
    setup_vars = ['header', 'filename', 'volume_size', 'color_channels', 'max_frame']

    def __init__(self, input, aspect_ratio=1., mem_buffer=False, load_setup=True, setup_dirs=['']):
        self.color_channels = 1
        self.downsample = 1
        
        if type(input) == str:
            self.filename = input
            self.source = Sparse4D(input, cache_blocks=True, preload=True)
            self.header = self.source.header
            tf = self.source[0]
            if tf.ndim == 4:
                self.color_channels = tf.shape[3]
                

        elif type(input) == ndarray:
            self.filename = None
            self.header = {}
            if len(input.shape) == 3: #Accept 3D arrays, but convert to 4D
                input.shape = (1,) + input.shape
            
            if input.dtype == 'u1':            
                self.source = input
            elif input.dtype in ('d', 'f'):
                self.source = (clip(input, 0, 1) * 255).astype('u1')
            else:
                raise ValueError('Input array type must be unsigned bytes or float/double (assumed 0-1 brightness scaling)')


        else:
            raise ValueError('Input type should be filename (s4d file) or array.')

        if self.filename: setup_dirs.insert(0, os.path.split(self.filename)[0])

        self.max_frame = len(self.source)
#        self.depth = depth
#        self.aspect_ratio = aspect_ratio
#        self.display_range = list(range(*display_range.indices(self.depth)))
#        self.display_depth = len(self.display_range)
        
        test = self.source[0]
        if self.downsample > 1:
            test = test[::self.downsample, ::self.downsample, ::self.downsample]
            
        self.aspect_ratio = 1.
        self.volume_size = test.shape
        #print self.volume_size
        if self.color_channels != 1: self.volume_size = self.volume_size[:-1]
        #print self.volume_size
        
        self.display_size = max(array(self.volume_size) * (aspect_ratio, 1., 1.))
        
        
        self.current_frame_num = -1
        
        if self.color_channels != 1:
            if self.color_channels == 4:
                self.internal_format = GL_RGBA
                self.format = GL_RGBA
            else:
                raise ValueError('3D texture must have 1 (I) or 4 (RGBA) color channels!')
        else:
            self.internal_format = GL_INTENSITY
            self.format = GL_LUMINANCE
        
        
        self.mem_buffer = True#mem_buffer
        if self.mem_buffer:
            self.gl_texture = glGenTextures(1)
            glBindTexture(GL_TEXTURE_3D, self.gl_texture)
#            glTexImage3D(GL_TEXTURE_3D, 0, GL_INTENSITY, self.volume_size[2], self.volume_size[1], self.volume_size[0], 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, zeros(self.volume_size, dtype='u1'))
            glTexImage3D(GL_TEXTURE_3D, 0, self.internal_format, self.volume_size[2], self.volume_size[1], self.volume_size[0], 0, self.format, GL_UNSIGNED_BYTE, test)#zeros(self.volume_size, dtype='u1'))
        else:    
            self.gl_texture = {}
            
    
            
        filter_vars = dict((var, getattr(self, var)) for var in self.setup_vars)
        
        has_3dsetup = try_eval(self.header.get('use_3dsetup_perspective', False))
        
        if has_3dsetup:
            #print self.header['3dsetup']
            self.setup = eval_vals(self.header['3dsetup'])
            #print self.header['3dsetup']
            self.setup_file = '-- embedded in s4d --'
        else:
            self.setup_file, self.setup = get_setup('*.3dsetup',
                    setup_dirs, filter_vars, default={'x_func':I3x, 'y_func':I3y, 'z_func':I3z})
        
        if self.setup_file:
            print "Using 3d setup file: %s" % self.setup_file
            for fn, func in zip(('x_func', 'y_func', 'z_func'), (I3x, I3y, I3z)):
                if fn not in self.setup:
                    print "Warning: '%s' not in setup file -- no correction applied on this axis." % fn
                    self.setup[fn] = func
        
        else:
            print "No filter matching setup file found; no perspective correction applied."
            
        #print self.setup
         
         #pcx_x = 0.986451, pcx_z  = 0.164055, pcx_xz = 0.0392977,
         #        pcy_y = 1.000000, pcy_xy =-0.006814, pcy_yz = 0.0409707,
         #        pcz_z = 1.000000, pcz_xz = 0.057537,
         #        size_x = 1.0, size_y = 1.5, size_z = 1.0, step_size=1./128   
            
        self.shader_vars = {
            'pcx_x' : float(self.setup.get('x_x',  1)),
            'pcx_y' : float(self.setup.get('x_y',  0)),
            'pcx_z' : float(self.setup.get('x_z',  0)),
            'pcx_xy': float(self.setup.get('x_xy', 0)),
            'pcx_xz': float(self.setup.get('x_xz', 0)),
            'pcy_y' : float(self.setup.get('y_y',  1)),
            'pcy_x' : float(self.setup.get('y_x',  0)),
            'pcy_z' : float(self.setup.get('y_z',  0)),
            'pcy_xy': float(self.setup.get('y_xy', 0)),
            'pcy_yz': float(self.setup.get('y_yz', 0)),
            'pcz_z' : float(self.setup.get('z_z',  1)),
            'pcz_x' : float(self.setup.get('z_x',  0)),
            'pcz_y' : float(self.setup.get('z_y',  0)),
            'pcz_xz': float(self.setup.get('z_xz', 0)),
            'pcz_yz': float(self.setup.get('z_yz', 0)),
            'step_size': float(2./self.volume_size[2]),
            'size_x': 1., #X is the reference size
            'size_y': float(self.volume_size[1]) / self.volume_size[2],
            'size_z': float(self.volume_size[0]) / self.volume_size[2],
        }
        
        #print self.shader_vars
        
    def __len__(self):
        return self.max_frame
    
    def gl_load_texture(self, frame_num):
        #tick()
        if self.mem_buffer:
            glBindTexture(GL_TEXTURE_3D, self.gl_texture)
        else:
            if frame_num in self.gl_texture: return True
            self.gl_texture[frame_num] = glGenTextures(1)
            glBindTexture(GL_TEXTURE_3D, self.gl_texture[frame_num])
        
        buffer = self.source[frame_num]

        if self.downsample > 1:
            buffer = buffer[::self.downsample, ::self.downsample, ::self.downsample]

        #print buffer.shape

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
        
        if self.mem_buffer:
#            glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, self.volume_size[2], self.volume_size[1], self.volume_size[0], GL_LUMINANCE, GL_UNSIGNED_BYTE, buffer)
            glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, self.volume_size[2], self.volume_size[1], self.volume_size[0], self.format, GL_UNSIGNED_BYTE, buffer)
        else:
            glTexImage3D(GL_TEXTURE_3D, 0, GL_INTENSITY, self.volume_size[2], self.volume_size[1], self.volume_size[0], 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, buffer)
        
        #tock('Send to opengl: ')

        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)
        
    def bind_frame(self, frame_num):
        if frame_num == self.current_frame_num:
            return True
        
        if self.mem_buffer:
            self.gl_load_texture(frame_num)
            tex_id = self.gl_texture
            
            
        else:
            if frame_num not in self.gl_texture:
                #print 'Loading %d...' % frame_num
                self.gl_load_texture(frame_num)
            tex_id = self.gl_texture[frame_num]
            
        glBindTexture(GL_TEXTURE_3D, tex_id)
        self.current_frame_num = frame_num
            
            
    def render(self, rot_matrix=eye(3), z_step=1.):
        global BRIGHTNESS, COLOR_SCALE, COLOR_CYCLE, EDGES #TEXTURE_SHEAR, TEXTURE_PERSPECTIVE_SCALE, 
        
        box_size = array(self.volume_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
        
        box_size *= EDGES[:, 1] - EDGES[:, 0]
        

        if USE_ALIGNED:
            V, T, C = aligned.ViewAlignedBox(box_size, EDGES[:, 1], EDGES[:, 0]).calc_planes(rot_matrix, z_step / self.display_size)
            
            #zp = (T[..., 2] - 0.5) * self.volume_size[0] / self.volume_size[2]

            if COLOR_SCALE:
                if COLOR_SCALE < 4:
                    h = T[:, COLOR_SCALE - 1].copy()
                else:
                    h = -V[:, 2].copy()
                    h -= h.min()
                    h /= h.max()
                    
                if COLOR_CYCLE: h *= 18
                else: h *= 4
                
                Color = ones((len(h), 4), dtype='f') * BRIGHTNESS
                Color[:, :3] *= h_rgb(h) * .8 + .2
                #Color = hstack((Color, ones((len(Color), 1), dtype='f')))
            
            #T[..., 0] += TEXTURE_SHEAR * zp
            
            #T[..., :2] = (T[..., :2] - 0.5) * (1 + zp[..., newaxis] * TEXTURE_PERSPECTIVE_SCALE) + 0.5
            
            x0, y0, z0 = ((T - 0.5) * 2.).T
            
            ys = self.volume_size[1] / float(self.volume_size[2])
            zs = self.volume_size[0] / float(self.volume_size[2])
            #print ys, zs
            
            y0 *= ys
            z0 *= zs
            
            
            T[:, 0] = self.setup['x_func'](x0, y0, z0) * 0.5 + 0.5
            T[:, 1] = self.setup['y_func'](x0, y0, z0)/ys * 0.5 + 0.5
            T[:, 2] = self.setup['z_func'](x0, y0, z0)/zs * 0.5 + 0.5
            
            if COLOR_SCALE:
                glEnableClientState(GL_COLOR_ARRAY)
                glColorPointer(4, GL_FLOAT, 0, Color)
            
            glEnableClientState(GL_VERTEX_ARRAY)    
            glVertexPointer(3, GL_FLOAT, 0, V)

            glEnableClientState(GL_TEXTURE_COORD_ARRAY)
            glTexCoordPointer(3, GL_FLOAT, 0, T)

            glDrawElements(GL_TRIANGLES, 3 * len(C), GL_UNSIGNED_INT, C)
            
            glDisableClientState(GL_VERTEX_ARRAY)
            glDisableClientState(GL_TEXTURE_COORD_ARRAY)
            glDisableClientState(GL_COLOR_ARRAY)
        else:
    #        print box_size
            corners = UNIT_BOX * .5 * box_size
            r_corners = dot(corners, rot_matrix)
            
            mx, my, mz = r_corners.min(0)
            Mx, My, Mz = r_corners.max(0)
            
            R = rot_matrix.T
            
            z = arange(mz, Mz, z_step/self.display_size)
            N = len(z)
            XY = array((((mx, my, 0), (Mx, my, 0), (Mx, My, 0), (mx, My, 0))))
            ZZ = zeros((N, 3))
            ZZ[:, 2] = z
            QC = XY + ZZ.reshape((N, 1, 3))#array([XY + (0, 0, zp) for zp in z])
            TC = (dot(QC, rot_matrix.T) + 0.5 * box_size) / box_size
                
    
            glEnableClientState(GL_VERTEX_ARRAY)
            glEnableClientState(GL_TEXTURE_COORD_ARRAY)
    
            glVertexPointer(3, GL_DOUBLE, 0, QC)
            glTexCoordPointer(3, GL_DOUBLE, 0, TC)
            
            glDrawArrays(GL_QUADS, 0, len(z) * 4)
    
            glDisableClientState(GL_VERTEX_ARRAY)
            glDisableClientState(GL_TEXTURE_COORD_ARRAY)
        
    def draw_box(self, rot_matrix=eye(3)):
#        pm = (-1, 1)
        box_size = array(self.volume_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
        box_size *= EDGES[:, 1] - EDGES[:, 0]
        corners = dot(UNIT_BOX * .5 * box_size, rot_matrix)
        
        glEnable(GL_LINE_SMOOTH)
        glBegin(GL_LINES)
        for n in (0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 1, 3, 4, 6, 5, 7, 0, 4, 1, 5, 2, 6, 3, 7):
            glVertex3fv(corners[n])          
        glEnd()
        glDisable(GL_LINE_SMOOTH)
        
    def draw_background(self, rot_matrix=eye(3)):
        box_size = array(self.volume_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
        box_size *= EDGES[:, 1] - EDGES[:, 0]
        corners = dot(UNIT_BOX * .5 * box_size, rot_matrix)
        
        glBegin(GL_QUADS)
        for i in BOX_FACES.flat:
            glVertex3fv(corners[i])
        glEnd()
#        print 'yo'
        
class FuncSetter(object):
    def __init__(self, prop_name, glut_func=None):
        self.prop_name = prop_name
        self.glut_func = glut_func
        
    def __get__(self, obj, obj_type):
        return getattr(obj, '__' + self.prop_name, None)
        
    def __set__(self, obj, val):
        try:
            if self.glut_func is not None: self.glut_func(val)
            else: print 'Warning: %s property is ignored.' % self.prop_name
        except:
            print "Failed to set callback for %s" % self.prop_name
        
        
try:
    MWF = glutMouseWheelFunc
except:
    print "Warning: mouse scroll wheel disabled."
    MWF = None
        


class GLUTWindow(object):
    draw_func = FuncSetter('draw_func', glutDisplayFunc)
    idle_func = FuncSetter('idle_func', glutIdleFunc)
    keyboard_func = FuncSetter('keyboard_func', glutKeyboardFunc)
    special_func = FuncSetter('special_func', glutSpecialFunc)
    motion_func = FuncSetter('motion_func', glutMotionFunc)
    passive_motion_func = FuncSetter('passive_motion_func', glutPassiveMotionFunc)
    mouse_func = FuncSetter('mouse_func', glutMouseFunc)
    mouse_wheel_func = FuncSetter('mouse_wheel_func', MWF)
      
    
    def __init__(self, window_name, width=1000, height=1000, fov=45, depth_test=False, min_z=0.01, max_z=10.):
        self.width = width
        self.height = height
        self.fov = fov
        self.depth_test = depth_test
        self.min_z = min_z
        self.max_z = max_z

#        glutInit(sys.argv)

        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
        
        glutInitWindowSize(width, height)
        self.window = glutCreateWindow(window_name)

        glutReshapeFunc(self.resize_GL)

        self.fbo = glGenFramebuffers(1)
        
        try:            
            self.depth_shader_frag = shaders.compileShader("""
                #extension GL_ARB_texture_rectangle : enable
                                                                 
                uniform sampler2DRect back_buffer;
                
                void main() {
                    gl_FragColor = vec4(gl_TexCoord[0].xyz, gl_FragCoord.z);
                }""", GL_FRAGMENT_SHADER)
            self.depth_shader = compileProgram(self.depth_shader_frag)

            self.inteference_shader_frag = shaders.compileShader("""
                #extension GL_ARB_texture_rectangle : enable
                                                                 
                uniform sampler2DRect back_buffer;
                
                void main() {
                    //float x = length(gl_Color.xyz - texture2DRect(back_buffer, gl_FragCoord.st).rgb);
                    float x = (gl_FragCoord.z - texture2DRect(back_buffer, gl_FragCoord.st).a) * 50.0;
                    vec3 int_color = sin(vec3(x * 6.15, x*7.55, x*8.51));
                    int_color = int_color * int_color;
                    gl_FragColor = vec4(int_color, 1.);
                    //gl_FragColor = texture2DRect(back_buffer, gl_FragCoord.st);
                    //gl_FragColor = vec4(0.01 * gl_FragCoord.s, 0.01 * gl_FragCoord.t, 1., 1.);
                }""", GL_FRAGMENT_SHADER)
            
            self.interference_shader = compileProgram(self.inteference_shader_frag)

            self.volume_shader_frag = shaders.compileShader("""
                #extension GL_ARB_texture_rectangle : enable
                                                                 
                uniform sampler2DRect back_buffer;
                uniform sampler3D volume;
                
                uniform float brightness;
                uniform float opacity;
                
                uniform float pcx_x, pcx_y, pcx_z, pcx_xy, pcx_xz, pcx_0;
                uniform float pcy_y, pcy_x, pcy_z, pcy_xy, pcy_yz, pcy_0;
                uniform float pcz_z, pcz_x, pcz_y, pcz_yz, pcz_xz, pcz_0;
                
                uniform float size_x, size_y, size_z;
                uniform float step_size;
                
                void main() {
                    vec4 c;
                    float a;
                    vec3 p1 = texture2DRect(back_buffer, gl_FragCoord.st).xyz;
                    vec3 p0 = gl_TexCoord[0].xyz;
                    vec4 intensity = vec4(0.0);
                    
                    vec3 d = p1 - p0;
                    float l = length(d);
                    
                    vec3 s = d / l * step_size;
                    vec3 p = p0 + 0.5 * s;
                    
                    //Transform into (-1, 1) coordinates
                    //!! Now assumed to be done in original call
                    //p = p * 2.0 - 1.0;
                    //s = s * 2.0;
                    
                    vec3 pp;
                    
                    vec3 tex_rescale = vec3(1./size_x, 1./size_y, 1./size_z);
                    
                    float max_length = 2.1 * sqrt(size_x*size_x + size_y*size_y + size_z*size_z);
                    //Overestimate, just to be safe.
                    
                    if (l < max_length) {
                        for (float ll = 0.5 * step_size; ll < l; ll += step_size) {
                            //c = texture3D(volume, p);
                            pp = vec3(
                                pcx_0  +  pcx_x * p.x  +  pcx_y * p.y  +  pcx_z * p.z  +  pcx_xz * p.x*p.z  +  pcx_xy * p.x*p.y,
                                pcy_0  +  pcy_y * p.y  +  pcy_x * p.x  +  pcy_z * p.z  +  pcy_xy * p.x*p.y  +  pcy_yz * p.y*p.z,
                                pcz_0  +  pcz_z * p.z  +  pcz_x * p.x  +  pcz_y * p.y  +  pcz_xz * p.x*p.z  +  pcz_yz * p.y*p.z
                            );
                            
                            c = texture3D(volume, (pp * tex_rescale + 1.0)*0.5);
                            
                            intensity += clamp(c*c*opacity, 0.0, 1.0) * (1.0-intensity.a); 
                            p += s;
                        }
                    }
            
                    //gl_FragColor = brightness * intensity;
                    gl_FragColor = sqrt(brightness * intensity);
                }""", GL_FRAGMENT_SHADER)            
            
            
            self.volume_shader = compileProgram(self.volume_shader_frag,
                    back_buffer=1, volume=0, brightness=2. / OPACITY, opacity=OPACITY)
            
        except Exception, err:
            print "Shader compile error: "
            for l in err: print l
            #QtCore.QEventLoop.exit()
            sys.exit(1)
            
            

            
        
    def init_GL(self):
#        glClearColor(0.0, 0.0, 0.0, 0.0) 
        glClearColor(BACKGROUND_COLOR, BACKGROUND_COLOR, BACKGROUND_COLOR, 0.0) 
#        glClearColor(1.0, 1.0, 1.0, 0.0) 

        if self.depth_test:
            glClearDepth(1.0) 
            glDepthFunc(GL_LESS) 
            glEnable(GL_DEPTH_TEST)

        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
        glEnable(GL_BLEND)
        #glBlendFunc(GL_ONE, GL_ONE)
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA)
        glShadeModel(GL_SMOOTH) 
            
        self.resize_GL(self.width, self.height)
        

    def resize_GL(self, width=None, height=None):
        if width is None: width = self.width
        if height is None: height = self.height
        if height == 0: height = 1
        self.width = width
        self.height = height
    
        glViewport(0, 0, width, height) 
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        
        aspect_ratio = width/float(height)
        if self.fov:
            gluPerspective(self.fov, aspect_ratio, 0.1, 100.0)
        else:
            glOrtho(-.5 * aspect_ratio, .5 * aspect_ratio, -.5, .5, 0.1, 100.0)
            
        glMatrixMode(GL_MODELVIEW)
        
        if hasattr(self, 'back_buffer'):
            glDeleteTextures(self.back_buffer)
            del self.back_buffer
            
        
        #side = 256
        
        glActiveTexture(GL_TEXTURE1)
        self.back_buffer = glGenTextures(1)
        glBindTexture(GL_TEXTURE_RECTANGLE, self.back_buffer)
        glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, None)
        glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glActiveTexture(GL_TEXTURE0)


    def fullscreen(self, fullscreen=False):
        if fullscreen:
            self.pre_height = self.height
            self.pre_width = self.width
            glutFullScreen()
        else:
            glutReshapeWindow(self.pre_width, self.pre_height)       

    def start(self):
        self.init_GL()

        glutMainLoop()
        
    def screenshot(self, fn=None, alpha=False):
        if alpha:
            channels = "RGBA"
            glc = GL_RGBA
        else:
            channels = "RGB"
            glc = GL_RGB
        
        glFlush(); glFinish() #Make sure drawing is complete
        
        glReadBuffer(GL_FRONT)
        img = glReadPixels(0, 0, self.width, self.height, glc, GL_UNSIGNED_BYTE)
        #print len(img), self.width, self.height, len(channels)
        img = frombuffer(img, dtype='u1').reshape((self.height, self.width, len(channels)))

        glFlush(); glFinish() #Make sure drawing is complete

        img = Image.fromarray(img[::-1])
        #img = Image.frombuffer(channels, (self.width, self.height), img, "raw", channels, 0, 0)
        if fn:
            img.save(fn)   
        return img
    
    
    

LAST_DRAW = None

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


valid_movie_options = {
    'r':array,
    '3d':bool,
    'frame':float,
    'z':float,
    'box':bool,
    'brightness':float,
    'z_step':float,
    'fov':int,
    'edges':array,
}

def color_cube():
    global EDGES
    
    X1, X0 = EDGES.T
    
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_COLOR_ARRAY)
    
    ub = (UNIT_BOX * (X1 - X0) + X0).astype('f')
    
    glVertexPointer(3, GL_FLOAT, 0, ub - 0.5 * (X0 + X1))
    glColorPointer(3, GL_FLOAT, 0, ub)
    glDrawElements(GL_QUADS, 4*len(BOX_FACES), GL_UNSIGNED_INT, BOX_FACES)
    
    glDisableClientState(GL_VERTEX_ARRAY)
    glDisableClientState(GL_COLOR_ARRAY)

def texture_cube():
    global EDGES
    
    X1, X0 = EDGES.T
    
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_TEXTURE_COORD_ARRAY)
    
    ub = (UNIT_BOX * (X1 - X0) + X0).astype('f')
    
    glVertexPointer(3, GL_FLOAT, 0, ub - 0.5 * (X0 + X1))
    glTexCoordPointer(3, GL_FLOAT, 0, ub)
    glDrawElements(GL_QUADS, 4*len(BOX_FACES), GL_UNSIGNED_INT, BOX_FACES)
    
    glDisableClientState(GL_VERTEX_ARRAY)
    glDisableClientState(GL_TEXTURE_COORD_ARRAY)


def box():
    global EDGES
    
    X1, X0 = EDGES.T
    
    glEnableClientState(GL_VERTEX_ARRAY)
#    glEnableClientState(GL_COLOR_ARRAY)
    
    ub = (UNIT_BOX * (X1 - X0) + X0).astype('f')
    glEnable(GL_LINE_SMOOTH)
    glVertexPointer(3, GL_FLOAT, 0, ub - 0.5 * (X0 + X1))
#    glColorPointer(3, GL_FLOAT, 0, ub)
    glDrawElements(GL_LINES, 2*len(BOX_EDGES), GL_UNSIGNED_INT, BOX_EDGES)
    glDisable(GL_LINE_SMOOTH)
    glDisableClientState(GL_VERTEX_ARRAY)
#    glDisableClientState(GL_COLOR_ARRAY)        

def movie_draw_scene():
    #Essentially just a copy of draw_scene, hacked up to hijack the playback
    #This could probably be done smarter, but oh well...
    
    try: #Catch exceptions and exit, as there is no keyboard control    
        global WINDOW, SEQUENCE, CURRENT_FRAME, DISPLAY_3D, Z_SHIFT, R, BRIGHTNESS, Z_STEP, DRAW_BOX, FRAME_COUNT, EXPORT_DIR, FOV, EDGES
        params = SEQUENCE.pop(0)
        
        if 'r' in params: R = normalize_basis(params['r'])
        if '3d' in params: DISPLAY_3D = params['3d']
        if 'frame' in params: CURRENT_FRAME = params['frame']
        if 'z' in params: Z_SHIFT = params['z']
        if 'box' in params: DRAW_BOX = params['box']
        if 'brightness' in params: BRIGHTNESS = params['brightness']
        if 'z_step' in params: Z_STEP = params['z_step']
        if 'fov' in params: FOV = params['fov']
        if 'edges' in params: EDGES = params['edges']
        
        
        if FOV != WINDOW.fov:
            WINDOW.fov = FOV
            WINDOW.resize_GL()

        

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)      # Clear The Screen And The Depth Buffer

        glLoadIdentity()                                        # Reset The View
        if not FOV:
            glScalef(1.25/Z_SHIFT, 1.25/Z_SHIFT, 1.25/Z_SHIFT)
            glTranslatef(0, 0, -10)
        else:
            glTranslatef(0.0, 0.0, -Z_SHIFT)                      # Move Into The Screen

        

        inner_draw()
        
        #DATA.bind_frame(FRAMES[int(CURRENT_FRAME)])            
        #Rv = R.copy()
        #Rv[1:, :] *= -1    
        #
        #if DISPLAY_3D:
        #    for a, r, g, b in ((EYE_SPLIT, 1, 0, 0), (-EYE_SPLIT, 0, 1, 1)):
        #        Rp = rot_y(a, Rv)
        #        
        #        if DRAW_BOX:
        #            glDisable(GL_TEXTURE_3D)
        #            glColor3f(.25*r, .25*g, .25*b)
        #            DATA.draw_box(Rp)
        #                    
        #        glEnable(GL_TEXTURE_3D)
        #        glColor3f(r*BRIGHTNESS, g*BRIGHTNESS, b*BRIGHTNESS)
        #        DATA.render(Rp, z_step=Z_STEP)
        #else:
        #    if DRAW_BOX:
        #        glDisable(GL_TEXTURE_3D)
        #        glColor3f(.25, .25, .25)
        #        DATA.draw_box(Rv)
        #
        #
        #    glEnable(GL_TEXTURE_3D)
        #    glColor3f(BRIGHTNESS, BRIGHTNESS, BRIGHTNESS)
        #    DATA.render(Rv, z_step=Z_STEP)            
        #
        #
        #glutSwapBuffers()
    
        fn = os.path.join(EXPORT_DIR, '%08d.tga' % FRAME_COUNT)
        WINDOW.screenshot(fn)
        print fn
        FRAME_COUNT += 1
    
    
    except:
        traceback.print_exc()
        sys.exit()

    if not SEQUENCE: sys.exit()
    
PRINT_FRAMERATE = False
FRAME_COUNT_START = None
FRAME_TIMES = []


def innermost_draw(Rv, r=True, g=True, b=True): # (:
    self = WINDOW
    global BRIGHTNESS, Z_STEP, DATA
    
    glPushMatrix()
   
    R4 = eye(4)
    R4[:3, :3] = Rv
    glMultMatrixd(R4)
    
    glDisable(GL_BLEND)
    glEnable(GL_CULL_FACE)
    
    shaders.glUseProgram(self.depth_shader) 
#    colorMask()

    glColorMask(True, True, True, True)

    glBindFramebuffer(GL_FRAMEBUFFER, self.fbo)
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_RECTANGLE, self.back_buffer, 0)


    flip_sides = False
    
    if flip_sides:
        glCullFace(GL_BACK)
    else:
        glCullFace(GL_FRONT)
        
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    texture_cube()
    
    glBindFramebuffer(GL_FRAMEBUFFER, 0)  

    glColorMask(r, g, b, True)

    
    if not flip_sides:
        glCullFace(GL_BACK)
    else:
        glCullFace(GL_FRONT)

        
    shaders.glUseProgram(self.volume_shader)
    set_uniforms(self.volume_shader, brightness=float(BRIGHTNESS), opacity=float(OPACITY), **DATA.shader_vars)
    if not PERSPECTIVE_CORRECT:
        set_uniforms(self.volume_shader,
                 pcx_x = 1.000000, pcx_y = 0., pcx_z = 0., pcx_xy = 0., pcx_xz = 0.,
                 pcy_y = 1.000000, pcy_x = 0., pcy_z = 0., pcy_xy = 0., pcy_yz = 0.,
                 pcz_z = 1.000000, pcz_x = 0., pcz_y = 0., pcz_xz = 0., pcz_yz = 0.)

#    set_uniforms(self.volume_shader, brightness=float(BRIGHTNESS))
                 #pcx_x = 0.986451, pcx_z  = 0.164055, pcx_xz = 0.0392977,
                 #pcy_y = 1.000000, pcy_xy =-0.006814, pcy_yz = 0.0409707,
                 #pcz_z = 1.000000, pcz_xz = 0.057537,
                 #size_x = 1.0, size_y = 1.5, size_z = 1.0, step_size=1./128)
    
    glActiveTexture(GL_TEXTURE1)
    glBindTexture(GL_TEXTURE_RECTANGLE, self.back_buffer)

    glActiveTexture(GL_TEXTURE0)

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glEnable(GL_DEPTH_TEST)
    glDepthFunc(GL_ALWAYS)    
    
    texture_cube()


    shaders.glUseProgram(0)
    
    if DRAW_BOX:
        lb = 0.2

        glEnable(GL_BLEND)



        
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glDepthFunc(GL_LEQUAL)        
        glColor4f(1, 1, 1, lb)
        box()

        glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_ONE)
        glDepthFunc(GL_GREATER)
        glColor4f(lb, lb, lb, 1)
        box()
        
        glDisable(GL_BLEND)        
        glDepthFunc(GL_LESS)
    
    glDisable(GL_DEPTH_TEST)
    
    glPopMatrix()
 
    glColorMask(True, True, True, True)

        
def draw_box(self, rot_matrix=eye(3)):
#        pm = (-1, 1)
    box_size = array(self.volume_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
    box_size *= EDGES[:, 1] - EDGES[:, 0]
    corners = dot(UNIT_BOX * .5 * box_size, rot_matrix)
    
    glEnable(GL_LINE_SMOOTH)
    glBegin(GL_LINES)
    for n in (0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 1, 3, 4, 6, 5, 7, 0, 4, 1, 5, 2, 6, 3, 7):
        glVertex3fv(corners[n])          
    glEnd()
    glDisable(GL_LINE_SMOOTH)

    glDisable(GL_CULL_FACE)

def inner_draw():
    global FRAMES, BRIGHTNESS, Z_STEP, EYE_SPLIT, DISPLAY_3D, DATA, R, CURRENT_FRAME

    
    cf = int(CURRENT_FRAME + 1E-3)

    Rv = R.copy()
    #Rv[1:, :] *= -1 #Invert axis
    Rv[1, :] *= -1
    

#            tick()
    DATA.bind_frame(FRAMES[cf])            
#            tock("Texture: ")

    
#            tick()



    if DISPLAY_3D:
        if GREEN_MAGENTA:
            argb = ((EYE_SPLIT, 0, 1, 0), (-EYE_SPLIT, 1, 0, 1))
        else:
            argb = ((EYE_SPLIT, 1, 0, 0), (-EYE_SPLIT, 0, 1, 1)) #Red Cyan
        
        for a, r, g, b in argb:
            Rp = rot_y(a, Rv)

#            glColorMask(r, g, b, True)
#            glClear(GL_DEPTH_BUFFER_BIT)            

            #glMatrixMode(GL_MODELVIEW)
            #glPushMatrix()
            #glTranslate(-.5*a, 0, 0)

            innermost_draw(Rp, r, g, b)
            #glPopMatrix()
            
        #glColorMask(True, True, True, True)
    else:
        innermost_draw(Rv)

    
    glutSwapBuffers()

def draw_scene():
        global FOV, LAST_TIME, CURRENT_FRAME, PLAYING, EYE_SPLIT, FRAMES, R, \
            LAST_DRAW, DISPLAY_3D, UNLOADED_FRAMES, AUTOROTATE, AUTOROTATE_SPEED, \
            BASENAME, SCREENSHOT, PRINT_FRAMERATE, FRAME_COUNT_START, FRAME_TIMES, \
            COLOR_SCALE, COLOR_CYCLES, EDGES, BACKGROUND_COLOR, GREEN_MAGENTA, PERSPECTIVE_CORRECT, OPACITY
        ##TEXTURE_SHEAR, TEXTURE_PERSPECTIVE_SCALE, 
        
        if LAST_TIME is None: LAST_TIME = time.time()

        if FOV != WINDOW.fov:
            WINDOW.fov = FOV
            WINDOW.resize_GL()

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)      # Clear The Screen And The Depth Buffer

        glLoadIdentity()                                        # Reset The View
        if not FOV:
            glScalef(1.25/Z_SHIFT, 1.25/Z_SHIFT, 1.25/Z_SHIFT)
            glTranslatef(0, 0, -10)
        else:
            glTranslatef(0.0, 0.0, -Z_SHIFT)                      # Move Into The Screen
        #R = rot_xyz(XR, YR, ZR)
        
        current_time = time.time()
        elapsed = current_time - LAST_TIME
        LAST_TIME = current_time
        
        if FRAME_COUNT_START is None: FRAME_COUNT_START = current_time
        
        if PLAYING:
            CURRENT_FRAME = (CURRENT_FRAME + FRAME_RATE * elapsed) % len(FRAMES)

        if AUTOROTATE:
            R = rot_y(AUTOROTATE_SPEED * elapsed, R)
            #R = dot(rot_y(AUTOROTATE_SPEED * elapsed), R)

        display_settings = (FOV, int(CURRENT_FRAME), EYE_SPLIT, map(tuple, R), DISPLAY_3D,
                            BRIGHTNESS, WINDOW.width, WINDOW.height,
                            MOUSE_CAPTURED, FOV, Z_SHIFT, DRAW_BOX, Z_STEP,
                            #TEXTURE_SHEAR, TEXTURE_PERSPECTIVE_SCALE,
                            BACKGROUND_COLOR, GREEN_MAGENTA,
                            COLOR_SCALE, COLOR_CYCLE, LINE_THICKNESS, list(EDGES.flatten()), PERSPECTIVE_CORRECT, OPACITY)

        if display_settings != LAST_DRAW or SCREENSHOT:
            inner_draw()

            if SCREENSHOT:
                WINDOW.screenshot(SCREENSHOT)    
                SCREENSHOT = None
            
            LAST_DRAW = display_settings
            #print display_settings
            
            if PRINT_FRAMERATE:
                FRAME_TIMES.append(time.time() - current_time)
                
                if current_time - FRAME_COUNT_START > 1:
                    print '%5.1f' % (1. / mean(FRAME_TIMES))
                    FRAME_COUNT_START = current_time
                    FRAME_TIMES = []

        else:
            time.sleep(1E-2)
            
def mouse_wheel_func(button, dir, x, y):
    global Z_SHIFT
    Z_SHIFT *= (1 + Z_RATE * dir)

def find_new_file(name_str, max=10000):
    for i in xrange(max):
        if not os.path.exists(name_str % i):
            return name_str % i
    else:
        raise RuntimeError('Failed to find unused filename for "%s"' % name_str)

MOVIE_FILE = None

def format_array(a, precision=3):
    fmt = '%% 1.%sf' % precision
    if len(a.shape) > 1:
        return '(' + ', '.join(format_array(aa) for aa in a) + ')'
    else:
        return '(' + ','.join(fmt % aa for aa in a) + ')'
    

LAST_MOVIE_FRAME = 0
LINE_THICKNESS = 1

EDGE_ADJUST = ('Q', 'A', 'W', 'S', 'E', 'D', 'R', 'F', 'T', 'G', 'Y', 'H')

LAST_MOVIE_SETTINGS = {}

def movie_key_pressed(*args):
    if args[0] == ESCAPE:
        sys.exit()       
    
    
def key_pressed(*args):
    global BRIGHTNESS, FRAME_RATE, PLAYING, DISPLAY_3D, EYE_SPLIT, FULLSCREEN, \
           WINDOW, DRAW_BOX, Z_SHIFT, R, AUTOROTATE, AUTOROTATE_SPEED, Z_STEP, \
           SCREENSHOT, MOVIE_FILE, CURRENT_FRAME, LAST_MOVIE_FRAME, FOV, \
           PRINT_FRAMERATE, USE_ALIGNED, COLOR_SCALE, COLOR_CYCLE, LINE_THICKNESS, \
           EDGES, FOV, BACKGROUND_COLOR, GREEN_MAGENTA, LAST_MOVIE_SETTINGS, \
           PERSPECTIVE_CORRECT, OPACITY
            #TEXTURE_SHEAR, TEXTURE_PERSPECTIVE_SCALE, 
    
    
    if args[0] == ESCAPE:
        sys.exit()       
    elif args[0] == 'L':
        LINE_THICKNESS = (LINE_THICKNESS % 10) + 1
        glLineWidth(LINE_THICKNESS);
    elif args[0] in ('-', '_'):
        OPACITY /= 2**.25
    elif args[0] in ('+', '='):
        OPACITY *= 2**.25
        #if BRIGHTNESS > 1: BRIGHTNESS = 1.
    elif args[0] in ('<', ','):
        FRAME_RATE /= 2**.5
    elif args[0] in ('>', '.'):
        FRAME_RATE *= 2**.5
    elif args[0] in (' '):
        PLAYING = not PLAYING
    elif args[0] in ('3'):
        DISPLAY_3D = not DISPLAY_3D
    elif args[0] in ('[', '{'):
        EYE_SPLIT /= 2**(1./8)
        print "Total eye shear: %.3f" % (EYE_SPLIT*2)
    elif args[0] in (']', '}'):
        EYE_SPLIT *= 2**(1./8)
        print "Total eye shear: %.3f" % (EYE_SPLIT*2)
    elif args[0] in ('\\', '|'):
        EYE_SPLIT *= -1
    elif args[0] == 'g':
        FOV *= sqrt(2.)
        Z_SHIFT /= sqrt(2.)
        print 'FOV = %.2f' % FOV
    elif args[0] == 'h':
        FOV /= sqrt(2.)
        Z_SHIFT *= sqrt(2.)
        print 'FOV = %.2f' % FOV
    elif args[0] == 'a':
        display_rot(y = -ROTATE_AMOUNT)
    elif args[0] == 'd':
        display_rot(y = ROTATE_AMOUNT)
    elif args[0] == 'w':
        display_rot(x = ROTATE_AMOUNT)
    elif args[0] == 's':
        display_rot(x = -ROTATE_AMOUNT)
    elif args[0] == 'q':
        display_rot(z = ROTATE_AMOUNT)
    elif args[0] == 'e':
        display_rot(z = -ROTATE_AMOUNT)
    elif args[0] == '\t':
        FULLSCREEN = not FULLSCREEN
        WINDOW.fullscreen(FULLSCREEN)
    elif args[0] == 'b':
        DRAW_BOX = not DRAW_BOX
    elif args[0] == 'z':
        Z_SHIFT /= 2.**.25
    elif args[0] == 'x':
        Z_SHIFT *= 2.**.25
    elif args[0] == 'i':
        R = eye(3)
    elif args[0] == 'I':
        EDGES = array([[0., 1]]*3)
    elif args[0] == 'o':
        CURRENT_FRAME = 0.;
    elif args[0] == 'f':
        print int(CURRENT_FRAME)
    elif args[0] == 'r':
        AUTOROTATE = not AUTOROTATE
    elif args[0] == 't':
        AUTOROTATE_SPEED /= 2**.25
    elif args[0] == 'y':
        AUTOROTATE_SPEED *= 2**.25
    elif args[0] == 'j':
        Z_STEP *= 2**.5
    elif args[0] == 'k':
        Z_STEP = 1.
    elif args[0] == 'l':
        Z_STEP /= 2**.5
    elif args[0] == 'p':
        if FOV:
            FOV = 0
        else:
            FOV = 45
    elif args[0] == 'P':
        PERSPECTIVE_CORRECT = not PERSPECTIVE_CORRECT
        print 'Perspective correction: %s' % str(PERSPECTIVE_CORRECT)
    elif args[0] == '2':
        PRINT_FRAMERATE = not PRINT_FRAMERATE
            
    elif args[0] == 'c':
        COLOR_SCALE = (COLOR_SCALE + 1) % 5
    elif args[0] == 'v':
        COLOR_CYCLE = not COLOR_CYCLE
            
    elif args[0] == '1':
        fn = find_new_file('%s_frame%03d_%%04d.png' % (BASENAME, int(CURRENT_FRAME)))
        print "Saving screenshot to " + fn
        SCREENSHOT = fn
        #WINDOW.screenshot(fn)

    #elif args[0] == '0':
    #    TEXTURE_SHEAR += 0.01
    #    print "TEXTURE_SHEAR = %0.3f" % TEXTURE_SHEAR
    #elif args[0] == '9':
    #    TEXTURE_SHEAR -= 0.01
    #    print "TEXTURE_SHEAR = %0.3f" % TEXTURE_SHEAR
    #
    #elif args[0] == '8':
    #    TEXTURE_PERSPECTIVE_SCALE += 0.01
    #    print "TEXTURE_PERSPECTIVE_SCALE = %0.3f" % TEXTURE_PERSPECTIVE_SCALE
    #elif args[0] == '7':
    #    TEXTURE_PERSPECTIVE_SCALE -= 0.01
    #    print "TEXTURE_PERSPECTIVE_SCALE = %0.3f" % TEXTURE_PERSPECTIVE_SCALE

    elif args[0] == '4':
        USE_ALIGNED = not USE_ALIGNED
        print 'USE_ALIGNED = %s' % repr(USE_ALIGNED)

    elif args[0] == 'u':
        EDGES = array([[0., 1]]*3)

    elif args[0] == 'm':
        cf = int(CURRENT_FRAME)
        
        movie_settings = {
            'R':format_array(R),
            'frame':cf,
            'z':Z_SHIFT,
            'brightness':BRIGHTNESS,
            '3d':DISPLAY_3D,
            'box':DRAW_BOX,
            'z_step':Z_STEP,
            'fov':FOV,
            'edges':format_array(EDGES)
        }
        
        settings = []
        
        for k, v in movie_settings.iteritems():
            if getattr(LAST_MOVIE_SETTINGS, k, None) != v:
                settings.append('    %s: %s' % (k, v))
        
        LAST_MOVIE_SETTINGS = movie_settings.copy()
       
        if MOVIE_FILE is None:
            MOVIE_FILE = find_new_file('%s_movie_%%04d.txt' % BASENAME)
            print "Created movie file: %s" % MOVIE_FILE
            open(MOVIE_FILE, 'w').write('#Autogenerated on: %s\n\nsource: %s\n\nsingle:\n%s\n' % (time.strftime("%Y %b %d %a %H:%M:%S %Z"), DATA.source.fn, '\n'.join(settings)))
            #, format_array(R[:2, :]), int(CURRENT_FRAME), Z_SHIFT, BRIGHTNESS, DISPLAY_3D, DRAW_BOX, Z_STEP, FOV, format_array(EDGES)))

        else:
            if LAST_MOVIE_FRAME == cf: frames = 10
            else: frames = abs(cf - LAST_MOVIE_FRAME)
            open(MOVIE_FILE, 'a').write('steps: %d\n%s\n' % (frames, '\n'.join(settings)))
        
        LAST_MOVIE_FRAME = cf
        
    elif args[0] in EDGE_ADJUST:
        i = EDGE_ADJUST.index(args[0])
        axis = i // 4
        side = (i%4) // 2
        other_side = 1 - side
        pm = (i%2)*2-1
        
        EDGES[axis, side] += pm * 0.05
        
        print '(' + ', '.join('%.2f' % n for n in tuple(EDGES.flatten())) + ')'
    
    elif args[0] == 'n':
        if BACKGROUND_COLOR in BACKGROUND_COLORS:
            i = BACKGROUND_COLORS.index(BACKGROUND_COLOR)
        else:
            i = 0
        
        BACKGROUND_COLOR = BACKGROUND_COLORS[(i+1)%len(BACKGROUND_COLORS)]
            
        glClearColor(BACKGROUND_COLOR, BACKGROUND_COLOR, BACKGROUND_COLOR, 0.0)
    
    elif args[0] == '#':
        GREEN_MAGENTA = not GREEN_MAGENTA

    
def mouse_func(button, state, x, y):
    global MOUSE_CAPTURED, MOUSE_X, MOUSE_Y
    
    
    #print button, state, x, y
    if button == GLUT_LEFT:
        if state == GLUT_DOWN:
#            print 'captured'
            MOUSE_CAPTURED = True
#            glutSetCursor(GLUT_CURSOR_NONE)
            MOUSE_X = x
            MOUSE_Y = y
        else:
#            print 'released'
            MOUSE_CAPTURED = False
#            glutSetCursor(GLUT_CURSOR_INHERIT)
        
def draw_gl_circle(center=zeros(3), x=(1, 0, 0), y=(0, 1, 0), r=1., np = 100):
    theta = arange(np) * (2*pi) / np
    C = cos(theta)
    S = sin(theta)
    x = array(x)
    y = array(y)
    for c, s in zip(C, S):
        glVertex3fv(center + r * (x * c + y * s))
        
    
def draw_gl_line():
    pass
        
def motion_func(x, y):
    global MOUSE_CAPTURED, MOUSE_X, MOUSE_Y, WINDOW, R, FOV, Z_SHIFT #XR, YR, ZR

    if MOUSE_CAPTURED:
        dx = x - MOUSE_X
        dy = y - MOUSE_Y
        MOUSE_X = x
        MOUSE_Y = y
        
        #print dx, dy
        #if FOV:
        #    w = WINDOW.height * tan(0.5 * FOV*pi/180) / Z_SHIFT * 1.45 #WTF is the 1.45?  I don't really know...
        #else:
        #    w = WINDOW.height / float(Z_SHIFT) *0.5
        w = WINDOW.height /3.
        #print (x, y, w)
        #w = 300.
        dx /= w
        dy /= w
        x = (x - WINDOW.width / 2.) / w
        y = (y - WINDOW.height / 2.) / w

        #print x, y, dx, dy

        #Prevent blowups... probably unnecessary, but doens't hurt
        if abs(x) < 1E-3: x = 1E-3
        if abs(y) < 1E-3: y = 1E-3
        
        r = sqrt(x**2 + y**2)
        #print r
        phi = arctan2(y, x)

        r_hat = array([cos(phi), sin(phi)])
        phi_hat = array([-sin(phi), cos(phi)]) 
        
        dr = dot(r_hat, (dx, dy))
        dphi = dot(phi_hat, (dx, dy))

        if r > 1:
            dphi /= r
            r == 1.
        
        r_xy = r_hat * dr + (1 - r) * dphi * phi_hat
        r_z  = r * dphi

        display_rot(y = r_xy[0], x = r_xy[1], z = -r_z)
        
#            YR += dx / w
#            XR += dy / w

#        print x - MOUSE_X, y - MOUSE_Y
#        glutWarpPointer(MOUSE_X, MOUSE_Y)
        
def special_pressed(*args):
    global CURRENT_FRAME, FRAMES
    
    if args[0] == GLUT_KEY_LEFT:
        CURRENT_FRAME = (CURRENT_FRAME - 1) % len(FRAMES)
    elif args[0] == GLUT_KEY_RIGHT:
        CURRENT_FRAME = (CURRENT_FRAME + 1) % len(FRAMES)
    elif args[0] == GLUT_KEY_UP:
        CURRENT_FRAME = (CURRENT_FRAME + 10) % len(FRAMES)
    elif args[0] == GLUT_KEY_DOWN:
        CURRENT_FRAME = (CURRENT_FRAME - 10) % len(FRAMES)



if __name__ == '__main__':
    main()
