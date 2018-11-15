#!/usr/bin/env python
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from numpy import *
import sys, os, shutil
import time
#from cine import Cine, Tiff
from sparse4d import Sparse4D
import argparse
import Image
import traceback

ESCAPE = '\033'

#XR = 0
#YR = 0
#ZR = 0
    
R = eye(3)

FRAME_RATE = 10
DISPLAY_3D = False
EYE_SPLIT = .02
CURRENT_FRAME = 0.
BRIGHTNESS = .5

PLAYING = False
LAST_TIME = None
FULLSCREEN = None

MOUSE_CAPTURED = False
MOUSE_X = 0
MOUSE_Y = 0
FOV = 45
Z_SHIFT = 2
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

def display_rot(x = 0, y = 0, z = 0):
    global R
    R = rot_xyz(x, y, z, R)

def if_int(str):
    if not str: return None
    else: return int(str)
    
def make_movie(source, sequence, export_dir, s=slice(None), window_kwargs={}):
    global DATA, WINDOW, BASENAME, SEQUENCE, FRAMES, EXPORT_DIR
    
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

    DATA = Image4D(source)
    FRAMES = range(*s.indices(len(DATA)))

    WINDOW.start()

def show(source, s=slice(None), window_name='4d viewer', mem_buffer=True, window_kwargs={}):
    global DATA, WINDOW, FRAMES
    
    glutInit(sys.argv)

    
    WINDOW = GLUTWindow(window_name, **window_kwargs)
    WINDOW.draw_func = draw_scene
    WINDOW.idle_func = draw_scene
    WINDOW.keyboard_func = key_pressed
    WINDOW.special_func = special_pressed
    WINDOW.mouse_func = mouse_func
    WINDOW.motion_func = motion_func
#    WINDOW.mouse_wheel_func = mouse_wheel_func #Not supported by standard GLUT ):

    glLightfv(GL_LIGHT0, GL_DIFFUSE, (1, 0, 1, 1))
    glLightfv(GL_LIGHT0, GL_POSITION, (1, 1, 1, 0))

    DATA = Image4D(source, mem_buffer=mem_buffer)
    FRAMES = range(*s.indices(len(DATA)))
     
    print '''----Keys----
wasdqe -> Rotate volume
zx -> Zoom
i -> Rest rotation
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
'''

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
    
    show(args.input, s, mem_buffer=args.mem_buffer)

    
    

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
    return rot_z(az, rot_y(ay, rot_x(ax, x)))


UNIT_BOX = array([(x, y, z) for x in (1, -1) for y in (1, -1) for z in (1, -1)])

TICK_TIME = None

def tick():
    global TICK_TIME
    TICK_TIME = time.time()
    
def tock(message=''):
    global TICK_TIME

    if TICK_TIME is None: return

    now = time.time()
    print message + '%.3f' % (now - TICK_TIME)
    TICK_TIME = now

class Image4D(object):
    def __init__(self, input, aspect_ratio=1., mem_buffer=False):
        if type(input) == str:
            self.source = Sparse4D(input, cache_blocks=True, preload=True)

        elif type(input) == ndarray:
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

        self.max_frame = len(self.source)
#        self.depth = depth
#        self.aspect_ratio = aspect_ratio
#        self.display_range = list(range(*display_range.indices(self.depth)))
#        self.display_depth = len(self.display_range)
        
        
        test = self.source[0]
        self.aspect_ratio = 1.
        self.volume_size = test.shape
        
        self.display_size = max(array(self.volume_size) * (aspect_ratio, 1., 1.))
        
        
        self.current_frame_num = -1
        
        self.mem_buffer = mem_buffer
        if mem_buffer:
            self.gl_texture = glGenTextures(1)
            glBindTexture(GL_TEXTURE_3D, self.gl_texture)
            glTexImage3D(GL_TEXTURE_3D, 0, 1, self.volume_size[2], self.volume_size[1], self.volume_size[0], 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, zeros(self.volume_size, dtype='u1'))
        else:    
            self.gl_texture = {}
            
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

        #print buffer.shape

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
        
        if self.mem_buffer:
#            glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, self.volume_size[2], self.volume_size[1], self.volume_size[0], GL_LUMINANCE, GL_UNSIGNED_BYTE, buffer)
            glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, self.volume_size[2], self.volume_size[1], self.volume_size[0], GL_ALPHA, GL_UNSIGNED_BYTE, buffer)
        else:
#            glTexImage3D(GL_TEXTURE_3D, 0, 1, self.volume_size[2], self.volume_size[1], self.volume_size[0], 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, buffer)
            glTexImage3D(GL_TEXTURE_3D, 0, 1, self.volume_size[2], self.volume_size[1], self.volume_size[0], 0, GL_ALPHA, GL_UNSIGNED_BYTE, buffer)
        
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
        box_size = array(self.volume_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
#        print box_size
        corners = UNIT_BOX * .5 * box_size
        r_corners = dot(corners, rot_matrix)
        
        mx, my, mz = r_corners.min(0)
        Mx, My, Mz = r_corners.max(0)
        
        R = rot_matrix.T
        
        z = arange(mz, Mz, z_step/self.display_size)
        N = len(z)
#        tick()
        XY = array((((mx, my, 0), (Mx, my, 0), (Mx, My, 0), (mx, My, 0))))
        ZZ = zeros((N, 3))
        ZZ[:, 2] = z
        QC = XY + ZZ.reshape((N, 1, 3))#array([XY + (0, 0, zp) for zp in z])
        TC = (dot(QC, rot_matrix.T) + 0.5 * box_size) / box_size
#        tock('Array: ')
        

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_TEXTURE_COORD_ARRAY)

        glVertexPointer(3, GL_DOUBLE, 0, QC)
        glTexCoordPointer(3, GL_DOUBLE, 0, TC)
        
        glDrawArrays(GL_QUADS, 0, len(z) * 4)

        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_TEXTURE_COORD_ARRAY)
        
        #glBegin(GL_QUADS)
        #for qc, tc in zip(QC, TC):
        #    for i in range(4):
        #        glTexCoord3fv(tc[i])
        #        glVertex3fv(qc[i])
        

                
        #glEnd()
#        tock('Draw: ')
        
    def draw_box(self, rot_matrix=eye(3)):
        pm = (-1, 1)
        box_size = array(self.volume_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
#        print box_size
        corners = dot(UNIT_BOX * .5 * box_size, rot_matrix)
        
        glEnable(GL_LINE_SMOOTH)
        glBegin(GL_LINES)
        for n in (0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 1, 3, 4, 6, 5, 7, 0, 4, 1, 5, 2, 6, 3, 7):
            glVertex3fv(corners[n])          
        glEnd()
        glDisable(GL_LINE_SMOOTH)
    
        
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
        glEnable(GL_DEPTH_TEST)

        
    def init_GL(self):
        glClearColor(0.0, 0.0, 0.0, 0.0) 

        if self.depth_test:
            glClearDepth(1.0) 
            glDepthFunc(GL_LESS) 
            glEnable(GL_DEPTH_TEST)

        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
        glEnable(GL_BLEND)
        #glBlendFunc(GL_ONE, GL_ONE)
        glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA)
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
        
    def screenshot(self, fn=None, A=False):
        if A:
            channels = "RGBA"
            glc = GL_RGBA
        else:
            channels = "RGB"
            glc = GL_RGB
        
        img = glReadPixels(0, 0, self.width, self.height, glc, GL_UNSIGNED_BYTE)
        #print len(img), self.width, self.height, len(channels)
        img = frombuffer(img, dtype='u1').reshape((self.height, self.width, len(channels)))
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
    'fov':int
}


def movie_draw_scene():
    #Essentially just a copy of draw_scene, hacked up to hijack the playback
    #This could probably be done smarter, but oh well...
    
    try: #Catch exceptions and exit, as there is no keyboard control    
        global WINDOW, SEQUENCE, CURRENT_FRAME, DISPLAY_3D, Z_SHIFT, R, BRIGHTNESS, Z_STEP, DRAW_BOX, FRAME_COUNT, EXPORT_DIR, FOV
        params = SEQUENCE.pop(0)
        
        if 'r' in params: R = normalize_basis(params['r'])
        if '3d' in params: DISPLAY_3D = params['3d']
        if 'frame' in params: CURRENT_FRAME = params['frame']
        if 'z' in params: Z_SHIFT = params['z']
        if 'box' in params: DRAW_BOX = params['box']
        if 'brightness' in params: BRIGHTNESS = params['brightness']
        if 'z_step' in params: Z_STEP = params['z_step']
        if 'fov' in params: FOV = params['fov']
        
        
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

        
        DATA.bind_frame(FRAMES[int(CURRENT_FRAME)])            
        Rv = R.copy()
        Rv[1:, :] *= -1    
    
        if DISPLAY_3D:
            for a, r, g, b in ((EYE_SPLIT, 1, 0, 0), (-EYE_SPLIT, 0, 1, 1)):
                Rp = rot_y(a, Rv)
                
                if DRAW_BOX:
                    glDisable(GL_TEXTURE_3D)
                    glColor3f(.25*r, .25*g, .25*b)
                    DATA.draw_box(Rp)
                            
                glEnable(GL_TEXTURE_3D)
                glColor3f(r*BRIGHTNESS, g*BRIGHTNESS, b*BRIGHTNESS)
                DATA.render(Rp, z_step=Z_STEP)
        else:
            if DRAW_BOX:
                glDisable(GL_TEXTURE_3D)
                glColor3f(.25, .25, .25)
                DATA.draw_box(Rv)
    
            glEnable(GL_TEXTURE_3D)
            glColor3f(BRIGHTNESS, BRIGHTNESS, BRIGHTNESS)
            DATA.render(Rv, z_step=Z_STEP)            
    
        
        glutSwapBuffers()
    
        fn = os.path.join(EXPORT_DIR, '%08d.png' % FRAME_COUNT)
        WINDOW.screenshot(fn)
        print fn
        FRAME_COUNT += 1
    
    
    except:
        traceback.print_exc()
        sys.exit()

    if not SEQUENCE: sys.exit()

def draw_scene():
        global FOV, LAST_TIME, CURRENT_FRAME, PLAYING, EYE_SPLIT, FRAMES, R, LAST_DRAW, DISPLAY_3D, UNLOADED_FRAMES, AUTOROTATE, AUTOROTATE_SPEED, BASENAME, SCREENSHOT
        
        
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
        
        if PLAYING:
            CURRENT_FRAME = (CURRENT_FRAME + FRAME_RATE * elapsed) % len(FRAMES)

        if AUTOROTATE:
            R = rot_y(AUTOROTATE_SPEED * elapsed, R)

        cf = int(CURRENT_FRAME)

        display_settings = (FOV, cf, EYE_SPLIT, map(tuple, R), DISPLAY_3D, BRIGHTNESS, WINDOW.width, WINDOW.height, MOUSE_CAPTURED, FOV, Z_SHIFT, DRAW_BOX, Z_STEP)

        if display_settings != LAST_DRAW or SCREENSHOT:

            if MOUSE_CAPTURED:
                glDisable(GL_TEXTURE_3D)
                glColor3f(.25, .25, .25)
                glBegin(GL_LINE_LOOP)
                draw_gl_circle(r=.5)
                glEnd()
                glBegin(GL_LINES)
                glVertex3f(.5, 0, 0)
                glVertex3f(-.5, 0, 0)
                glVertex3f(0, .5, 0)
                glVertex3f(0, -.5, 0)
                glEnd()

#            tick()
            DATA.bind_frame(FRAMES[cf])            
#            tock("Texture: ")
            
            Rv = R.copy()
            Rv[1:, :] *= -1
            
#            tick()

            def inner_draw():
                glDisable(GL_TEXTURE_3D)
                glDisable(GL_BLEND)
                glColor4f(1, 1, 1, 1)
                glEnable(GL_LIGHT0)
                glEnable(GL_LIGHTING)
                glutSolidSphere(0.25, 10, 10)
                glDisable(GL_LIGHTING)
                glEnable(GL_BLEND)

                if DRAW_BOX:
                    glDisable(GL_TEXTURE_3D)
                    glColor4f(.25, .25, .25, .25)
                    DATA.draw_box(Rv)


                glEnable(GL_TEXTURE_3D)
                glColor4f(1, 1, 1, BRIGHTNESS)
                DATA.render(Rv, z_step=Z_STEP) 

            if DISPLAY_3D:
                for a, r, g, b in ((EYE_SPLIT, 1, 0, 0), (-EYE_SPLIT, 0, 1, 1)):
                    Rp = rot_y(a, Rv)
                    glColorMask(r, g, b, True)
                    glClear(GL_DEPTH_BUFFER_BIT)
                    
                    glDisable(GL_TEXTURE_3D)
                    glDisable(GL_BLEND)
                    glColor3f(1, 1, 1)
                    glEnable(GL_LIGHT0)
                    glEnable(GL_LIGHTING)
                    glutSolidSphere(0.25, 10, 10)
                    glDisable(GL_LIGHTING)
                    glEnable(GL_BLEND)
                    
                    if DRAW_BOX:
                        glDisable(GL_TEXTURE_3D)
                        #glColor3f(.25*r, .25*g, .25*b)
                        glColor3f(.25, .25, .25)
                        DATA.draw_box(Rp)
                                
                    glEnable(GL_TEXTURE_3D)
                    #glColor3f(r*BRIGHTNESS, g*BRIGHTNESS, b*BRIGHTNESS)
                    glColor3f(BRIGHTNESS, BRIGHTNESS, BRIGHTNESS)
                    DATA.render(Rp, z_step=Z_STEP)
                    
                glColorMask(True, True, True, True)
            else:
                inner_draw()           
        
#            tock('Frame time: ')
#            tick()
            
            glutSwapBuffers()

            if SCREENSHOT:
                WINDOW.screenshot(SCREENSHOT)    
                SCREENSHOT = None
            
            LAST_DRAW = display_settings
            #print display_settings

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

def key_pressed(*args):
    global BRIGHTNESS, FRAME_RATE, PLAYING, DISPLAY_3D, EYE_SPLIT, FULLSCREEN, \
           WINDOW, DRAW_BOX, Z_SHIFT, R, AUTOROTATE, AUTOROTATE_SPEED, Z_STEP, \
           SCREENSHOT, MOVIE_FILE, CURRENT_FRAME, LAST_MOVIE_FRAME, FOV
    
    if args[0] == ESCAPE:
        sys.exit()       
    elif args[0] in ('-', '_'):
        BRIGHTNESS /= 2**.5
    elif args[0] in ('+', '='):
        BRIGHTNESS *= 2**.5
        if BRIGHTNESS > 1: BRIGHTNESS = 1.
    elif args[0] in ('<', ','):
        FRAME_RATE /= 2**.5
    elif args[0] in ('>', '.'):
        FRAME_RATE *= 2**.5
    elif args[0] in (' '):
        PLAYING = not PLAYING
    elif args[0] in ('3'):
        DISPLAY_3D = not DISPLAY_3D
    elif args[0] in ('[', '{'):
        EYE_SPLIT /= 2**.25
    elif args[0] in (']', '}'):
        EYE_SPLIT *= 2**.25
    elif args[0] in ('\\', '|'):
        EYE_SPLIT *= -1
    elif args[0] == 'a':
#        YR -= ROTATE_AMOUNT
        display_rot(y = -ROTATE_AMOUNT)
    elif args[0] == 'd':
#        YR += ROTATE_AMOUNT
        display_rot(y = ROTATE_AMOUNT)
    elif args[0] == 'w':
#        XR += ROTATE_AMOUNT
        display_rot(x = ROTATE_AMOUNT)
    elif args[0] == 's':
#        XR -= ROTATE_AMOUNT
        display_rot(x = -ROTATE_AMOUNT)
    elif args[0] == 'q':
#        XR -= ROTATE_AMOUNT
        display_rot(z = ROTATE_AMOUNT)
    elif args[0] == 'e':
#        XR -= ROTATE_AMOUNT
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
    elif args[0] == 'f':
        print CURRENT_FRAME
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
            
    elif args[0] == '1':
        fn = find_new_file('screenshot_%%04d_%s.png' % BASENAME)
        print "Saving screenshot to " + fn
        SCREENSHOT = fn
        #WINDOW.screenshot(fn)

    elif args[0] == 'm':
        cf = int(CURRENT_FRAME)
        if MOVIE_FILE is None:
            MOVIE_FILE = find_new_file('%s_movie_%%04d.txt' % BASENAME)
            print "Created movie file: %s" % MOVIE_FILE
            open(MOVIE_FILE, 'w').write('''#Autogenerated on: %s
source: %s

single:
    R: %s
    frame: %d
    z: %f
    brightness: %f
    3d: %s
    box: %s
    z_step: %f
    fov: %f
    
''' % (time.strftime("%Y %b %d %a %H:%M:%S %Z"), DATA.source.fn, format_array(R[:2, :]), int(CURRENT_FRAME), Z_SHIFT, BRIGHTNESS, DISPLAY_3D, DRAW_BOX, Z_STEP, FOV))

        else:
            if LAST_MOVIE_FRAME == cf: frames = 10
            else: frames = abs(cf - LAST_MOVIE_FRAME)
            open(MOVIE_FILE, 'a').write('''steps: %d
    R: %s
    frame: %d
    z: %f

''' % (frames, format_array(R[:2, :]), cf, Z_SHIFT))
        
        LAST_MOVIE_FRAME = cf    
    
def mouse_func(button, state, x, y):
    global MOUSE_CAPTURED, MOUSE_X, MOUSE_Y
    
    
#    print button, state, x, y
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
        if FOV:
            w = WINDOW.height * tan(0.5 * FOV*pi/180) / Z_SHIFT * 1.45 #WTF is the 1.45?  I don't really know...
        else:
            w = WINDOW.height / float(Z_SHIFT) *0.5
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
        