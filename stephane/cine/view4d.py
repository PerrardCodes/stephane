#!/usr/bin/env python
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from numpy import *
import sys, os
import time
from cine import Cine, Tiff
import argparse
import Image

ESCAPE = '\033'

#XR = 0
#YR = 0
#ZR = 0
    
R = eye(3)

FRAME_RATE = 10
DISPLAY_3D = False
EYE_SPLIT = .02
CURRENT_FRAME = 0.
BRIGHTNESS = 0.25/sqrt(2.)

PLAYING = False
LAST_TIME = None
FULLSCREEN = None

MOUSE_CAPTURED = False
MOUSE_X = 0
MOUSE_Y = 0
FOV = 45
Z_SHIFT = 1.5
Z_RATE = 0.25
DRAW_BOX = True

AUTOROTATE = False 
AUTOROTATE_SPEED = 2 * pi / 5 #0.2 Hz

Z_STEP = 1.

ROTATE_AMOUNT = 10 * pi / 180

SCREENSHOT = None
BACKGROUND_COLOR = 0.1

def display_rot(x = 0, y = 0, z = 0):
    global R
    R = rot_xyz(x, y, z, R)

def if_int(str):
    if not str: return None
    else: return int(str)
        
def main():
    global DATA, WINDOW, FRAMES, UNLOADED_FRAMES, LOAD_START, BASENAME

    glutInit(sys.argv)
    
    parser = argparse.ArgumentParser(description='View a 4D image.')
    parser.add_argument('input', metavar='input', type=str, help='input file (cine or tiff)')
    parser.add_argument('-r', dest='range', type=str, default=":", help='range of frames to display, in python slice format')
    parser.add_argument('-d', dest='depth', type=int, default='240', help='stacks per volume')
    parser.add_argument('-b', dest='brightness', type=float, default=10., help='brightness (multiply image)')
    parser.add_argument('-c', dest='clip', type=float, default=1E-2, help='brigthness to clip to 0')
    parser.add_argument('-s', dest='series', type=int, default=1, help='number of series (channels) in the image, assumed to be alternating frames')
    parser.add_argument('-o', dest='offset', type=int, default=0, help='offset of series to display [0, series-1]')
    parser.add_argument('-D', dest='displayframes', type=str, default="5:-30", help='range of z frames to display in volume')
    parser.add_argument('-x', dest='shear_x', type=float, default=0.0, help='x shear correction')
    parser.add_argument('-y', dest='shear_y', type=float, default=0.0, help='y shear correction')
    
    args = parser.parse_args() 
    #print args.input
    #print args.depth
    #print range(*slice(*[int(x) for x in args.range.split(':')]).indices(100))
    #sys.exit()    window = GLUTWindow('4D Viewer', 1000, 1000)

    #c = Cine('threefold_60k_ratio1_resave.cine')
    #ACTUAL_Z = 103

    display_range = slice(*[if_int(x) for x in args.displayframes.split(':')])
    DATA = Image4D(args.input, args.depth, brighten=args.brightness, clip=args.clip, series=args.series, offset=args.offset, display_range=display_range)
    BASENAME = os.path.splitext(args.input)[0]
    
    FRAMES = range(*slice(*[if_int(x) for x in args.range.split(':')]).indices(len(DATA)))
    UNLOADED_FRAMES = FRAMES[::-1]
    
    WINDOW = GLUTWindow('4D viewer')
    WINDOW.draw_func = draw_scene
    WINDOW.idle_func = draw_scene
    WINDOW.keyboard_func = key_pressed
    WINDOW.special_func = special_pressed
    WINDOW.mouse_func = mouse_func
    WINDOW.motion_func = motion_func
    WINDOW.mouse_wheel_func = mouse_wheel_func
     
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
'''

    print "Loading %d frames..." % (len(FRAMES)),
    sys.stdout.flush()
    LOAD_START = time.time()
    
    WINDOW.start()
    
    

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
    now = time.time()
    print message + '%.3f' % (now - TICK_TIME)
    TICK_TIME = now

class Image4D(object):
    def __init__(self, input, depth, image_size=None, bit_depth=None, brighten=None, aspect_ratio = 1., max_frame=None, clip=None, series=1, offset=0, display_range=slice(None)):
        if type(input) is str:
            if os.path.splitext(input)[1].lower().startswith('.cin'):
                self.source = Cine(input)
                bit_depth = self.source.real_bpp
            elif os.path.splitext(input)[1].lower().startswith('.tif'):
                sys.setrecursionlimit(10**5)
                self.source = Tiff(input)
#                bit_depth = self.source.real_bpp
            else:
                raise RuntimeError('Unsupported file type.')
        else:
            self.source = input

        self.max_frame = (len(self.source)//series)//depth
        if max_frame is not None: self.max_frame = min(self.max_frame, max_frame)
        if brighten == 1: brighten=None
        self.brighten = brighten
        self.frame_buffer = [None] * self.max_frame
        self.depth = depth
        self.aspect_ratio = aspect_ratio
        self.clip = clip
        self.series = series
        self.offset = offset
        self.display_range = list(range(*display_range.indices(self.depth)))
        self.display_depth = len(self.display_range)
        if clip == 0: clip = None
        
        
        test = self.source[0]
        if image_size is None: self.image_size = test.shape
        else: self.image_size = image_size
        
        if bit_depth is None: self.bit_depth = test.itemsize * 8
        else: self.bit_depth = bit_depth

        self.frame_size = (self.display_depth, ) + self.image_size
        
        self.display_size = max(array(self.frame_size) * (aspect_ratio, 1., 1.))
        
        self.gl_texture = {}
            
        
    def get_image(self, image_num):
        img = self.source[image_num * self.series + self.offset]

        if self.brighten or self.clip:
            img = img.astype('f') / (2**self.bit_depth-1) - self.clip
            img = array(clip(img * self.brighten * 255, 0, 255), dtype='u1')
        else:
            if self.bit_depth > 8:
                img = array(img // (2**(self.bit_depth - 8)), dtype='u1')
            elif self.bit_depth < 8:
                img = array(img * 2**(8 - self.bit_depth), dtype='u1')

        return img
    
    def get_frame(self, frame_num, store=False):
        buffer = zeros(self.frame_size)
        
        #print 'Loading frame %d' % frame_num
        i0 = frame_num * self.depth
        for i, j in enumerate(self.display_range):
            buffer[i] = self.get_image(j + i0)
            
        if store: self.frame_buffer[frame_num] = buffer
        return buffer
    
    def __len__(self):
        return self.max_frame
    
    def del_frame(self, frame_num):
        self.frame_buffer[frame_num] = None
        
    def gl_load_texture(self, frame_num):
        if frame_num in self.gl_texture: return True
        
        #tick()
        buffer = self.frame_buffer[frame_num]
        if buffer is None: buffer = self.get_frame(frame_num)
        #tock('Load from disk: ')
        
        #tick()
        self.gl_texture[frame_num] = glGenTextures(1)
        glBindTexture(GL_TEXTURE_3D, self.gl_texture[frame_num])
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
        glTexImage3D(GL_TEXTURE_3D, 0, 1, self.image_size[1], self.image_size[0], self.display_depth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, buffer)
        #tock('Send to opengl: ')

        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)
        
    def bind_frame(self, frame_num):
        if frame_num not in self.gl_texture:
            #print 'Loading %d...' % frame_num
            self.gl_load_texture(frame_num)
            
        glBindTexture(GL_TEXTURE_3D, self.gl_texture[frame_num])
            
            
    def render(self, rot_matrix=eye(3), z_step=1., shear_x=0.0, shear_y=0.0):
        box_size = array(self.frame_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
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
        box_size = array(self.frame_size[::-1], dtype='d') * (1., 1., self.aspect_ratio) / self.display_size
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

        
    def init_GL(self):
        glClearColor(BACKGROUND_COLOR, BACKGROUND_COLOR, BACKGROUND_COLOR, 0.0) 

        if self.depth_test:
            glClearDepth(1.0) 
            glDepthFunc(GL_LESS) 
            glEnable(GL_DEPTH_TEST)

        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_ONE, GL_ONE)
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
        gluPerspective(self.fov, width/float(height), 0.1, 100.0)
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

def draw_scene():
        global LAST_TIME, CURRENT_FRAME, PLAYING, EYE_SPLIT, FRAMES, R, LAST_DRAW, DISPLAY_3D, UNLOADED_FRAMES, AUTOROTATE, AUTOROTATE_SPEED, BASENAME, SCREENSHOT
        
        
        if LAST_TIME is None: LAST_TIME = time.time()

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)      # Clear The Screen And The Depth Buffer

        glLoadIdentity()                                        # Reset The View
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

        display_settings = (cf, EYE_SPLIT, map(tuple, R), DISPLAY_3D, BRIGHTNESS, WINDOW.width, WINDOW.height, MOUSE_CAPTURED, FOV, Z_SHIFT, DRAW_BOX, Z_STEP)

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

            DATA.bind_frame(FRAMES[cf])            
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
        
            #tick()
            
            glutSwapBuffers()
            #tock('Buffer swap: ')

            if SCREENSHOT:
                WINDOW.screenshot(SCREENSHOT)    
                SCREENSHOT = None
            
            LAST_DRAW = display_settings
            #print display_settings

        else:
            if UNLOADED_FRAMES:
                f = UNLOADED_FRAMES.pop()
                DATA.gl_load_texture(f)
                print '%d' % f,
                sys.stdout.flush()
                
                if not UNLOADED_FRAMES: print "Done in %.1fs." % (time.time() - LOAD_START)
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

def key_pressed(*args):
    global BRIGHTNESS, FRAME_RATE, PLAYING, DISPLAY_3D, EYE_SPLIT, FULLSCREEN, \
           WINDOW, DRAW_BOX, Z_SHIFT, R, AUTOROTATE, AUTOROTATE_SPEED, Z_STEP, \
           SCREENSHOT
    
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
    elif args[0] == '1':
        fn = find_new_file(time.strftime('%Y_%m_%d') + '_screenshot_%%04d_%s.png' % BASENAME)
        print "Saving screenshot to " + fn
        SCREENSHOT = fn
        
        #WINDOW.screenshot(fn)
    
    
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
            w = WINDOW.height * tan(0.5 * FOV*pi/180) / Z_SHIFT * 1.45 #WTF is the 1.45?  I don't really know...
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
        