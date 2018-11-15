#!/usr/bin/env python
from view4d import *

Z_SHIFT = 1.8
ANGLES = list(arange(100) * 2 * pi / 100.)
DISPLAY_3D = True

def draw_scene():
        global LAST_TIME, CURRENT_FRAME, PLAYING, EYE_SPLIT, FRAMES, R, \
               LAST_DRAW, DISPLAY_3D, UNLOADED_FRAMES, AUTOROTATE, \
               AUTOROTATE_SPEED, BASENAME, SCREENSHOT, ANGLES, FRAME_COUNT, FRAME_NAME
        
        
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

        if ANGLES:
            R = rot_y(ANGLES.pop())
        else:
            sys.exit()

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
    
        glutSwapBuffers()
        
        WINDOW.screenshot(FRAME_NAME % FRAME_COUNT)
        FRAME_COUNT += 1
        

glutInit(sys.argv)

parser = argparse.ArgumentParser(description='View a 4D image.')
parser.add_argument('input', metavar='input', type=str, help='input file (cine or tiff)')
parser.add_argument('-r', dest='range', type=str, default=":", help='range of frames to display, in python slice format')
parser.add_argument('-d', dest='depth', type=int, default='103', help='stacks per volume')
parser.add_argument('-b', dest='brightness', type=float, default=1., help='brightness (multiply image)')
parser.add_argument('-c', dest='clip', type=float, default=0., help='brigthness to clip to 0')
parser.add_argument('-s', dest='series', type=int, default=1, help='number of series (channels) in the image, assumed to be alternating frames')
parser.add_argument('-o', dest='offset', type=int, default=0, help='offset of series to display [0, series-1]')
parser.add_argument('-D', dest='displayframes', type=str, default=":", help='range of z frames to display in volume')

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


OUTPUT_DIR = find_new_file(BASENAME + '_video%03d')
os.mkdir(OUTPUT_DIR)
FRAME_NAME = os.path.join(OUTPUT_DIR, '%05d.png')

FRAMES = range(*slice(*[if_int(x) for x in args.range.split(':')]).indices(len(DATA)))
UNLOADED_FRAMES = FRAMES[::-1]

WINDOW = GLUTWindow('4D viewer')
WINDOW.draw_func = draw_scene
WINDOW.idle_func = draw_scene
#WINDOW.keyboard_func = key_pressed
#WINDOW.special_func = special_pressed
#WINDOW.mouse_func = mouse_func
#WINDOW.motion_func = motion_func
#WINDOW.mouse_wheel_func = mouse_wheel_func
FRAME_COUNT = 0

#    print "Loading %d frames..." % (len(FRAMES)),
#    sys.stdout.flush()
#    LOAD_START = time.time()

WINDOW.start()

