#!/usr/bin/env python
from scipy import ndimage
from pylab import *
import sys, os
import Image

ending = '-cropped'

for img in sys.argv[1:]:
    base, ext = os.path.splitext(img)
    if base.endswith(ending): continue
    
    print img

    img = array(Image.open(img))
    if img.shape[2] == 4:
        print '  (has alpha channel, which will be ignored)'
        img = img[..., :3]
    
    border_color = img[0, 0]
    
    labels, features = ndimage.label((img == border_color).all(-1))

    outer_label = labels[0, 0]
    
    border = (labels == outer_label)
    
    bsx = border.sum(0)
    bx = where(bsx < bsx.max())[0]
    mx, Mx = min(bx), max(bx) + 1

    bsy = border.sum(1)
    by = where(bsy < bsy.max())[0]
    my, My = min(by), max(by) + 1

    w = Mx - mx
    h = My - my
    
    offset = w - h
    mx += offset // 2
    Mx = mx + h
    
    print '  Cropped to (%d, %d), (%d, %d)' % (mx, Mx, my, My)
    
    alpha = 255 - 255 * border.astype('u1')
    #border[min(by):max(by)+1, min(bx):max(bx)+1] -= 0.5
    
    
    
    new_img = zeros((My-my, Mx-mx, 4), dtype='u1')
    new_img[..., 3] = alpha[my:My, mx:Mx]
    new_img[..., :3] = img[my:My, mx:Mx]
    
    #imshow(new_img)
    #show()
    ofn = base + ending + ext
    Image.fromarray(new_img).save(ofn)
    print '  -> %s' % ofn
    
    #sys.exit()
