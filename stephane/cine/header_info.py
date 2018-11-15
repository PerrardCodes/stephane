#!/usr/bin/env python
import sparse4d
import sys

for fn in sys.argv[1:]:
    try:
        i = sparse4d.Sparse4D(fn)
        print "--- %s ---" % fn
        
        for key, val in i.header.iteritems():
            print '%20s: %s' % (key, val)
            
        print
    
    except:
        print "--- Couldn't open '%s' ---" % fn
    
    
