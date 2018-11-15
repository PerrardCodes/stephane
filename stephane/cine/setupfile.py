#!/usr/bin/env python
import numpy as N
import datetime
import sys, glob, os

basic_namespace = {'datetime':datetime.datetime}
for name in ('sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2', 'pi'): basic_namespace[name] = getattr(N, name)

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
            
    
if __name__ == '__main__':
    #fn, vars = get_setup('*.3dsetup', filter_vars={'volume_size':(384, 384, 384), 'filename':'whatever.s4d', 'header':{}}, verbose=1)
    #fn, vars = get_setup('*.3dsetup', filter_vars={'volume_size':(256, 512, 256), 'filename':'whatever.s4d', 'header':{}}, verbose=1)
    fn, vars = get_setup('*.3dsetup', filter_vars={'volume_size':(256, 256, 256), 'filename':'whatever.s4d', 'header':{}}, verbose=1)
    
    for k in sorted(vars.keys()):
        print '%10s = %s' % (k, repr(vars[k]))