# -*- coding: utf-8 -*-





def author(s,delimiter=' and ',**kwargs):
    num = s.count(delimiter) # count number of authors
    #s = s.replace('{\^e}','ê'.decode('utf-8'))
    if num>1:
        s=s.replace(delimiter,', ',num-1)
        s=s.replace(delimiter,' & ')
    return s
    
def volume(s,typ='console'):
    if typ == 'console':
        return '\033[1m'+s+'\033[0m'
    if typ == 'html':
        return '<b>'+s+'</b>'           

def journal(s,typ='console'):
    #possible type are console, html
    if typ == 'console':
        return '\x1B[3m'+s+'\x1B[23m'
    if typ == 'html':
        return '<i>'+s+'</i>'        
        
def title(s,**kwargs):
    return s
    
def year(s,**kwargs):
    return '('+s+')'
    
def number(s,**kwargs):
    return s
    
def display(d,typ='console'):
    keys = ['author','title','journal','volume','number','year']
    delimiters = [', ',', ','&nbsp', ', ', ' ','']
    functions = [author,title,journal,volume,number,year]
    
    s=u''
    for key,delim,fun in zip(keys,delimiters,functions):
        if key in d.keys():
                s = s + fun(d[key],typ=typ) + delim  
    return s
#        return '<li>'+s+'</li> <br>' #print(s+'<br>')  