# -*- coding: utf-8 -*-
import numpy as np
import glob

import matplotlib.pyplot as plt
import numpy as np
import stephane.tools.rw_data as rw_data

import scipy as scipy
import scipy.special as special

from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from scipy.spatial import SphericalVoronoi
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import pickle
import os


def read_list(filename,Hdelimiter=',',Ddelimiter=','):

    f = open(filename,'r')
    Header,fline=rw_data.read_begin(f,Hdelimiter)
    #Name the fields from the last line of the Header
    Names=rw_data.parseline(Header[-1],Hdelimiter)
    
    Names=[name for name in Names]    
    Data={name:[] for name in Names}  #Just keep the first letter for now (!) should be improved for the first list of characters without space

    print(Names)
    for i,line in enumerate(f):
        List=rw_data.parseline(line,Ddelimiter)
            
        for j,name in enumerate(Names[:-1]):
            Data[name].append(float(List[0]))
        Data[Names[-1]].append([float(l) for l in List[j+1:]])
    f.close()
    return Header,Data
    
    
def read_bubbledata(filelist):
    D={}
    Header,D = rw_data.read_dataFile(filelist[0],Hdelimiter=' ',Ddelimiter=' ')

    for i in range(1,4):
        Header,Data = read_list(filelist[i],Hdelimiter=' ',Ddelimiter=' ')
        for key in Data.keys():
            if not (key in D.keys()):
                D[key]=Data[key]
    print(D.keys())

    return D
    
def process(filename,nmax=6,n=None,nsum=3):    
    Data = {}

    Data['file_interface']=os.path.basename(filename)
    with open(filename, 'rb') as handle:
        data = pickle.load(handle)

    Data['zeta']=data   
    Data['n'] = len(Data['zeta']['p'])
    if n is None or n>Data['n']:
        n = Data['n']
    
    normlist = []
    for p in Data['zeta']['p']:
        norm = np.squeeze(np.linalg.norm(p,axis=1))
        normlist.append(norm)

    Data['zeta']['norm']=normlist
    Alist = Data['zeta']['A']
    #plist = Data['zeta']['p']
    philist = Data['zeta']['phi']
    thetalist = Data['zeta']['theta']

    tc = 2.99
    R0 = 8
    zeta = []
    E = []
    an = []
    time = (Data['zeta']['t'][:n]-Data['zeta']['t'][0])/tc

    Ent = []#np.zeros((nmax,n))
    
#    print(np.sum(Alist[0]))
    for (A,norm,Phi,Theta) in zip(Alist[:n],normlist[:n],philist[:n],thetalist[:n]):
        a = coef_harmonics_fromnorm(norm,A,Theta,Phi,nmax=nmax)
        En = energies(a,nmax=nmax)
    
        an.append(a)

    #avec A de somme unitaire !
        Rmean = np.sum(A*norm)
        zeta.append(np.sqrt(np.sum(A*(norm-Rmean)**2)/R0**2))
        E.append(np.sqrt(4*np.pi*np.sum(En[2:nsum])/R0**2))
        Ent.append(np.sqrt(4*np.pi*En/R0**2))
    
    Ent = np.asarray(Ent)
    
    Data['zeta']['an']=an

    E = np.asarray(E)
    zeta = np.asarray(zeta)
    
    return time,zeta,E,Ent,an
    
    
def PolygonArea(corners):
    #implementation of shoelace formula. Works on convex polygons : Voronoi diagram are convex !
     
    n = len(corners) #  # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

def Area_voronoi(sv,display=False):
    #currently compute the solid angle using the projection on the sphere
    # It neglects the orientation of dA with respect to the radial vector r
    # can't be used to evaluate properly the surface Area
    
    A = []
    for i,region in enumerate(sv.regions):
        u = sv.points[i] # ||u|| = 1 !
    
        #we build the vector normal to the surface
        v1 = np.cross(u,[1,0,0]) # we compute a vector in the plane perpendicular to u
        v2 = np.cross(u,[0,1,0]) # we compute a vector in the plane perpendicular to u
        v3 = np.cross(u,[0,0,1]) # we compute a vector in the plane perpendicular to u
        V = [v1,v2,v3]
        v = V[np.argmax(np.linalg.norm(V))]
        norm = np.sqrt(np.sum(v**2))
        
        v = v/norm # we normalise v !!
        w = np.cross(u,v) # second generator vector in the plane, u & v are unit vectors, so is w
    
        corners = sv.vertices[region]
        c = np.transpose([np.dot(corners,u)]) # coef of the projection of all corners along u

        corners_proj = corners - c*u #we project on the plane (v,w) by removing the u component of each corner    
        corners_proj = np.asarray(np.transpose([np.dot(corners_proj,v),np.dot(corners_proj,w)]))  # projection on the plane (v,w)
        A.append(PolygonArea(corners_proj))

        if display:
            plt.plot(corners_proj[:,0],corners_proj[:,1],'k-')
            plt.plot([corners_proj[0,0],corners_proj[-1,0]],[corners_proj[0,1],corners_proj[-1,1]],'k-')
        
    A = np.asarray(A)
    
    return A
    
def reprojection(sv,f,phi0=np.pi):
    vertices = np.zeros((len(sv.vertices),3))
    for i,vertice in enumerate(sv.vertices):
        re,thetae,phie = cart2sphere(vertice)
        r_estimate = f(thetae,phie+phi0)
        vertices[i,:] = vertice*r_estimate[0]

    sv.vertices = vertices
    return sv

def sphere2cart(r,theta,phi,dim=2):
    if dim==1:
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        return x,y,z
    if dim==2:
        x = r*np.outer(np.sin(theta),np.cos(phi))
        y = r*np.outer(np.sin(theta),np.sin(phi))
        z = r*np.outer(np.cos(theta),np.ones(np.size(phi)))
        return x,y,z
        
    return None
def cart2sphere(U):
    x = U[0]
    y = U[1]
    z = U[2]
    
    r = np.sqrt(x**2+y**2+z**2)
    phi = np.arctan2(y,x)
    theta = np.arccos(z/r)
    
    return r,theta,phi

def compute_voronoi(points):
    sv = scipy.spatial.SphericalVoronoi(points,1,[0,0,0])
    sv.sort_vertices_of_regions()  #important !!! 
    return sv
        
def compute_weights(points):
    sv = scipy.spatial.SphericalVoronoi(points,1,[0,0,0])
    sv.sort_vertices_of_regions()  #important !!! 

    A = Area_voronoi(sv)
    #print(np.sum(A)/4/np.pi) # total area, A ~ 4 pi
    A = A/np.sum(A)*4*np.pi
    return A

def unit_sphere(No=100):
    u = np.linspace(0, 2 * np.pi, No)
    v = np.linspace(0, np.pi, No)
    x = np.outer(np.sin(v),np.cos(u))
    y = np.outer(np.sin(v),np.sin(u))
    z = np.outer(np.cos(v),np.ones(np.size(u)))
    return x,y,z
    
def display_voronoi(sv,ax,spherical=True):
    sv.sort_vertices_of_regions()
    points = sv.points
#    fig = plt.figure(1)
#    fig.set_size_inches(15,7.5)
#    ax1 = fig.add_subplot(121, projection='3d')
#    ax2 = fig.add_subplot(122, projection='3d')
    #ax3 = fig.add_subplot(133, projection='3d')
    if spherical:
        x,y,z = unit_sphere()
        ax.plot_surface(x, y, z, color='y', alpha=0.1)
# plot generator points
#    ax2.scatter(points[:, 0], points[:, 1], points[:, 2], marker='.')
    # plot Voronoi vertices
#    ax2.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], marker='x')
    for region in sv.regions:
        random_color = colors.rgb2hex(np.random.rand(3))
        polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
        polygon.set_color(random_color)
        ax.add_collection3d(polygon)
    
#    return ax1,ax2
#    graphes.legende('$x$','$y$','Voronoi Diagram on the sphere')
#    graphes.save_fig(1,folder+'Sphere_tiling_i0'+str(i0),frmt='png')

def display_voronoi2(sv,f):
    sv.sort_vertices_of_regions()
    fig = plt.figure(1)
    fig.set_size_inches(15,7.5)

    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')
    #ax3 = fig.add_subplot(133, projection='3d')
    
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = -np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))

    ax2.plot_surface(x, y, z, color='y', alpha=0.1)
#    ax3.plot_surface(x, y, z, color='y', alpha=0.1)

    # plot generator points
#    ax2.scatter(points[:, 0], points[:, 1], points[:, 2], marker='.')
    # plot Voronoi vertices
#    ax2.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], marker='x')

    for region in sv.regions:
        random_color = colors.rgb2hex(np.random.rand(3))
        polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
        polygon.set_color(random_color)
        ax2.add_collection3d(polygon)
    
    return ax1,ax2
    
def real_ylm(m,i,Phi,Theta):
    # real harmonics function for surface displacement field

    if m<0:
        y = np.sqrt(2)*(-1)**m*np.imag(special.sph_harm(-m,i,Phi,Theta))
    if m>0:
        y = np.sqrt(2)*(-1)**m*np.real(special.sph_harm(m,i,Phi,Theta))        
    if m==0:
        y = np.real(special.sph_harm(0,i,Phi,Theta)) # il s'agit déjà un réel, juste pour le typage
    return y

def test_shoelace():
    R0 = 2

    for N in range(2,20):
        corners = [[] for i in range(N)]
        for i in range(N):    
            theta = 2*np.pi*i/N
            corners[i] = [R0*np.cos(theta),R0*np.sin(theta)]
            #plt.plot(corners[i][0],corners[i][1],'x')

        Ath = N/2*np.sin(2*np.pi/N)*R0**2
        print(N,PolygonArea(corners),Ath)

# test rms deformation vs mode 2 amplitude
def mesh(No,dim=1,phi0=0):
    epsphi = 2*np.pi/No
    epstheta = np.pi/No
    phi = np.linspace(-phi0,2*np.pi-epsphi-phi0,No)
    theta = np.linspace(epstheta,np.pi-epstheta,No)

    Phi,Theta = np.meshgrid(phi,theta)
    
    if dim==1:
        Phi = np.squeeze(np.reshape(Phi,No**2))
        Theta = np.squeeze(np.reshape(Theta,No**2))
    return Theta,Phi


def shape(a,m,n,Theta=[],Phi=[],R0=[0,0,0],No=30,harm0=1,unitary=False):
    if Theta==[]:
        Theta,Phi = mesh(No)
    #singe harmonic shape
    Y = a*np.real(special.sph_harm(m,n,Phi,Theta))

    (x0,y0,z0) = R0
    x = (harm0+Y)*np.sin(Theta)*np.cos(Phi)+x0
    y = (harm0+Y)*np.sin(Theta)*np.sin(Phi)+y0
    z = (harm0+Y)*np.cos(Theta)+z0

    p = np.transpose([x-x0,y-y0,z-z0])
    #p = np.concatenate((p,[[0,0,1]],[[0,0,-1]]))
#    Theta = np.concatenate((Theta,[0],[np.pi]))
#    Phi = np.concatenate((Phi,[0],[0]))

    norm = np.linalg.norm(p,axis=1,keepdims=True)
    pnorm = p / norm #project on the sphere
    
    sv = scipy.spatial.SphericalVoronoi(pnorm,1,[0,0,0]) #compute the Voronoi
    sv.sort_vertices_of_regions() #sort the mess out #important !!! 

    A = Area_voronoi(sv) #compute the areas
    A = A/np.sum(A)*4*np.pi #normalize solid angles to 1
    
    norm = np.squeeze(norm)
    if unitary:
        R0 = (np.sum(A*norm**3)/4/np.pi)**(1/3) #volume conservation
        p = p/R0 #conserve the volume !!

    return p,sv,A

def shapes(alist,mlist,nlist,Theta=[],Phi=[],R0=[0,0,0],No=30,harm0=1,unitary=False):
    if Theta==[]:
        Theta,Phi = mesh(No)
    #singe harmonic shape
    
    Y = 0*np.real(special.sph_harm(0,0,Phi,Theta))
    for (a,m,n) in zip(alist,mlist,nlist):
        Y += a*np.real(special.sph_harm(m,n,Phi,Theta))

    (x0,y0,z0) = R0
    x = (harm0+Y)*np.sin(Theta)*np.cos(Phi)+x0
    y = (harm0+Y)*np.sin(Theta)*np.sin(Phi)+y0
    z = (harm0+Y)*np.cos(Theta)+z0

    p = np.transpose([x-x0,y-y0,z-z0])
    #p = np.concatenate((p,[[0,0,1]],[[0,0,-1]]))
#    Theta = np.concatenate((Theta,[0],[np.pi]))
#    Phi = np.concatenate((Phi,[0],[0]))

    norm = np.linalg.norm(p,axis=1,keepdims=True)
    pnorm = p / norm #project on the sphere
    
    sv = scipy.spatial.SphericalVoronoi(pnorm,1,[0,0,0]) #compute the Voronoi
    sv.sort_vertices_of_regions() #sort the mess out #important !!! 

    A = Area_voronoi(sv) #compute the areas
    A = A/np.sum(A)*4*np.pi #normalize solid angles to 4 pi
    
    norm = np.squeeze(norm)
    if unitary:
        R0 = (np.sum(A*norm**3)/4/np.pi)**(1/3) #volume conservation
        p = p/R0 #conserve the volume !!

    return p,sv,A,Y
    

def random_points(N):
    Angles = np.reshape(np.random.rand(N*2),(N,2))
    Phi = Angles[:,0]*2*np.pi
    Theta = Angles[:,1]*np.pi

    return (Theta,Phi)

def random_decomposition(nmax=5,a0=1,epsilon=0.1):
    a_n = {}
    
    for i in range(nmax):
        for m in range(-i,i+1):
            a_n[(i,m)] = np.random.rand()*epsilon
    a_n[(0,0)] = a0
    
    return a_n
#a_n[(1,0)] = 1

def shape_from_coef(Theta,Phi,a_n):
    Y_n = []
    for (i,m) in a_n.keys():
        Y_n.append(a_n[(i,m)]*real_ylm(m,i,Phi,Theta))
    
    Y_n = np.asarray(Y_n)
    Y_n = np.sum(Y_n,axis=0)

    x = Y_n*np.sin(Theta)*np.cos(Phi)
    y = Y_n*np.sin(Theta)*np.sin(Phi)
    z = Y_n*np.cos(Theta)

    p = np.transpose([x,y,z])
    norm = np.linalg.norm(p,axis=1,keepdims=True)
    pnorm = p / norm

    return (p,pnorm)

def coef_harmonics(p,A,Theta,Phi,nmax=5):
    a = dict()
    for i in range(nmax):
        for m in range(-i,i+1):
            if i==0:
                Y = np.linalg.norm(p,axis=1)
            else:
                Y = np.linalg.norm(p,axis=1)-a[(0,0)]*real_ylm(0,0,Phi,Theta)

            a[(i,m)] = np.sum(Y*real_ylm(m,i,Phi,Theta)*A)#*np.pi*4/No#
    return a
    
def coef_harmonics_fromnorm(pnorm,A,Theta,Phi,nmax=10):
    a = dict()
    for i in range(nmax):
        for m in range(-i,i+1):
            if i==0:
                print(toto)
                Y = pnorm
            else:
                Y = pnorm-a[(0,0)]*real_ylm(0,0,Phi,Theta)

            a[(i,m)] = np.sum(Y*real_ylm(m,i,Phi,Theta)*A)#*np.pi*4/No#
    return a
    
def energies(a,nmax=5):
    #return the energy in each n mode
    E = np.zeros(nmax)
    for i in range(nmax):
        for m in range(-i,i+1):
            if (i,m) in a.keys():
                E[i] = E[i]+a[(i,m)]**2
            else:
                pass
    return E
    
def compare(a_n,a,nmax=5,display=True):
    errors= dict()
    for i in range(nmax):
        for m in range(-i,i+1):
            if (i,m) in a_n.keys():
                err = (a[(i,m)]-a_n[(i,m)])
                if display:
                    print(str((i,m))+' : '+str(err)+' % '+str(a[(i,m)]))
            else:
                err = a[(i,m)]
                if display:
                    print(str((i,m))+' : '+str(a[(i,m)])+'')
            errors[(i,m)] = err
        
    return errors
    
def test_decomposition(No=900,epsilon=0.1,nmax=5):
    
    (Theta,Phi) = mesh(int(np.sqrt(No)))
#    (Theta,Phi) = random_points(No)

    a_n = random_decomposition(nmax=nmax,a0=1,epsilon=epsilon)

    (p,pnorm) = shape_from_coef(Theta,Phi,a_n)

    sv = scipy.spatial.SphericalVoronoi(pnorm,1,[0,0,0])
    sv.sort_vertices_of_regions()  #important !!! 

    A = Area_voronoi(sv)
    A = A/np.sum(A)*4*np.pi #correct solid angle computation
    #print(np.sum(A)/4/np.pi) # total area, A ~ 4 pi

    a = coef_harmonics(p,A,Theta,Phi,nmax=5)
    
    err = compare(a_n,a,nmax=5,display=False)

    return err
    
def tests_decomposition(N,No=1000,epsilon=0.1,a0=1,nmax=5):
    errors = dict()
    for i in range(N):
        if np.mod(i+1,10)==0:
            print(i+1)
            err = test_decomposition(No=No,epsilon=epsilon,nmax=nmax)
            
            for (i,m) in err.keys():
                if (i,m) in errors.keys():              
                    errors[(i,m)].append(err[(i,m)])
                else:
                    errors[(i,m)] = [err[(i,m)]]
    for (i,m) in errors.keys():
        err = errors[(i,m)]
        print('Error '+str((i,m)) +': '+str(np.mean(err))+' ± '+str(np.std(err)))

    return errors

def test_orthonormal(No=30):
    #test of orthononormal properties
    Theta,Phi = mesh(No)
    p,sv,A = shape(0,0,0,Theta=[],Phi=[],R0=[0,0,0],No=No,harm0=1,unitary=False)

    P = []
    for i in range(5):
        for m in range(-i,i+1):
            S = np.sum(real_ylm(m,i,Phi,Theta)*np.conjugate(real_ylm(m,i,Phi,Theta))*A)
            print(i,m,np.real(S))
            P.append(np.real(S))

    print(np.mean(P),np.std(P))
        #print(np.imag(S))
        
        
def toto_decomposition(No=30):
    Theta,Phi = mesh(No)

    nlist = []
    mlist = []
    alist = []
    avalues = {}
    for n in range(6):
        for m in range(-n,n+1):
            nlist.append(n)
            mlist.append(m)
            avalues[(n,m)]= 0.1*np.random.rand(1)
            alist.append(avalues[(n,m)])
            
    N = len(nlist)
    alist= 0.1*np.random.rand(N)
    
    p,sv,A,Y = shapes(alist,mlist,nlist,Theta=[],Phi=[],R0=[0,0,0],No=30,harm0=1,unitary=False)

    a = dict()
    for i in range(5):
        for m in range(-i,i+1):
            if i==0:
                Y = np.linalg.norm(p,axis=1)
            else:
                Y = np.linalg.norm(p,axis=1)-a[(0,0)]*real_ylm(0,0,Phi,Theta)
    
            a[(i,m)] = np.sum(Y*real_ylm(m,i,Phi,Theta)*A)#*np.pi*4/No#

            print((i,m),a[(i,m)],avalues[(i,m)])
