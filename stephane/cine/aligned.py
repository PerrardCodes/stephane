#!/usr/bin/env python
from numpy import *

try:
    from scipy import weave
except:
    USE_C_CODE = False
else:
    USE_C_CODE = True
    
import time

UNIT_BOX = array([(x, y, z) for z in arange(2) for y in arange(2) for x in arange(2)], dtype='f')
UNIT_BOX_C = UNIT_BOX - 0.5

BOX_EDGES = array([
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7),
    (0, 2),
    (1, 3),
    (4, 6),
    (5, 7),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7)
], dtype='u4')

BOX_EDGE_NORMALS = (UNIT_BOX_C[BOX_EDGES[:, 0]] + UNIT_BOX_C[BOX_EDGES[:, 1]]) / sqrt(2.)

BOX_FACES = array([
    (1, 3, 7, 5),
    (0, 4, 6, 2),
    (2, 6, 7, 3),
    (0, 4, 5, 1),
    (4, 5, 7, 6),
    (0, 2, 3, 1)
], dtype='u4')

BOX_NORMALS = array([
    ( 1, 0, 0),
    (-1, 0, 0),
    ( 0, 1, 0),
    ( 0,-1, 0),
    ( 0, 0, 1),
    ( 0, 0,-1)
], dtype='f')


class ViewAlignedBox(object):
    def __init__(self, edges=(1., 1., 1.), texture_edges=(1., 1., 1.), texture_offsets=(0., 0., 0.)):
        self.edges = asarray(edges, dtype='f')
        self.texture_edges = asarray(texture_edges, dtype='f')
        self.texture_offsets = asarray(texture_offsets, dtype='f')
        
        #self.texture_extent = asarray(texture_extent, dtype='f')
        
    def calc_planes(self, rot=eye(3), dz=0.2):
        disp_X = dot((UNIT_BOX_C) * self.edges, rot)
        disp_T = UNIT_BOX * (self.texture_edges - self.texture_offsets) + self.texture_offsets
        
        edge_N = dot(BOX_EDGE_NORMALS, rot)
#        angles = nan_to_num(arctan2(edge_N[:, 0], edge_N[:, 1])) #nan_to_num is really slow (17ms!?) and anyway arctan2 handles (0, 0) gracefully
        angles = arctan2(edge_N[:, 0], edge_N[:, 1])
        edge_order = argsort(angles)

        edges = BOX_EDGES[edge_order]
        
        z_max = floor(max(disp_X[:, 2]) / dz) * dz
        Z = arange(-z_max, z_max+dz*0.5, dz)
        
        
        if USE_C_CODE:
            N = len(Z)
            #Assume a maximum of 6 edges per slice on average
            V = zeros((N * 7, 3), dtype='f')
            T = zeros((N * 7, 3), dtype='f')
            C = zeros((N * 6, 3), dtype='u4')
            return_len = zeros(2, dtype='u4')
            
            code = r'''
                unsigned int i, j, k, e1, e2, nV=0, nF=0, nc=0, ic;
                float z, z1, z2, x, v, t;
                
                unsigned int c[12];
                
                //Cycle over z planes
                for(i = 0; i < N; i++)
                {
                    z = Z[i];
                    nc = 0;
                    
                    //Center vertex
                    ic = nV;
                    nV ++;
                    for(k = 0; k < 3; k++)
                    {
                        V[ic*3+k] += 0;
                        T[ic*3+k] += 0;
                    }                    
                    
                    //Cycle over cube edges
                    for (j = 0; j < 12; j++)
                    {
                        e1 = edges[j*2  ] * 3;
                        e2 = edges[j*2+1] * 3;
                        z1 = disp_X[e1+2];
                        z2 = disp_X[e2+2];
                        
                        x = (z - z1) / (z2 - z1);
                        
                        //If this plane intersects an edge
                        if((x > 0.) && (x <= 1.))
                        {
                            
                            c[nc] = nV;
                            
                            //Add a point
                            for(k = 0; k < 3; k++)
                            {
                                v = disp_X[e1+k] * (1 - x) + disp_X[e2+k] * x;
                                t = disp_T[e1+k] * (1 - x) + disp_T[e2+k] * x;
                                V[nV*3+k] = v;
                                T[nV*3+k] = t;
                                V[ic*3+k] += v;
                                T[ic*3+k] += t;
                            }
                            
                            nc ++;
                            nV ++;
                        }
                    }
                    
                    if (nc > 0) //Should always be true, but make sure...
                    {
                        //Cycle over intersection points, making triangles linked to a center point
                        //The center point is not strictly necessary, but it improves the perspective correction quality.
                        for(k = 0; k < 3; k++)
                        {
                            V[ic*3+k] /= float(nc);
                            T[ic*3+k] /= float(nc);
                        }                         
                        
                        for(j = 0; j < nc; j++)
                        {
                            k = (j + 1) % nc;
                            
                            C[nF*3] = ic;
                            C[nF*3+1] = c[j];
                            C[nF*3+2] = c[k];
                            nF ++;
                        }
                    }
                }
                
                
                return_len[0] = nV;
                return_len[1] = nF;
                
                
            '''
            
#            start = time.time()
            weave.inline(code, ('V', 'T', 'C', 'N', 'return_len', 'disp_X', 'disp_T', 'Z', 'edges'))
#            print time.time() - start
            
            V = V[:return_len[0]]
            T = T[:return_len[0]]
            C = C[:return_len[1]]
            
        else:
            V = []
            T = []
            C = []
            
            for z in Z:
                c = []
                for e1, e2 in edges:
                    z1 = disp_X[e1, 2]
                    z2 = disp_X[e2, 2]
                    
                    x = (z - z1) / (z2 - z1)
                    
                    if x > 0 and x <= 1:
                        c.append(len(V))
                        T.append(disp_T[e1] * (1-x) + disp_T[e2] * x)
                        V.append(disp_X[e1] * (1-x) + disp_X[e2] * x)
                        
                for i in range(1, len(c) - 1):
                    C.append([c[0], c[i], c[i+1]])
                    
            V, T, C = map(asarray, (V, T, C))
            
        return V, T, C
    
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
        
if __name__ == '__main__':
    import mesh, glviewer
    
    test = ViewAlignedBox((0.5, 0.5, 1.0))
    
    R = normalize_basis(random.rand(2, 3))
    
    V, T, C = test.calc_planes(R, dz=0.05)

    glviewer.show([mesh.Mesh(V, C), mesh.Mesh(T, C)])
#    V, T, C = test.calc_planes(dz=0.001)

#    print len(C)
    #print V
    #print T
    #print C