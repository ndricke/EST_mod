import numpy as np
import random
import math

##Other algorithms for generating a spherical covering
#def sph2cart(r,v,p):
#    x = r*np.sin(v)*np.cos(p)
#    y = r*np.sin(v)*np.sin(p)
#    z = r*np.cos(v)
#    return x,y,z
#
#def distSph1(N,r):
#    count = 0
#    a = 4.0*np.pi*r*r/N #approximate area per point on the sphere
#    d = a**0.5
#    Mv = int(np.pi/d)
#    dv = np.pi/Mv
#    dp = a/dv
#
#    xcoord = []; ycoord = []; zcoord = [] #list of coordinates
#    for i in range(Mv):
#        v = np.pi*(i+0.5)/Mv
#        Mp = int(2.0*np.pi*np.sin(v/dp))
#        for j in range(Mp):
#            p = 2.0*np.pi*j/Mp
#            x,y,z = sph2cart(r,v,p)
#            xcoord.append(x);ycoord.append(y);zcoord.append(z)
#            count += 1
#    print "Total numer of points:",count
#    return xcoord,ycoord,zcoord

#Although random, this one seems to work fine for doing a pretty good even covering
#It doesn't seem to have problematic fluctuations in covering
def fiboSph(samples=1,radius=1.):#,randomize=True):

##I think this would make some of the coordinates a little random, but it doesn't seem necessary
    rnd = 1.
    radius2 = pow(radius,2)
#    if randomize:
#        rnd = random.random() * samples

    xcoord = []; ycoord = []; zcoord = [] #list of coordinates
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = radius*(((i * offset) - 1) + (offset / 2));
        r = radius*math.sqrt(radius2 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        xcoord.append(x);ycoord.append(y);zcoord.append(z)

    return xcoord,ycoord,zcoord

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    N = int(sys.argv[1])
    r = float(sys.argv[2])
    xc,yc,zc = fiboSph(N,r)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xc, yc, zc) #, c=c, marker=m)
    plt.show()


















