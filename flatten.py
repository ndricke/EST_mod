import sys
import numpy
import matplotlib.pylab as plt
import matplotlib

bohr = 0.529177

# settings
bond_width = 2
bond_cut = 1.55 / bohr

atom_size = 50

# preamble
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Computer Modern Roman'
plt.rcParams['legend.numpoints'] = 1

arrayify = lambda x: numpy.array(map(float,x))
fname = sys.argv[1]
f = open(fname)
f.readline()
f.readline()

s = f.readline().split()
n = int(s[0]) #number of atoms
origin = arrayify(s[1:4]) #coordinates for the origin

#generates a 3x3 array of the voxel vectors in cube file
vox = numpy.zeros((3,3))
s = f.readline().split()
nx = int(s[0])
vox[0,:] = arrayify(s[1:4])
s = f.readline().split()
ny = int(s[0])
vox[1,:] = arrayify(s[1:4])
s = f.readline().split()
nz = int(s[0])
vox[2,:] = arrayify(s[1:4])

xl = numpy.linalg.norm(vox[0,:])
yl = numpy.linalg.norm(vox[1,:])

ivox = numpy.linalg.inv(vox)

atoms = []
for i in range(n):
    s = f.readline().split()
    r = arrayify(s[2:5])
    cr = r.copy()
    r -= origin
    r = r.dot(vox.T)
    r[0] /= xl/0.529177
    r[1] /= yl/0.529177
    atoms.append({'nr':int(s[0]), 'r':r, 'cr':cr})

# build bonds from distance search

#bond_omit = {}
#bond_omit['RO_200.cub'] = [(
bonds = []
for i in range(n-1):
    for j in range(i+1,n):
        if numpy.linalg.norm(atoms[i]['cr']-atoms[j]['cr']) < bond_cut:
            bonds.append((i,j))

#It looks like we are actually loading this juggernaut of a file into memory
data = numpy.zeros(nx*ny*nz)
idx = 0
for line in f:
    for datum in map(float, line.split()):
        data[idx] = datum
        idx += 1
print idx

data = data.reshape([nx,ny,nz])
data = numpy.abs(data)
#data = data*data
data = data.sum(axis=2)
print data.max(), data.min()

floor = 0.0
#ceiling = 0
#floor = 0
ceiling = 1.5
#ceiling = 0.3


xgrid = numpy.array(range(nx))*xl*0.529177
ygrid = numpy.array(range(ny))*yl*0.529177

for x in range(nx):
    for y in range(ny):
        data[x,y] = min(data[x,y], ceiling)
        data[x,y] = max(data[x,y], floor)

atomx = numpy.array([atom['r'][0] for atom in atoms])
atomy = numpy.array([atom['r'][1] for atom in atoms])

minx, maxx = atomx.min(), atomx.max()
miny, maxy = atomy.min(), atomy.max()

# Attempt to build fancy map
def cmap_map(function,cmap):
    from numpy import array
    import matplotlib
    """ Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         
        step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : array(cmap(step)[0:3])
    old_LUT = array(map( reduced_cmap, step_list))
    old_LUT[0,:] = 1
    new_LUT = array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

#mymap = cmap_map(lambda x: x/2+.5, plt.get_cmap('jet'))
#mymap = plt.get_cmap('jet')
#mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mymap', [
#    (0, 'lime'),
#    (0.25, 'white'),
#    (.5, 'skyblue'),
#    (.75, 'blue'),
#    (1, 'red')])
#mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mymap', [
#    (0, 'firebrick'),
#    #(.125, 'firebrick'),
#    (0.25, 'white'),
#    (.5, 'mediumpurple'),
#    (1, 'skyblue')])
mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mymap', [
    (0, 'white'),
    (.25, 'skyblue'),
    (.5,'violet'),
    (.75, 'orange'),
    (1.0, 'red')])
#mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mymap', [
#    (0, 'white'),
#    (0.15, 'purple'),
#    (0.3, 'cyan'),
#    (0.65, 'yellow'),
#    (0.90, 'red'),
#    (1.0,'darkred')])

plt.contourf(xgrid, ygrid, data.T, levels=numpy.linspace(floor,ceiling,100), cmap=mymap)
plt.colorbar(ticks = numpy.arange(0,1.6,0.3), orientation="horizontal",fraction=0.0598, pad = 0.08)
#plt.colorbar(ticks = numpy.arange(0,.31,0.1), orientation="horizontal",fraction=0.0598, pad = 0.08)


# build color list
cmap = {6: 'gray',
        1: 'white',
        8: 'red',
        7: 'blue',
        9: 'pink'}
colors = [cmap[atom['nr']] for atom in atoms]

# Draw bonds
for bond in bonds:
    r1 = atoms[bond[0]]['r'] 
    r2 = atoms[bond[1]]['r'] 
    plt.plot([r1[0], r2[0]], [r1[1],r2[1]], '-', color='.25', lw=bond_width, zorder=1)

plt.scatter(atomx, atomy, c=colors, s=atom_size, zorder=2)


# TODO:
# - same limits every plot
# - square aspect ratio
plt.xlim(minx-2,maxx+2)
plt.ylim(miny-2,maxy+2)
print plt.xlim()
print plt.ylim()
plt.xlim(6.6398969977882984, 31.557921302998952)
plt.ylim(4.0850987544726207, 21.823688494451083)

plt.axes().set_aspect('equal', 'box')
plt.xticks([])
plt.yticks([])
plt.tight_layout()
plt.show()

