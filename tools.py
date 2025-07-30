import numpy as np

# __all__ = ["xyzRead", "xyzWrite"]


def xyzRead(fname):
    fin = open(fname, "r")
    line1 = fin.readline().split()
    natoms = int(line1[0])
    comments = fin.readline()[:-1]
    coords = np.zeros([natoms, 3], dtype="float64")
    atomtypes = []
    for x in coords:
        line = fin.readline().split()
        atomtypes.append(line[0])
        x[:] = list(map(float, line[1:4]))

    return natoms, atomtypes, coords;

def getCharge(element):
    f = open("mol.txt")
    atomicnum = [line.split()[1] for line in f if line.split()[0] == element]
    f.close()
    return int(atomicnum[0])

def coulombMatrix(fname):
    natoms, atomtypes, coords = xyzRead(fname)
    i=0 ; j=0    
    colM = np.zeros((natoms,natoms))
    chargearray = np.zeros((natoms,1))
    charge = [getCharge(symbol)  for symbol in atomtypes]
    for i in range(0,natoms):
        colM[i,i]=0.5*charge[i]**2.4   # Diagonal term described by Potential energy of isolated atom
        for j in range(i+1,natoms):
            # print(str(i)+" - "+str(j))
            dist= np.linalg.norm(coords[i,:] - coords[j,:])   
            colM[j,i] = charge[i]*charge[j]/dist   #Pair-wise repulsion 
            colM[i,j] = colM[j,i]
    return colM
 
# def coulombMatrixSort(fname):
#     uCoulomb = coulombMatrix(fname)
#     sumLine = np.array([sum(x**2) for x in uCoulomb])
#     sCoulomb = uCoulomb[np.argsort(sumLine)[::-1,],:]
#     return sCoulomb
#     # return sCoulomb.ravel()

def eigenCoulomb(fname, num):
    sCoulomb = coulombMatrix(fname)
    sCoulomb = sCoulomb.astype(int)
    # print("\n\n")
    # print(sCoulomb)
    # print("Sym: ")
    # print(check_symmetric(sCoulomb))
    # print("\n\n")
    eigValues = -np.sort(-np.linalg.eigvals(sCoulomb))
    # eigValues = np.linalg.eigvals(sCoulomb)
    # print(eigValues)
    return eigValues[0:num]

def evolveParticles(num_p):
    particles = np.random.rand(num_p,3)
    for t in range(1000):
        step = np.zeros((num_p,3))
        dmin = 2
        for i in range(0,num_p-1):
            for j in range(i+1,num_p):
                diff = particles[i,:] - particles[j,:]
                norm = np.linalg.norm(diff)
                dist = np.exp(-norm*4)
                step[i,:] += dist*diff / norm
                step[j,:] -= dist*diff / norm
                if norm<dmin: dmin=norm


        particles += step
        soma = 0
        for i in range(num_p):
            norm = np.linalg.norm(particles[i,:])
            particles[i,:] /= norm
            soma += norm
        print("T: "+ str(t)+" -> "+str(soma/num_p)+" Step: "+str(dmin))

    return particles

def distMin(xyz, num_p):
    dist = np.ones((num_p))*10
    for i in range(num_p):
        for j in range(num_p):
            if j==i: continue
            norm = np.linalg.norm(xyz[i,:] - xyz[j,:])
            if norm < dist[i]:
                dist[i] = norm
    return dist


# 

# def write_xyz(fout, coords, title="", atomtypes=("A",)):
#     """ write a xyz file from file handle

#     Writes coordinates in xyz format. It uses atomtypes as names. The list is
#     cycled if it contains less entries than there are coordinates,

#     One can also directly write xyz data which was generated with read_xyz.

#     >>> xx = read_xyz("in.xyz")
#     >>> write_xyz(open("out.xyz", "w"), *xx)

#     Parameters
#     ----------
#     fout : an open file
#     coords : np.array
#         array of coordinates
#     title : title section, optional
#         title for xyz file
#     atomtypes : iteratable
#         list of atomtypes.

#     See Also
#     --------
#     read_xyz

#     """
#     fout.write("%d\n%s\n" % (coords.size / 3, title))
#     for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
#         fout.write("%s %.18g %.18g %.18g\n" % (atomtype, x[0], x[1], x[2]))
#  
