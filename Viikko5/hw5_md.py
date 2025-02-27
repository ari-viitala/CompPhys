import numpy as np

#For plotting
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#Reduce all particle positions to [0,L]^dim box
#Mostly to aid visualization if needed
def particle_reduced_position(r, L):
    rr = r.copy()
    rr -= np.floor(rr/L)*L
    return rr

def visualize_particles_3d(r, L, title=""):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)
    ax.scatter(r[:,0], r[:,1], r[:,2], c='b', marker='o')
    plt.axis(xmin=0.0, xmax=L[0], ymin=0.0, ymax=L[1])
    ax.set_xlim(0.0, L[0])
    ax.set_ylim(0.0, L[1])
    ax.set_zlim(0.0, L[2])
    #plt.show()

#We need three essential pieces for MD:
#1. Periodic boundary conditions, provided here by minimum image convention
#2. Forces from distances where PBC has been applied
#3. Integration of equations of motion. Provided here by Velocity Verlet

#Below all of these provided. Note that nonvectorized versions are quite slow,
#but they are provided for clarity. They do exactly the same thing as vectorized versions
def reduced_distance_nonvectorized(_r, L):
    N = _r.shape[0]
    invL = 1.0/L
    r_ij = np.zeros((N,N,_r.shape[1]))
    for i in range(1,N):
        for j in range(i):
            r = _r[i] - _r[j]
            r -= np.round(r*invL)*L
            r_ij[i,j,:] = r
            r_ij[j,i,:] = -r
    return r_ij

#Calculate reduced distance with minimum image convention
#https://en.wikipedia.org/wiki/Periodic_boundary_conditions
def reduced_distance(r, L):
    dim = r.shape[1]
    N = r.shape[0]
    invL = 1.0/L

    A = np.tile(r, N).reshape(N,N,dim)
    AT = np.transpose(A, (1,0,2))
    dr = A - AT
    dr -= np.round(dr*invL)*L
    return dr


def do_velocity_verlet_integration(r, v, f, dt, L):
    r = r + dt*v + (0.5*dt**2)*f
    v += (0.5*dt)*f
    dr = reduced_distance(r,L)
    f,potential = md_force(dr)
    v += (0.5*dt)*f
    return r,v,f,dr,potential

#TODO: implement LJ force with epsilon = sigma = 1
#dr_ijk is N x N x 3 array (N is number of particles), which contain pair-wise
#separation r_i - r_j where minimum image convention has been applied.
#Please see implementation of reduced_distance for more information
def md_force(dr):
    N = dr.shape[0]
    f = np.zeros((N,3))
    dists = np.nansum(dr**2, axis = 2).reshape((N,N,1))

    #if we dont take the square root at the dists then we can raise the denominators
    #only to powers 7 and 4 when calculating the force
    f =  48 * np.nansum(dr * (1 / dists**7 - 1 / (2 * dists**4)), axis = 1)
    potential_energy = 2 * np.nansum(1 / dists**6 - 1 / dists**3)
    return f, potential_energy


def run_md_3d(r,L,T=0.0, dt=5e-3, thermalization_time=2.0, production_time=10.0):
    N = r.shape[0]
    dim = r.shape[1]
    v = np.zeros((N, dim))

    vs = np.random.normal(0, 1, (N, dim))
    norms = np.linalg.norm(vs, axis = 1).reshape(N, 1)

    v = vs / norms * np.sqrt(3 * T)
    print("Initial temp: {}".format(np.sum(v**2) / N / 3))
    #v = np.zeros((N, dim))

    #TODO: Implement initialization to given temperature T!

    rho = N/np.prod(L)
    print("Density: {}".format(rho))
    print("Box edge length: {}".format(L[0]))
    print("Box diameter: {}".format(np.linalg.norm(L)))

    num_thermalization_steps = int(thermalization_time/dt)
    num_production_steps = int(production_time/dt)

    #Velocity verlet integration
    dr = reduced_distance(r,L)
    f, potential_energy = md_force(dr)
    for i in range(num_thermalization_steps):
        r,v,f,dr,potential_energy = do_velocity_verlet_integration(r, v, f, dt, L)

        #Cooling the system during the theramalization to maybe help it stay together
        if i % 50 == 0:
            avg_speed = np.sqrt(v**2).mean()
            v *= T / avg_speed

        if i % 100 == 0:
            print("Thermalization: {}/{}".format(i*dt, thermalization_time))

    #Initialize quantities we want to sample
    num_samples = num_production_steps
    t = np.linspace(0, production_time, num=num_production_steps)
    Ekin = np.zeros((num_samples))
    Epot = np.zeros((num_samples))
    rdf = np.zeros((N, N, num_samples))
    #TODO: Implement RDF sampling. You can e.g. sample all the distances at all times
    #and postprocess them after simulation

    #Velocity verlet integration and sampling
    for i in range(num_production_steps):
        r,v,f,dr,potential_energy = do_velocity_verlet_integration(r, v, f, dt, L)
        #sample index
        si = i
        Ekin[si] = 0.5*np.sum(v**2)/N
        Epot[si] = potential_energy/N

        #calculating the distances in the system
        dists = np.sqrt(np.sum(dr**2, axis = 2))

        #saving the upper triangle of the distance matrix
        rdf[:,:,si] = np.triu(dists, 1)


        # if i % 10 == 0:
        #     visualize_particles_3d(r, L)

        if i % 100 == 0:
            print("Production: {}/{}".format(i*dt, production_time))

    result = {}
    result["t"] = t
    result["Ekin"] = Ekin
    result["Epot"] = Epot
    result["Final"] = r
    result["rdf"] = rdf

    return result

def rdf(results, bins_count):
    """Function for calculating the radial distribution function"""

    #take the mean of the distances over time
    mean = np.mean(results, axis = 2)

    #make a histogram of the data and removing zeros
    hist = np.histogram(mean[np.nonzero(mean)], bins = bins_count)

    #the heights of the histogram bars
    heights = hist[0]

    bins = hist[1]

    binw = bins[1] - bins[0]

    #the locations of the height values
    rs = bins[:-1] + binw / 2

    #the nomalization term for g(r)
    norm = (N - 1) / (rho * (4 * np.pi) * np.trapz(heights, rs))

    #using the normalization term and the spherical coordinate transform to obtain the correct g(r)
    gr = norm * heights / rs**2

    #calculating the integral to verify that this g(r) fulfills the normalization requirement
    integral = np.trapz(gr * rs**2, rs) * rho * 4 * np.pi

    print("Trapz value: " + str(integral))

    #reuturn the function and the x-axis values
    return gr, rs

#By default initialize one particle in 4x4x4 box doing nothing
if __name__ == "__main__":
    rho = 0.8442
    N = 108
    #Position array is always array of size N x 3
    #where N is the number of particles and 3 is the dimension

    #making the box
    L = np.array([1.0,1.0,1.0]) * np.power(N / rho, 1/3.)

    #making the r-vector
    r = np.zeros((125,3))

    x = 5
    y = 5
    z = 5

    #populating the grid using for-loops
    for i in range(0, x):
        for j in range(0, y):
            for k in range(0, z):
                idx = (i * y + j) * z + k
                r[idx,:] = np.array([i, j, k], dtype = float)
                #print(r.shape)

    #remocing random elements to achieve the desired density
    while(len(r) > N):
        r = np.delete(r, np.random.randint(0, len(r)), axis = 0)
        #print(r.shape)

    #scaling the distances
    r *= L[0] / 5
    #r += (np.random.rand(N, 1) - 0.5) * 0.01
    #r = np.random.rand(N, 3) * L[0]
    #L = np.array([4.0,4.0,4.0])
    #r = np.array([[0,0,0], [0,0,1.2]], dtype = float)
    #r[0] = L/2.0
    #r[0] = L/1.4

    #Simulation for the b)-part
    result = run_md_3d(r,L, T = 0.728, thermalization_time=2, production_time=5.0)

    #Simulation for the C-part
    #result = run_md_3d(r,L, T = 2, thermalization_time=2, production_time=5.0)

    plt.figure()
    plt.plot(result["Ekin"], label="Kinetic energy")
    plt.plot(result["Epot"], label="Potential energy")
    plt.plot(result["Ekin"] + result["Epot"], label="Total energy")
    #plt.plot(result["Ekin"] * 2. / 3)
    plt.legend()

    #calculating the rdf and plotting it
    gr, rs = rdf(result["rdf"], 70)

    plt.figure()
    plt.plot(rs, gr)
    plt.xlabel("r")
    plt.ylabel("$g(r)$")
    print("Final temp: {}".format(result["Ekin"][-1] * 2 / 3.0))
    visualize_particles_3d(r,L)
    visualize_particles_3d(result["Final"], L)
    plt.show()
