import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

#spatial grid
N = 100
h = 1 / (N - 1)
x = np.linspace(0,1, N)

#initial state
u0 = 1.831 * np.exp(-10 * (x - 0.5)**2)
u0 -= u0[0]

#time
t_steps = 100
time = np.linspace(0, 0.1, 100)[1:]
t_h = time[1] - time[0]

#the three point stencil
stencil =  -1 / h**2 * (2 * np.eye(N) - np.diag(np.ones(N-1), 1) - np.diag(np.ones(N-1), -1))  

#boundary conditions
stencil[0] = 0
stencil[-1] = 0
stencil[0, 0] = 1
stencil[-1, -1] = 1

#vector to store the values
u = [u0]
#for each timestep
for i, t in enumerate(time):
    #solve the system for the next state
    u1 = np.linalg.solve(np.eye(N) - t_h * stencil, u[i])
    #append to the result vecotr
    u.append(u1)

#convert to numpy array
u= np.array(u)

fig = plt.figure(1, (10, 7)) 
ax = fig.add_subplot(111, projection='3d') 

time_ax = [0] + list(time)
X, Y = np.meshgrid(x, time_ax) # Plot the surface 
ax.plot_surface(X, Y, u, color='b') 
plt.show() 

#result vector
u_cn = [u0]
#the implicit side
first = np.linalg.solve(np.eye(N) - 0.5 * t_h * stencil, np.eye(N))
#the explicit side
second = np.eye(N) + 0.5 * t_h * stencil
m = first @ second
for i, t in enumerate(time):
    #the next step is just matrix multiplication
    u1 = m.dot(u_cn[i])
    u_cn.append(u1)
    
u_cn = np.array(u_cn)

fig = plt.figure(1, (10, 7)) 
ax = fig.add_subplot(111, projection='3d') 

time_ax = [0] + list(time)
X, Y = np.meshgrid(x, time_ax) # Plot the surface 
ax.plot_surface(X, Y, u_cn, color='b') 
plt.show()

def sten(N):
    D = 4 * np.diag(np.ones(N)) - np.diag(np.ones(N-1), k = -1) - np.diag(np.ones(N-1), k = 1) 
    
    #boundary conditions
    D[0] = D[-1] = 0
    D[0,0] = D[-1,-1] = 1
    
    #identity matrix for off diagonals
    I = np.diag(np.ones(N))
    IO = np.copy(I)
    
    #boundary conditions
    IO[0] = IO[-1] = 0
    
    s = np.diag(np.zeros(N**2))
    
    #filling the stencil with block matrices
    for i in range(0, N):
            if i == 0:
                #first row
                s[:N,0:N] = I
            elif i == N - 1:
                #last row
                s[i*N:,N*i:] = I
            else:
                #other rows
                s[i*N:(i+1)*N,N*(i-1):N*(i+2)] = np.hstack((-IO,D,-IO))
    
    #scaling
    s *= 1 / (1 / (N-1)**2)
    s += 2*np.eye(N*N)
    return s