
# coding: utf-8

# # PHYS-E0412 Computational Physics -Homework 7
# Ari Viitala
# 432568

# ### a) Dimension 1. Implement a solver for the one-dimensional problem

# In[2]:


import numpy as np
import matplotlib.pyplot as plt


# In[3]:


def oned_solve(N):
    h = 1 / (N - 1)
    
    a = np.diag(np.ones(N))
    b = np.diag(np.ones(N-1), k = 1)
    c = np.diag(np.ones(N-1), k = -1)
    
    s = (- b + 2 * a - c)
    
    s[0,0] = 1
    s[0, 1] = 0
    s[-1,-1] = 1
    s[-1, -2] = 0

    
    s *= 1 / h**2
    
    x  = np.linspace(0, 1, N)
    b = (x - 0.5)**3 - 2*(x - 0.5)
    
    b[0] = 0
    b[-1] = 0
    
    return  np.linalg.solve(s, b)


# In[4]:


def analytical(N):
    x = np.linspace(0, 1, N)
    return -1./20 * (x-0.5) ** 5 + x**3 / 3. - x**2 / 2. +163 * x / 960. - 1./640 


# In[7]:


x = np.linspace(0,1,100)
plt.figure(1, (12, 8))
plt.plot(x, analytical(100), label = "Analytical")
for i in range(2, 10):
    N = i **2
    x  = np.linspace(0, 1, N)
    plt.plot(x, oned_solve(N), label = "N = " + str(i ** 2))

plt.legend()
plt.show()


# In[8]:


def error(N):
    h = 1 / (N - 1)
    return np.sqrt(h) * np.linalg.norm(analytical(N) -oned_solve(N))
    


# In[29]:


N = np.array(list(range(10, 1000, 50)))

errors = [error(i) for i in N]
plt.loglog(1 / (N-1), errors, marker = "x")
plt.show()


# In[17]:


np.polyfit(np.log(1 / (N-1)), np.log(errors), 1)


# In[20]:


def neumann_oned_solve(N):
    h = 1 / (N - 1)
    
    a = np.diag(np.ones(N))
    b = np.diag(np.ones(N-1), k = 1)
    c = np.diag(np.ones(N-1), k = -1)
    
    s = (- b + 2 * a - c)
    
    s[0,0] = 1
    s[0, 1] = 0
    s[-1,-1] = 1
    s[-1, -2] = -1

    s *= 1 / h**2
    
    x  = np.linspace(0, 1, N)
    b = (x - 0.5)**3 - 2*(x - 0.5)
    
    b[0] = 0
    b[-1] = 0
    
    return  np.linalg.solve(s, b)


# In[21]:


def neumann_analytic(N):
    x = np.linspace(0, 1, N)
    return -(1 / 20 * x**5 - 1/8* x**4 - 5/24* x**3 + 7/16* x **2)


# In[22]:


x = np.linspace(0,1,100)
plt.plot(x, neumann_analytic(100), label = "Analytical")
for i in range(1, 10):
    N = i * 10
    x  = np.linspace(0, 1, N)
    plt.plot(x, neumann_oned_solve(N), label = "N = " + str(i * 10))

plt.legend()
plt.show()


# In[23]:


def neumann_error(N):
    h = 1 / (N - 1)
    return np.sqrt(h) * np.linalg.norm(neumann_analytic(N) -neumann_oned_solve(N))
    


# In[31]:


N = np.array(list(range(10, 1000, 50)))

neumann_errors = [neumann_error(i) for i in N]
errors = [error(i) for i in N]
plt.loglog(1/(N-1), errors, marker = "x")
plt.loglog(1/(N-1), neumann_errors, marker = "x")
plt.show()


# In[25]:


np.polyfit(np.log(1 / (N-1)), np.log(neumann_errors), 1)


# In[72]:


N = 10

D = 4 * np.diag(np.ones(N)) - np.diag(np.ones(N-1), k = -1) - np.diag(np.ones(N-1), k = 1) 
I = -np.diag(np.ones(N))


# In[79]:


s = np.diag(np.zeros(N**2))

for i in range(0, N):
        if i == 0:
            s[:N,0:2*N] = np.hstack((D,I))
        elif i == N - 1:
            s[i*N:,N*(i-1):] = np.hstack((I,D))
        else:
            print(np.shape(s[i*N:(i+1)*N,N*(i-1):N*(i+2)]))
            s[i*N:(i+1)*N,N*(i-1):N*(i+2)] = np.hstack((I,D,I))
            
s *= 1 / (1 / (N-1)**2)


# In[96]:


b = np.zeros((N,N))
x = y = np.linspace(0,1,N)


# In[97]:


for i in range(N):
    for j in range(N):
        b[i, j] = np.exp(-((x[i] - 0.5)**2 + (y[j]-0.5)**2) / 18)
b[0,:] = 0
b[:,0] = 0
b[-1,:] = 0
b[:,-1] = 0
b = b.flatten()


# In[98]:


z = np.linalg.solve(s, b).reshape((N,N))


# In[99]:


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(x, y)
# Plot the surface
ax.plot_surface(X, Y, z, color='b')

plt.show()
