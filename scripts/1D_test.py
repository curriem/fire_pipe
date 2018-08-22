
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

kappa = 0.5e5 # cm ^ 2 / sec, local radiation and turbulent heat diffusion
L = 10.0e5 # len of pipe, cm
robert_filter_coeff = 0.25
btime = 1000.0 # s
rtime = 5000.0 # linear radiation damping time s
sigma = 5.6704e-5 # erg * cm^-2 * s^-1 * K^-4

# advecting velocity
u = 3.5 # 2 # cm / sec

Hc = 0.0 # heat of combustion J / ms not used yet - revise(Q~ Hc * u ^ alpha)
# add nondim 'l combustion parameter fcn of Hc, u, slope...

Pe = u * L / kappa # Bulk Peclet  # a controlling nondimensional parameter

# period of integration:
year = 360.0 * 84000.0
day = 84000.
T = 5.0 * day # seconds

# finite differencing parameters:
nx = 100 # number of grid points along - pipe(along x direction)
dx = L / (nx - 1) # size of grid cell
dt = 0.002 * day # time step(seconds)
nt = np.floor(T / dt) # number of time steps

Pelocal = u * dx / kappa # Local Peclet  # a controlling nondimensional parameter

# Initial conditions:
C = np.zeros((nx,3))
F=np.ones(nx)
for i in np.arange(0,nx):
    for j in np.arange(0,3):
        C[i][j] = np.exp(-(i+1 - nx / 2) ** 2 / (nx / 8) ** 2)

pipe2d = np.zeros((int(nt), int(nx))) # concentration field
time = np.zeros((int(nt)))
F[0] = 0.0 # no fuel at the ends
F[nx-1] = 0.0


# Integrate in time: n is one of:
# 0 - past(known)
# 1 - present(known)
# 2 - future(to be calculated)

for n in np.arange(0,nt):
# Loop over all space points and find future value
    for i in np.arange(1,nx-2):
        Q = 0.0
        if C[i][1] > 0.7: # here 0.7 = ignition temperature

            if F[i] > 0.1: # heat source if there is fuel
                Q = 0.0015 # should depend on wind, fuel type etc

            F[i] = F[i] - dt * F[i] / btime # fuel mass decay end
            print(F[i])

        C_x = (C[i+1][1] - C[i-1][1]) / (2 * dx)

        C_xx = (C[i+1][0] - 2 * C[i][0] + C[i -1][0]) / (dx * dx)

        C[i][2] = C[i][0] + 2 * dt * (-u * C_x + kappa * C_xx + Q - (sigma * C[i][0]**4))

    # cyclic boundary coditions for C at ends of pipe:
    C[0, 2] = 0 # C(nx - 1, 3)
    C[nx-1][2] = C[nx - 2][2] # 0. # 1 # C(2, 3)

    # Robert filter:
    C[:, 1] = C[:, 1]+robert_filter_coeff * (C[:, 2] - 2 * C[:, 1] + C[:, 0])
    # prepare for next time step:
    C[:, 0]=C[:, 1]
    C[:, 1]=C[:, 2]

    # store solution for plotting:
    time[int(n)] = n+1 * dt
    pipe2d[int(n),:]=C[:, 2]



"""""""""
plt.figure(1) # open figure 1 if it is not already open
# plot concentration at center of pipe as function of time:
plt.plot(time, pipe2d[:, int(nx / 2)])
plt.title('C(t)')
plt.xlabel('time (sec)')
plt.ylabel('Concentration at center')
plt.grid()
"""

plt.figure(2)
# contour tracer concetration as function of time and space:
plt.contour(pipe2d)
#plt.text(0,450,"Local radiation and turbulent heat diffusion (kappa; cm^2/s): "+str(kappa))
#plt.text(40,400,"Advecting Velocity (cm/s): "+str(u))
#plt.text(40,350,"Peclet number: "+str(Pe))
#plt.text(40,300,"Linear radiation damping time (s): "+str(rtime))
plt.xlabel('Distance (cm)')
plt.ylabel('Time (s)')
plt.title('Pipe Temperature Progression Map')
plt.ylim(0,650)
plt.grid()

# Second figure for animation
fig=plt.figure(3)
ax = plt.axes(xlim=(0,100), ylim=(0,10))
x = np.arange(0, len(pipe2d[1]))
y = pipe2d[1]
h, = ax.plot(x, y, lw=2)
plt.xlabel='Pipe Position (cm)'
plt.ylabel='Temperature Concentration'

def init():
    h.set_data([], [])
    return h,

def animate(i):
    x = np.arange(0, len(pipe2d[1]))
    y = pipe2d[i,:]
    h.set_data(x, y)
    plt.title('t='+str(i)+' (s)')
    plt.xlabel = 'Pipe Position (cm)'
    plt.ylabel = 'Temperature Concentration'
    return h,

anim = animation.FuncAnimation(fig, animate,init_func=init, frames=200, blit=False)
anim.save('firepipe1D_radiation.mp4', fps=20)

plt.show()
