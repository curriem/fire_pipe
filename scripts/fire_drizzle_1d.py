import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

'''
0 = Cell is unburnable [U].
1 = Cell is flammable, but not currently ignited [N].
2 = Cell is flammable and is ignited, but fuel is not yet consumed [I].
3 = All fuel in cell has been consumed by the fire [C].

'''

correction = False
num_neighbors = 8

def init_grid(len_side, lake):
    state_grid = np.ones((len_side, 1))
    bucket_grid = np.zeros((len_side, 1))
    if lake:
        state_grid = np.load('lake2.npy')

    # start a point fire
    state_grid[50, 0] = 2
    bucket_grid[50, 0] = 1.

    return bucket_grid, state_grid

def init_elev_testcase(len_side):
    elev_grid = np.zeros((len_side, len_side))
    for n in range(int(len_side/2)):
        elev_grid[n:len_side-n, n:len_side-n] = n
    return elev_grid


def init_elev_grid(len_side, flat=False):
    len_side_div2 = len_side/2
    if flat: 
        elev_grid = np.zeros((len_side, 1))
    else:
        # gaussian mountain
        x, y = np.meshgrid(np.linspace(-len_side_div2, len_side_div2, len_side),
                           np.linspace(-len_side_div2, len_side_div2, len_side))
        d = np.sqrt(x*x+y*y)
        sigma, mu = 30.0, 1
        elev_grid = np.exp(-((d-mu)**2 / (2.0 * sigma**2)))
        elev_grid *= 100
    '''
    plt.figure()
    plt.imshow(elev_grid)
    plt.colorbar()
    plt.show()
    assert False
    '''
    return elev_grid


def init_wind_mag_grid(len_side, U_eq):
    wind_grid = U_eq*np.ones((len_side, 1))
    return wind_grid


def init_wind_dir_grid(len_side, wind90):
    wind_grid = np.pi * np.ones((len_side, 1))
    if wind90:
        wind_grid[:100, :] = np.pi/2
    return wind_grid


def init_fuel_grid(len_side, R0):
    fuel_grid = R0 * np.ones((len_side, 1))
    return fuel_grid


def rate_of_spread_xy(R0, Ebar, alpha, theta, ti):

    '''
    # old formulation:
    ros_x = (1) * R0*((1-Ebar)
                         / (1-Ebar*np.cos(theta-alpha))) * np.cos(-theta)
    ros_y = (1) * R0*((1-Ebar)
                         / (1-Ebar*np.cos(theta-alpha))) * np.sin(-theta)

    return np.abs(ros_x), np.abs(ros_y)
    '''
    # new formulation:
    ros = R0 * (1-Ebar) / (1-Ebar*np.cos(alpha-theta))

    ros_x = ros * np.cos(theta)
    ros_y = ros * np.sin(theta)

    return np.abs(ros_x), np.abs(ros_y)


def rate_of_spread(R0, U_eq, alpha, phi, ti):
    # from the Farsite paper: the original paper's units were wrong probably
    LW = 0.936*np.exp(0.2566*U_eq) + 0.461*np.exp(-0.1548*U_eq) - 0.397
    Ebar = np.sqrt(1 - 1/LW**2)
    theta = phi - alpha

    if correction:

        # corrected rate of spread
        k1 = 1
        k2 = 2
        k3 = 0
        k4 = 0
        k5 = 0

        # k values for Ebar = 0.5:
        k1 = 1.
        k2 = 1.15
        k3 = 1.850
        k4 =  1.7
        k5 = 1.55

        R0_prime = k1*R0
        Ebar_prime = np.sqrt(1 - 1/LW**k2)
        if (k3 == 0) and (k4 == 0):
            theta_prime = theta
        else:
            if np.abs(theta) <= np.pi/2:
                theta_prime = np.arctan(np.tan(theta) / LW**k3)
            elif np.abs(theta) > np.pi/2:
                theta_prime = np.arctan(np.tan(theta) / LW**k4)

            else:
                assert False, 'weird theta'
        ros = (R0_prime * (1-Ebar_prime) / (1-Ebar_prime*np.cos(theta_prime))
               - R0_prime*((1-Ebar_prime)/(1+Ebar_prime)
                           - (1-Ebar)/(1+Ebar))
               * (np.abs(theta_prime)/np.pi)**k5)

        # assert Ebar_prime == Ebar, '%.03f =/= %.03f' % (Ebar_prime, Ebar)
        # assert R0_prime == R0, '%.03f =/= %.03f' % (R0_prime, R0)
        # assert theta_prime == theta, '%.03f =/= %.03f' % (theta_prime, theta)

    else:
        ### HARD CODED DELETE LATER ####
        Ebar = 0.5
        ros = R0 * (1-Ebar) / (1-Ebar*np.cos(theta))

    return ros


def phi_corrections(slope, beta=50000, beta_op=1, B=1, C=0, U=1, E=1):
    phi_s = np.sign(slope)*5.275 * beta**-0.3 * np.tan(slope)**2
    phi_w = C*(3.281*U)**B * (beta/beta_op)**-E
    if np.isnan(phi_s):
        print 'error:', phi_s
        print slope
        assert False
    return phi_s, phi_w




def horizontal_component(elev_grid_neighbors):

    neighbor_angles = [[3*np.pi/4, np.pi/2, np.pi/4],
                       [np.pi, 0, 0],
                       [5*np.pi/4, 3*np.pi/2, 7*np.pi/4]]

    pix_of_interest = elev_grid_neighbors[1, 1]

    delt_elev_grid = elev_grid_neighbors - pix_of_interest

    delt_x_grid = np.array([[-1, 0, 1],
                            [-1, 0, 1],
                            [-1, 0, 1]])

    delt_y_grid = np.array([[1, 1, 1],
                            [0, 0, 0],
                            [-1, -1, -1]])

    x_ang_grid = np.tan(delt_elev_grid / delt_x_grid)

    y_ang_grid = np.tan(delt_elev_grid / delt_y_grid)

    h_correction = np.sqrt(np.cos(neighbor_angles)**2 * np.cos(x_ang_grid)**2
                           + np.sin(neighbor_angles)**2 * np.cos(y_ang_grid)**2)

    return h_correction


def calculate_slope_neighbors(elev_neighbors):
    pix_of_interest_elev = elev_neighbors[1]

    delt_elev_neighbors = elev_neighbors - pix_of_interest_elev

    distance_neighbors = [1, 0, 1]

    slope_neighbors = delt_elev_neighbors / distance_neighbors
    return slope_neighbors


def fire_colormap():
    fire_cmap = colors.ListedColormap(['blue', 'limegreen',
                                       'red', 'k'])
    bounds = [0, 1, 2, 3, 4]
    norm = colors.BoundaryNorm(bounds, fire_cmap.N)

    return fire_cmap, norm


def plot_frame(frame, savepath):
    fire_cmap, norm = fire_colormap()

    plt.figure()
    plt.imshow(frame, cmap=fire_cmap, norm=norm)
    plt.axis('off')
    plt.colorbar()
    plt.savefig(savepath)
    plt.close()



def run_code(flat, pyramid, lake, wind90, save_name):

    # params
    U_eq = 0.78
    R0 = 1
    timesteps = 140
    len_side = 100

    neighbor_angles = [0, np.pi]
    neighbor_locations = np.array([(+1, 0),
                                   (-1, 0)])
    neighbor_distances = np.array([1, 1])

    bucket_grid, state_grid = init_grid(len_side, lake)
    wind_mag_grid = init_wind_mag_grid(len_side, U_eq)
    fuel_grid = init_fuel_grid(len_side, R0)
    wind_dir_grid = init_wind_dir_grid(len_side, wind90)
    if flat:
        elev_grid = init_elev_grid(len_side, flat=True)
    if pyramid:
        elev_grid = init_elev_testcase(len_side)

    plot_frame(state_grid, '../plots/im000.png')
    timeseries = []
    for ti in range(timesteps):
        print 'working on timestep', ti

        # set burnt inds with this loop:
        I_inds_x, I_inds_y = np.where(state_grid == 2)
        for n in range(len(I_inds_x)):
            if np.all(state_grid[I_inds_x[n]-1:I_inds_x[n]+2,
                                 I_inds_y[n]] != 1):
                # if the pixel is surrounded by all burnt, unburnable,
                # or ignited pixels, set it to burnt
                state_grid[I_inds_x[n], I_inds_y[n]] = 3

        # continue lighting with this loop:
        I_inds_x, I_inds_y = np.where(state_grid == 2)
        for n in range(len(I_inds_x)):
            elev_neighbors = elev_grid[I_inds_x[n]-1:I_inds_x[n]+2,
                                       I_inds_y[n]]

            slope_neighbors = calculate_slope_neighbors(elev_neighbors)
            #except ValueError:
            #    continue
    #        h_correction = horizontal_component(elev_neighbors)

            # OLD WAY
            for neighbor_num, phi in enumerate(neighbor_angles):
                ros = rate_of_spread(fuel_grid[I_inds_x[n], I_inds_y[n]],
                                     wind_mag_grid[I_inds_x[n], I_inds_y[n]],
                                     wind_dir_grid[I_inds_x[n], I_inds_y[n]],
                                     phi, ti)

                print ros
                # correct for horizontal view
                # THE H CORRECTION CODE IS NOT RIGHT!!!!
                # ros *= h_correction[neighbor_locations[neighbor_num, 0],
                #                     neighbor_locations[neighbor_num, 1]]
                phi_s, phi_w = \
                phi_corrections(slope_neighbors[neighbor_locations[neighbor_num, 0]
                                                + 1])
                ros = ros * (1 + phi_s + phi_w)
                bucket_grid[I_inds_x[n]
                            + neighbor_locations[neighbor_num, 0],
                            I_inds_y[n]
                            + neighbor_locations[neighbor_num, 1]] \
                    += (ros / neighbor_distances[neighbor_num])
                new_I_inds_x, new_I_inds_y = np.where((bucket_grid >= 1.)
                                                      & (state_grid ==
                                                         1))
                state_grid[new_I_inds_x, new_I_inds_y] = 2
            '''

            # NEW WAY, vectorized
            neighbor_nums = np.arange(0, len(neighbor_angles))
            ros = rate_of_spread(fuel_grid[I_inds_x[n], I_inds_y[n]],
                                 wind_mag_grid[I_inds_x[n], I_inds_y[n]],
                                 wind_dir_grid[I_inds_x[n], I_inds_y[n]],
                                 neighbor_angles, ti)
            bucket_grid[I_inds_x[n]
                        + neighbor_locations[neighbor_nums, 0],
                        I_inds_y[n]
                        + neighbor_locations[neighbor_nums, 1]] \
                += (ros / neighbor_distances[neighbor_nums])
            new_I_inds_x, new_I_inds_y = np.where((bucket_grid >= 1.)
                                                  & (state_grid == 1))
            state_grid[new_I_inds_x, new_I_inds_y] = 2
        '''
        temp = np.copy(state_grid)
        timeseries.append(temp)

        plot_frame(state_grid, '../plots/im%s.png' % str(ti+1).zfill(3))

    timeseries = np.array(timeseries)
    np.save('../model_timeseries/%s.npy' % save_name, timeseries)


if __name__ == '__main__':
    run_code(flat=True, pyramid=False, lake=False, wind90=False,
             save_name='flat_1d')



