import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

##############################################################################################


class Grid1D:
    def __init__(self, domain=(-10.0, 10.0), dx=20/999, **kwargs):
        count = int((domain[1] - domain[0])/dx + 1)
        self.x = np.linspace(domain[0], domain[1], count)
        self.dx = dx
        self.y = np.zeros(count, dtype=complex)

    def __add__(self, other):
        other_y = other.y if isinstance(other, Grid1D) else other
        result = Grid1D(domain=(self.x[0], self.x[-1]), dx=self.dx)
        result.y = self.y + other_y
        return result

    def sec_deriv(self):
        d2ydx2 = np.zeros_like(self.y)
        # central difference for the interior points
        d2ydx2[1:-1] = (self.y[:-2] - 2*self.y[1:-1] + self.y[2:])/(self.dx**2)
        # forward difference for the first point
        d2ydx2[0] = (self.y[2] - 2*self.y[1] + self.y[0])/(self.dx**2)
        # backward difference for the last point
        d2ydx2[-1] = (self.y[-1] - 2*self.y[-2] + self.y[-3])/(self.dx**2)
        return d2ydx2

    def rk4_step(self, dydt, dt, boundary="fixed", *args, **kwargs):
        k1 = dydt(self, *args, **kwargs)
        k2 = dydt(self + k1*dt/2, *args, **kwargs)
        k3 = dydt(self + k2*dt/2, *args, **kwargs)
        k4 = dydt(self + k3*dt, *args, **kwargs)
        self.y += (k1 + 2*k2 + 2*k3 + k4)*dt/6
        if boundary == "fixed":
            self.y[0] = 0
            self.y[-1] = 0
        elif boundary == "periodic":
            pass
        else:
            pass

    def inner_prod(grid_1, grid_2):
        return np.sum(np.conjugate(grid_1.y)*grid_2.y*grid_1.dx)


##############################################################################################


class Particle:
    def __init__(self, mass=1.0, h_bar=1.0, **kwargs):
        self.mass = float(mass)
        self.h_bar = float(h_bar)
        self.psi = Grid1D(**kwargs)
        self.potentials = []
        self.time = 0.0
        self.running = True

    def from_initial(func, *args, **kwargs):
        particle = Particle(**kwargs)
        particle.psi.y = func(particle.psi.x, *args, **kwargs)
        particle.normalize()
        return particle

    def add_initial(self, func, *args, **kwargs):
        self.psi.y += func(self.psi.x, *args, **kwargs)
        self.normalize()
        return self

    def add_potential(self, func, *args, **kwargs):
        def v(t, x): return func(t, x, *args, **kwargs)
        self.potentials.append(v)
        return self

    def potential(self):
        V = np.zeros_like(self.psi.x)
        for pot in self.potentials:
            V += pot(self.time, self.psi.x)
        return V

    def normalize(self):
        norm = np.sqrt(np.sum(np.abs(self.psi.y)**2 * self.psi.dx))
        self.psi.y /= norm

    def step(self, dt, before_step=None, **kwargs):
        if self.running:
            V = np.zeros_like(self.psi.x)
            for pot in self.potentials:
                V += pot(self.time, self.psi.x)

            def dydt(psi):
                momen = 1j*self.h_bar*psi.sec_deriv()/(2*self.mass)
                pot = 1j*self.potential()*psi.y/self.h_bar
                return momen - pot

            before_step(self) if before_step is not None else None
            self.psi.rk4_step(dydt, dt, **kwargs)
            self.normalize()
            self.time += dt

    def imaginary_step(self, dt, **kwargs):
        self.step(-1j*dt, **kwargs)

    def pause_play(self):
        self.running = not self.running

    def animate(self, dt=1/3600, anim_length=10, anim_speed=1, x_lim=None, y_lim=None, notebook=False, **kwargs):
        fig = plt.figure("Schrodinger Simulation")
        ax = fig.add_subplot()

        # set the limits for the x and y axes
        x_lim = (self.psi.x[0], self.psi.x[-1]) if x_lim is None else x_lim
        y_lim = (-np.max(np.abs(self.psi.y)),
                 np.max(np.abs(self.psi.y))) if y_lim is None else y_lim
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)

        # plot the wavefunction lines
        mag_line, = ax.plot(self.psi.x, np.abs(self.psi.y),
                            label="$\sqrt{P}$", linestyle="-")
        real_line, = ax.plot(self.psi.x, self.psi.y.real,
                             label="$Re(\psi)$", linestyle="--")
        imag_line, = ax.plot(self.psi.x, self.psi.y.imag,
                             label="$Im(\psi)$", linestyle="--")

        # initialize the potential fill
        pot_fill = ax.fill_between(
            self.psi.x, 0, self.potential()/700, color='green', alpha=0.2)

        # set labels and legend
        ax.legend(loc="upper right")
        ax.set_xlabel("x")
        ax.set_ylabel("$\psi$")

        # initialize the animation
        def init():
            return mag_line, real_line, imag_line, pot_fill

        # update function for the animation
        def update(frame):
            nonlocal pot_fill
            self.step(dt, **kwargs)
            mag_line.set_ydata(np.abs(self.psi.y))
            real_line.set_ydata(self.psi.y.real)
            imag_line.set_ydata(self.psi.y.imag)

            # update the potential fill if it changes over time
            pot_fill.remove()
            pot_fill = ax.fill_between(
                self.psi.x, 0, self.potential()/700, color='green', alpha=0.2)

            return mag_line, real_line, imag_line, pot_fill

        # create the animation
        ani = anim.FuncAnimation(fig, update, frames=range(
            int(anim_length/dt)), init_func=init, blit=True, interval=1000*dt/anim_speed)

        if notebook:
            return ani
        else:
            plt.show()

##############################################################################################
# Builtin Initial Conditions


def uniform(psi_x, x_area=(-5, 5)):
    cond = (psi_x >= x_area[0]) & (psi_x <= x_area[1])
    psi_y = np.where(cond, 1+0j, 0+0j)/(x_area[1] - x_area[0])
    return psi_y


def plane_wave(psi_x, p=0, h_bar=1.0):
    """
    An eigenstate of the momentum operator and the (zero-potential) energy operator.
    """
    const = 1/np.sqrt(2*np.pi*h_bar)
    psi_y = const*np.exp(1j*p*psi_x/h_bar)
    return psi_y


def wave_packet(psi_x, x_0=0.0, p_0=0.0, sigma_x=0.2, h_bar=1.0):
    const = 1/(2*np.pi*sigma_x**2)**(1/4)
    gauss = np.exp(-((psi_x - x_0)**2)/(4*sigma_x**2))
    momen = np.exp(1j*p_0*psi_x/h_bar)
    psi_y = const*gauss*momen
    return psi_y

##############################################################################################
# Builtin Potentials


def simple_harmonic(t, x, m=1.0, omega=10.0):
    return m*(omega**2)*(x**2)/2


def barrier(t, x, x_0, width=1, height=300):
    cond = (x >= x_0 - width/2) & (x <= x_0 + width/2)
    return np.where(cond, height, 0)

##############################################################################################


if __name__ == "__main__":
    particle = Particle.from_initial(
        wave_packet, x_0=0, p_0=20, sigma_x=0.1).add_potential(barrier, 1.5)
    particle.animate(x_lim=(-4, 4))
