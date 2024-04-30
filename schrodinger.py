import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

##############################################################################################


class Grid1D:
    def __init__(self, domain=(-20.0, 20.0), dx=20/999, **kwargs):
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
        return particle

    def add_initial(self, func, *args, **kwargs):
        self.psi.y += func(self.psi.x, *args, **kwargs)
        return self

    def add_potential(self, func, *args, **kwargs):
        def v(t, x): return func(t, x, *args, **kwargs)
        self.potentials.append(v)
        return self

    def potential(self):
        V = np.zeros_like(self.psi.x)
        for pot in self.potentials:
            V += pot(self.time.real, self.psi.x)
        return V

    def hamiltonian(self):
        mom = -(self.h_bar**2)*self.psi.sec_deriv()/(2*self.mass)
        pot = self.potential()*self.psi.y
        return mom + pot

    def normalize(self):
        norm = np.sqrt(np.sum(np.abs(self.psi.y)**2 * self.psi.dx))
        self.psi.y /= norm
        return self

    def step(self, dt, before_step=None, damp_eigen=False, **kwargs):
        dt = -1j*abs(dt) if damp_eigen else abs(dt)
        if self.running:
            def dydt(psi):
                momen = 1j*self.h_bar*psi.sec_deriv()/(2*self.mass)
                pot = 1j*self.potential()*psi.y/self.h_bar
                return momen - pot

            self.normalize()
            before_step(self) if before_step is not None else None
            self.psi.rk4_step(dydt, dt, **kwargs)
            self.time += dt

    def find_eigen(self, count=4, dt=1/3600, damp_time=2.0, **kwargs):
        self.eigenstates = [None]*count
        self.eigenvalues = [None]*count
        self.time = 0.0
        initial_psi_y = np.copy(self.psi.y)
        for n in range(count):
            def remove_prev(particle):
                for i in range(n):
                    overlap = np.sum(np.conjugate(
                        particle.eigenstates[i])*particle.psi.y*particle.psi.dx)
                    particle.psi.y -= overlap*particle.eigenstates[i]
            while abs(self.time) <= damp_time:
                self.step(dt, before_step=remove_prev,
                          damp_eigen=True, **kwargs)
            self.eigenstates[n] = np.copy(self.psi.y)
            self.eigenvalues[n] = np.sum(np.conjugate(
                self.eigenstates[n])*self.hamiltonian()*self.psi.dx).real
            self.psi.y = initial_psi_y
            self.time = 0.0

    def pause_play(self):
        self.running = not self.running

    def animate(self, dt=1/3600, anim_length=5, anim_speed=1, x_lim=None, y_lim=None, notebook=False, v_scale=700, title="Schrodinger Simulation", **kwargs):
        fig = plt.figure(title)
        ax = fig.add_subplot()
        plt.title(title)

        # set the limits for the x and y axes
        x_lim = (self.psi.x[0], self.psi.x[-1]) if x_lim is None else x_lim
        y_lim = (-np.max(np.abs(self.psi.y)),
                 np.max(np.abs(self.psi.y))) if y_lim is None else y_lim
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)

        # plot the wavefunction lines
        mag_line, = ax.plot(self.psi.x, np.abs(self.psi.y),
                            label="$\\sqrt{P}$", linestyle="-")
        real_line, = ax.plot(self.psi.x, self.psi.y.real,
                             label="$Re(\\psi)$", linestyle="--")
        imag_line, = ax.plot(self.psi.x, self.psi.y.imag,
                             label="$Im(\\psi)$", linestyle="--")
        pot_line, = ax.plot(self.psi.x, self.potential() /
                            v_scale, label="V", color="black")

        # initialize the potential fill
        pot_fill = ax.fill_between(self.psi.x, min(0, np.min(self.potential(
        )/v_scale)), self.potential()/v_scale, color="green", alpha=0.2)

        # show the time in the top left
        time_text = ax.text(0.02, 0.97, "", transform=ax.transAxes, fontsize=10,
                            verticalalignment="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

        # set labels and legend
        ax.legend(loc="upper right")
        ax.set_xlabel("x")
        ax.set_ylabel("$\\psi$")

        # initialize the animation
        def init():
            return mag_line, real_line, imag_line, pot_line, pot_fill, time_text

        # update function for the animation
        def update(frame):
            nonlocal pot_fill
            nonlocal time_text
            self.step(dt, **kwargs)
            mag_line.set_ydata(np.abs(self.psi.y))
            real_line.set_ydata(self.psi.y.real)
            imag_line.set_ydata(self.psi.y.imag)
            pot_line.set_ydata(self.potential()/v_scale)
            pot_fill.remove()
            pot_fill = ax.fill_between(self.psi.x, min(0, np.min(self.potential(
            )/v_scale)), self.potential()/v_scale, color="green", alpha=0.2)
            time_text.set_text("Time: {:.2f}".format(self.time))
            return mag_line, real_line, imag_line, pot_line, pot_fill, time_text

        if dt != 0.0:
            frame_count = int(anim_length/dt)
        else:
            frame_count = 0

        # create the animation
        ani = anim.FuncAnimation(fig, update, frames=range(
            frame_count), init_func=init, blit=True, interval=1000*dt/anim_speed)

        if notebook:
            return ani
        else:
            plt.show()

    def plot_eigen(self, x_lim=None, y_lim=None, eigenvalues=True, imag=False, prob=False, v_scale=1, title="Potential & Eigenstates"):
        fig = plt.figure(title)
        ax = fig.add_subplot()
        plt.title(title)
        x_lim = [-2, 2] if x_lim is None else x_lim
        y_lim = [-1.5, 1.5*self.eigenvalues[-1] -
                 0.5*self.eigenvalues[-2]+1.5] if y_lim is None else y_lim
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)
        for n, psi_y in reversed(list(enumerate(self.eigenstates))):
            ax.plot(self.psi.x, psi_y.real +
                    self.eigenvalues[n], label="E = {:.2f}".format(self.eigenvalues[n]))
            if imag:
                ax.plot(self.psi.x, psi_y.imag +
                        self.eigenvalues[n], linestyle="--", color="red")
            if prob:
                ax.plot(self.psi.x, np.abs(psi_y) +
                        self.eigenvalues[n], linestyle="--", color="black")
        ax.plot(self.psi.x, self.potential()/v_scale, label="V", color="black")
        ax.fill_between(self.psi.x, min(0, np.min(self.potential()/v_scale)),
                        self.potential()/v_scale, color="green", alpha=0.2)
        ax.set_xlabel("x")
        ax.set_ylabel("$\\psi, E$")
        if eigenvalues:
            ax.legend(loc="upper right")
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


def simple_harmonic(t, x, m=1.0, omega=5.0):
    return m*(omega**2)*(x**2)/2


def barrier(t, x, x_0, width=1, height=300):
    cond = (x >= x_0 - width/2) & (x <= x_0 + width/2)
    return np.where(cond, height, 0)

def coulomb(t, x, x_0=0.0, q_squared=-1.0, epsilon=0.001):
    return q_squared/(4*np.pi*epsilon*np.abs(x-x_0))

##############################################################################################


if __name__ == "__main__":
    particle = Particle.from_initial(
        wave_packet, x_0=0, p_0=20, sigma_x=0.1).add_potential(barrier, 1.5)
    particle.animate(x_lim=(-4, 4))
