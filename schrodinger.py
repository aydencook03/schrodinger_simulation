import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

##############################################################################################


class Grid1D:
    def __init__(self, domain=(-10.0, 10.0), dx=20/999):
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

    def rk4_step(self, dydt, dt, bound="zero", *args, **kwargs):
        k1 = dydt(self, *args, **kwargs)
        k2 = dydt(self + k1*dt/2, *args, **kwargs)
        k3 = dydt(self + k2*dt/2, *args, **kwargs)
        k4 = dydt(self + k3*dt, *args, **kwargs)
        self.y += (k1 + 2*k2 + 2*k3 + k4)*dt/6
        if bound == "zero":
            self.y[0] = 0
            self.y[-1] = 0
        elif bound == "periodic":
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

    def wave_packet(x_0, p_0, sigma_x, **kwargs):
        particle = Particle(**kwargs)
        const = 1/(2*np.pi*sigma_x**2)**(1/4)
        gauss = np.exp(-((particle.psi.x - x_0)**2)/(4*sigma_x**2))
        momen = np.exp(1j*p_0*particle.psi.x/particle.h_bar)
        particle.psi.y = const*gauss*momen
        return particle

    def add_potential(self, func, period=None, *args, **kwargs):
        def v(t, x): return func(t, x, *args, **kwargs)
        self.potentials.append(v)

    def normalize(self):
        norm = np.sqrt(np.sum(np.abs(self.psi.y)**2 * self.psi.dx))
        self.psi.y /= norm

    def step(self, dt, **kwargs):
        V = np.zeros_like(self.psi.x)
        for pot in self.potentials:
            V += pot(self.time, self.psi.x)

        def dydt(psi):
            return 1j*self.h_bar*psi.sec_deriv()/(2*self.mass) - 1j*V*psi.y/self.h_bar

        self.psi.rk4_step(dydt, dt, **kwargs)
        self.normalize()
        self.time += dt

    def imaginary_step(self, dt):
        pass

    def animate(self, dt=1/3600, anim_length=5, anim_speed=1, x_lim=None, y_lim=None, **kwargs):
        fig = plt.figure("Schrodinger Simulation")
        ax = fig.add_subplot()

        mag_line, = ax.plot(self.psi.x, np.abs(self.psi.y),
                            label="$\sqrt{P}$", linestyle="-")
        real_line, = ax.plot(self.psi.x, self.psi.y.real,
                             label="$Re(\psi)$", linestyle="--")
        imag_line, = ax.plot(self.psi.x, self.psi.y.imag,
                             label="$Im(\psi)$", linestyle="--")

        ax.legend(loc="upper right")
        ax.set_xlabel("x")
        ax.set_ylabel("$\psi$")

        x_lim = (self.psi.x[0], self.psi.x[-1]) if x_lim is None else x_lim
        y_lim = (-np.max(np.abs(self.psi.y)),
                 np.max(np.abs(self.psi.y))) if y_lim is None else y_lim

        def init():
            ax.set_xlim(*x_lim)
            ax.set_ylim(*y_lim)
            return mag_line, real_line, imag_line

        def update(frame):
            self.step(dt, **kwargs)
            mag_line.set_ydata(np.abs(self.psi.y))
            real_line.set_ydata(self.psi.y.real)
            imag_line.set_ydata(self.psi.y.imag)
            return mag_line, real_line, imag_line

        ani = anim.FuncAnimation(fig, update, frames=range(
            int(anim_length/dt)), init_func=init, blit=True, interval=int(1000*dt/anim_speed))
        plt.show()

##############################################################################################


if __name__ == "__main__":
    particle = Particle.wave_packet(x_0=0, p_0=10, sigma_x=0.1)
    particle.animate(x_lim=(-1, 1))
