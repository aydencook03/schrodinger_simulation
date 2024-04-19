import numpy as np

##############################################################################################


class Grid1D:
    def __init__(self, domain=(-1.0, 1.0), count=1000):
        self.x = np.linspace(domain[0], domain[1], count)
        self.dx = (domain[1] - domain[0]) / (count - 1)
        self.y = np.zeros(count, dtype=complex)

    def __add__(self, other):
        result = Grid1D(domain=(self.x[0], self.x[-1]), count=self.x.size)
        result.y = self.y + other.y
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

    def rk4_step(self, dydt, dt, *args, **kwargs):
        k1 = dydt(self, *args, **kwargs)
        k2 = dydt(self + k1*dt/2, *args, **kwargs)
        k3 = dydt(self + k2*dt/2, *args, **kwargs)
        k4 = dydt(self + k3*dt, *args, **kwargs)
        self.y += (k1 + 2*k2 + 2*k3 + k4)*dt/6

##############################################################################################


class Particle:
    def __init__(self, mass=1.0, h_bar=1.0, **kwargs):
        self.mass = float(mass)
        self.h_bar = float(h_bar)
        self.psi = Grid1D(**kwargs)
        self.potentials = []
        self.time = 0.0

    def wave_packet():
        pass

    def add_potential(self, func, period=None, *args, **kwargs):
        def v(t, x): return func(t, x, *args, **kwargs)
        self.potentials.append(v)

    def normalize(self):
        norm = np.sqrt(np.sum(np.abs(self.psi.y)**2 * self.psi.dx))
        self.psi.y /= norm

    def step(self, dt):
        V = np.zeros_like(self.psi.x)
        for pot in self.potentials:
            V += pot(self.time, self.psi.x)

        def dydt(psi):
            return 1j*self.h_bar*psi.sec_deriv()/(2*self.mass) - 1j*V*psi.y/self.h_bar

        self.psi.rk4_step(dydt, dt)
        self.normalize()
        self.time += dt

    def imaginary_step(self, dt):
        pass

##############################################################################################


class Grid2D:
    def __init__(self):
        pass

##############################################################################################


class System:
    def __init__(self, mass_1=1.0, mass_2=1.0, stats="fermion", h_bar=1.0, *args, **kwargs):
        self.mass_1 = float(mass_1)
        self.mass_2 = float(mass_2)
        self.h_bar = float(h_bar)
        self.stats = stats
        self.psi = Grid2D(*args, **kwargs)
        self.potentials = []
        self.time = 0.0

    def add_potential(self, function, *args, **kwargs):
        pass

    def step(self, dt):
        pass

    def imaginary_step(self, dt):
        pass

##############################################################################################
