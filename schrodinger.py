import numpy as np

##############################################################################################


class Grid1D:
    def __init__(self, domain=(-1.0, 1.0), count=1000, *args, **kwargs):
        self.x = np.linspace(domain[0], domain[1], count)
        self.dx = (self.x[-1] - self.x[0]) / (self.x.size - 1)
        self.y = np.zeros(self.x.size, dtype=complex)

##############################################################################################


class Particle:
    def __init__(self, mass=1.0, h_bar=1.0, *args, **kwargs):
        self.mass = float(mass)
        self.h_bar = float(h_bar)
        self.psi = Grid1D(*args, **kwargs)
        self.potentials = []
        self.time = 0.0

    def add_potential(self, function, period=None, *args, **kwargs):
        pass

    def step(self, dt):
        pass

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
