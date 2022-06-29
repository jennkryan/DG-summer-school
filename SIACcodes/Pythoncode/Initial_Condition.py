# -*- coding: utf-8 -*-
"""
@author: Laura C. Kühle

"""
import numpy as np


class InitialCondition(object):
    def __init__(self, left_bound, right_bound, config):
        self._left_bound = left_bound
        self._right_bound = right_bound

        self._reset(config)

    def _reset(self, config):
        self._interval_len = self._right_bound-self._left_bound

    def get_name(self):
        return self.__class__.__name__

    def induce_adjustment(self, value):
        pass

    def randomize(self, config):
        pass

    def calculate(self, x):
        while x < self._left_bound:
            x = x + self._interval_len
        while x > self._right_bound:
            x = x - self._interval_len
        return self._get_point(x)

    def _get_point(self, x):
        pass


class Sine(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
<<<<<<< HEAD
        self._factor = config.pop('factor', 1)
=======
        self._factor = config.pop('factor', 2)
>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e
        self._stretch_factor = config.pop('stretch_factor', 1)
        self._height_adjustment = config.pop('height_adjustment', 0)

    def randomize(self, config):
        factor = config.pop('factor', np.random.uniform(low=-100, high=100))
        config = {'factor': factor}
        self._reset(config)

    def _get_point(self, x):
<<<<<<< HEAD
        # currently doesn't always use the correct self._factor
=======
>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e
        return self._height_adjustment + self._stretch_factor * np.sin(self._factor * np.pi * x)


class Box(InitialCondition):
    def _get_point(self, x):
        if x < -1:
            x = x + 2
        if x > 1:
            x = x - 2
        if (x >= -0.5) & (x <= 0.5):
            return 1
        else:
            return 0

    def is_smooth(self):
        return False


class FourPeakWave(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Set additional necessary parameter
        self._alpha = 10
        self._delta = 0.005
        self._beta = np.log(2) / (36 * self._delta**2)
        self._a = 0.5
        self._z = -0.7

    def is_smooth(self):
        return False

    def _get_point(self, x):
        if (x >= -0.8) & (x <= -0.6):
            return 1/6 * (self._gaussian_function(x, self._z-self._delta)
                          + self._gaussian_function(x, self._z+self._delta)
                          + 4 * self._gaussian_function(x, self._z))
        if (x >= -0.4) & (x <= -0.2):
            return 1
        if (x >= 0) & (x <= 0.2):
            return 1 - abs(10 * (x-0.1))
        if (x >= 0.4) & (x <= 0.6):
            return 1/6 * (self._elliptic_function(x, self._a-self._delta)
                          + self._elliptic_function(x, self._a+self._delta)
                          + 4 * self._elliptic_function(x, self._a))
        return 0

    def _gaussian_function(self, x, z):
        return np.exp(-self._beta * (x-z)**2)

    def _elliptic_function(self, x, a):
        return np.sqrt(max(1 - self._alpha**2 * (x-a)**2, 0))


class Linear(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
        self._factor = config.pop('factor', 1)

    def randomize(self, config):
        factor = config.pop('factor', np.random.uniform(low=-100, high=100))
        config = {'factor': factor}
        self._reset(config)

    def _get_point(self, x):
        return self._factor * x


class LinearAbsolut(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
        self._factor = config.pop('factor', 1)

    def is_smooth(self):
        return False

    def randomize(self, config):
        factor = config.pop('factor', np.random.uniform(low=-100, high=100))
        config = {'factor': factor}
        self._reset(config)

    def _get_point(self, x):
        return self._factor * abs(x)


class DiscontinuousConstant(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
        self._x0 = config.pop('x0', 0)
        self._left_factor = config.pop('left_factor', 1)
        self._right_factor = config.pop('right_factor', 0.5)

    def _get_point(self, x):
        return self._left_factor * (x <= self._x0) + self._right_factor * (x > self._x0)


class Polynomial(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
        self._factor = config.pop('factor', 1)
        self._exponential = config.pop('exponential', 2)

    def randomize(self, config):
        factor = config.pop('factor', np.random.uniform(low=-100, high=100))
        exponential = config.pop('exponential', np.random.randint(2, high=6))
        config = {'factor': factor, 'exponential': exponential}
        self._reset(config)

    def _get_point(self, x):
        return self._factor * (x ** self._exponential)


class Continuous(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
        self._factor = config.pop('factor', 1)

    def randomize(self, config):
        factor = config.pop('factor', np.random.uniform(low=-100, high=100))
        config = {'factor': factor}
        self._reset(config)

    def _get_point(self, x):
        return self._factor


class HeavisideOneSided(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
        self._factor = config.pop('factor', -1)

    def is_smooth(self):
        return False

    def randomize(self, config):
        factor = config.pop('factor', np.random.choice([-1, 1]))
        config = {'factor': factor}
        self._reset(config)

    def _get_point(self, x):
        return self._factor - 2 * self._factor * np.heaviside(x, 0)


class HeavisideTwoSided(InitialCondition):
    def _reset(self, config):
        super()._reset(config)

        # Unpack necessary configurations
        self._left_factor = config.pop('left_factor', 1)
        self._right_factor = config.pop('right_factor', 2)
        self._adjustment = config.pop('adjustment', 0)

    def is_smooth(self):
        return False

    def induce_adjustment(self, value):
        self._adjustment = value

    def randomize(self, config):
        left_factor = config.pop('left_factor', np.random.choice([-1, 1]))
        right_factor = config.pop('right_factor', np.random.choice([-1, 1]))
        adjustment = config.pop('adjustment', np.random.uniform(low=-1, high=1))
        config = {'left_factor': left_factor, 'right_factor': right_factor, 'adjustment': adjustment}
        self._reset(config)

    def _get_point(self, x):
        return self._left_factor\
               - self._left_factor * np.heaviside(x - self._adjustment, 0)\
               - self._right_factor * np.heaviside(x + self._adjustment, 0)
