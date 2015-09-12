#!/usr/bin/env python

"""
  Multiscale Turing Patterns using scipy.weave to implement fast box blur in C

  Written by Daniel, not Tripzilch.
"""

from pylab import *
from numpy import *
from numpy.random import random
from scipy.weave import inline, converters

class Filter(object):
  def __init__(self, activator, inhibitor, increment, weight, symmetry):
    super(Filter, self).__init__()
    self.activator = activator
    self.inhibitor = inhibitor
    self.increment = increment
    self.weight = weight
    self.symmetry = symmetry

  def _blur(self, source, end_target, scale):
    target = zeros(shape(source))
    dimx, dimy = shape(source)
    code_template = """
    for (int %(a)s = 0; %(a)s < dim%(a)s; ++%(a)s) {
      double value = 0.0;
      for (int %(b)s = 0; %(b)s <= dim%(b)s + 2*scale; ++%(b)s) {
        value += source(%(aa)s%(b)s %% dim%(b)s%(bb)s);
        if (%(b)s >= 2*scale) {
          target(%(aa)s(%(b)s - scale) %% dim%(b)s%(bb)s) = value / (2*scale+1);
          value -= source(%(aa)s(%(b)s - 2*scale) %% dim%(b)s%(bb)s);
        }
      }
    }
    """
    names = ['dimx', 'dimy', 'source', 'target', 'scale']
    code  = code_template % {"a": "x", "b": "y", "aa": "x, ", "bb": ""}
    inline(code, names, type_converters = converters.blitz)
    source, target = target, end_target
    code += code_template % {"a": "y", "b": "x", "aa": "", "bb": ", y"}
    inline(code, names, type_converters = converters.blitz)


  def apply(self, source):
    dimx, dimy = shape(source)
    activator_field = zeros((dimx, dimy))
    inhibitor_field = zeros((dimx, dimy))
    intermediate1 = zeros((dimx, dimy))
    intermediate2 = zeros((dimx, dimy))
    self._blur(source, intermediate1, self.activator)
    self._blur(intermediate1, intermediate2, self.activator)
    self._blur(intermediate2, activator_field, self.activator)
    self._blur(source, intermediate1, self.inhibitor)
    self._blur(intermediate1, intermediate2, self.inhibitor)
    self._blur(intermediate2, inhibitor_field, self.inhibitor)
    return (activator_field - inhibitor_field)*self.weight

class MSTP(object):
  def __init__(self, dimx, dimy, filters):
    super(MSTP, self).__init__()
    self.dimx = dimx
    self.dimy = dimy
    self.filters = filters
    self.field = random((dimx, dimy))*2.0 -1.0

  def _update_increments(self, v, increment, variation, increment_field):
    dimx, dimy = self.dimx, self.dimy
    code = """
    for (int x = 0; x < dimx; ++x) {
      for (int y = 0; y < dimy; ++y) {
        double var = v(x, y);
        if (std::abs(variation(x, y)) > std::abs(var)) {
          variation(x, y) = var;
          increment_field(x, y) = increment*(var > 0 ? 1.0 : -1.0);
        }
      }
    }
    """
    names = ['dimx', 'dimy', 'variation', 'v', 'increment', 'increment_field']
    inline(code, names, type_converters = converters.blitz)

  def step(self, counts):
    for i in xrange(counts):
      dimx, dimy = self.dimx, self.dimy
      variations = [f.apply(self.field) for f in self.filters]
      variation = ones((dimx, dimy))*1E10
      increment_field = zeros((dimx, dimy))
      for v, f in zip(variations[::-1], self.filters[::-1]):
        self._update_increments(v, f.increment, variation, increment_field)

      self.field += increment_field
      reshaped_field = reshape(self.field, dimx*dimy)
      self.field -= min(*reshaped_field)
      self.field /= 0.5*max(*reshaped_field)
      self.field -= 1

filters = list()
filters.append(Filter(25,  50, 2*0.05, 1, 3))
filters.append(Filter(10,  20, 2*0.04, 1, 3))
filters.append(Filter(5,   10, 2*0.03, 1, 3))
filters.append(Filter(2,   5,  2*0.02, 1, 3))
filters.append(Filter(1,   3,  2*0.01, 1, 3))

mstp = MSTP(300, 300, filters)
f = figure()
ax = f.add_subplot(111)
ion()
while True:
  mstp.step(1)
  ax.clear()
  ax.imshow(mstp.field)
  draw()
