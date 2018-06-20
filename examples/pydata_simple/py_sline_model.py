#!/usr/bin/env python

import math

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
# New ModelBuilder object, defines a new model NcHIPrimExample 
# implementing the base Ncm.Model abstract class.
# 
mb = Ncm.ModelBuilder.new (Ncm.Model, "NcPySLineModel", 
                           "A simple python example model")

#
# New parameter m to describe the slope 
# Allowed interval: [0, 5]
# Default scale: 0.1
# Absolute tolerance: 0
# Default value: 2
#
mb.add_sparam ("m", "m", 0.0, 5.0, 0.1, 0.0, 2.0, Ncm.ParamType.FREE)

#
# New parameter b to describe the intercept
# Allowed interval: [-10, 10]
# Default scale: 0.1
# Absolute tolerance: 0
# Default value: 1
#
mb.add_sparam ("b", "b", -10.0, 10.0, 0.1, 0.0, 1.0, Ncm.ParamType.FREE)

#
# Creates a new GObject, it is not a Python object yet!
#
GNcPySLineModel = mb.create ()

#
# Register the new object in the GObject type system by creating a new instance
#
GObject.new (GNcPySLineModel)

#
# Gets the Python version of the object (.pytype) and register
# it as a PyGObject object.
#
NcPySLineModel = GNcPySLineModel.pytype
GObject.type_register (NcPySLineModel)

#
# Creating a new class implementing our object NcPySLineModel
#
class PySLineModel (NcPySLineModel):
  #
  # Defining some property which is not part of the model paramers.
  # All model parameters must be defined by the ModelBuilder.
  #
  some_property = GObject.Property (type = str)

  #
  # Calling the father's constructor
  #  
  def __init__ (self):
    NcPySLineModel.__init__ (self)

  #
  # Method to calculate the y(x)
  #
  def f_x (self, x):
    return math.exp (self.props.m * x) * self.props.b

#
# Register our new Python class PyNcPySLineModel
#
GObject.type_register (PySLineModel)
