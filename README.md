# Pylattice

## Why it Exists

I wanted to learn how X-Ray Diffraction works, and there's no better way to do that than by simulating it.

## How to Use It

Let's say for some reason you'd rather use this than one of the many infinitely better, more extensible, and more usable tools that other people have made. Or maybe you're me in the future trying to figure out how the hell this ever worked. Anyway, here's a simple example (for Sodium Chloride)

```python
from matplotlib import pyplot as p
from lattice import *

# Set up the crystal lattice, in this case using a handy-dandy
# FCC subclass of Lattice(). The lattice constant should be in
# Angstroms
lattice = FCC(5.64)

# Now we set up our atomic basis. Specify atoms and their
# corresponding locations, along with a lattice constant
# (in Anstroms) to multiply all distances by
basis = Basis([('Cl',[0,0,0]),
               ('N',[0.5,0.5,0.5])],
              l_const=5.64)

# Yes, I know it's really a lattice convolved with a basis,
# but this will make condensed matter professors everywhere
# happy
crystal = lattice + basis

# powder_XRD outputs a fancy-lookin' simulated spectrum. Don't
# trust the intensities, but the angles should be spot-on
p.plot(*crystal.powder_XRD(1.54))
```
