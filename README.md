# Pylattice

## What it Does

Pylattice calculates powder x-ray diffraction spectra of crystalline materials. You define the crystal structure and the wavelength of light, and it does the rest. It's designed to easily mesh with python scripts and programs that need to calculate powder XRD spectra.

## Dependencies

Pylattice only depends on [numpy](http://www.numpy.org). If you want to plot the calculated spectra, you will need [matplotlib](http://matplotlib.org) as well.

## Installation

Simply put the "lattice" folder somewhere on your pythonpath.

## Basic Use

Pylattice has simple syntax and just tries to get the job done. Below is an example that calculates a powder XRD spectrum for NaCl:

```python
from matplotlib import pyplot as p
from lattice import *

# Set up the crystal structure
lattice = FCC(5.63)
basis = Basis([('Cl',[0,0,0]),
               ('Na',[0.5,0.5,0.5])],
              l_const=5.63)
crystal = lattice + basis

# Plot a simulated XRD with copper radiation
scattering_data = powder_XRD(crystal, 1.5405)
angles, values = spectrumify(scattering_data)
p.plot(angles, values)

# Add some more info to the plot
p.title(r'Simulated Powder XRD of NaCl, $\lambda = 1.5405$')
p.xlabel(r'$2\theta$')
p.ylabel(r'Scattering Intensity per Cubic Angstrom')
p.show()
```

#### Defining a Crystal

* Define a lattice, either manually with __Lattice__  or using one of the helper functions __FCC__,__BCC__,__Cubic__,or __Hexagonal__.
* Define a basis using __Basis__, providing a list of atoms and locations
* Define the crystal as the sum of the lattice and basis

#### Calculating the Spectrum

* __powder_xrd__(wavelength) returns a dictionary mapping angles to intensities. With the keyword argument "get_mults=True", it will also return a dictionary mapping angles to multiplicities

* __spectrumify__ turns the dictionary of scattering data into a fake scattering spectrum that can easily understood by humans.
