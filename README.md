# Pylattice

## What it Does

pylattice calculates x-ray diffraction spectra for most crystalline materials. You tell it the crystal structure and the wavelength of light, and it does the rest. It's designed to be used in other scripts that need to programatically calculate x-ray diffraction spectra.

## Dependencies

pylattice only depends on [numpy](http://www.numpy.org). If you want to run the examples which plot spectra, you will need [matplotlib](http://matplotlib.org) as well.

## Installation

Simply put the "lattice" folder somewhere on your pythonpath.

## How to Use It

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

__powder_xrd__ returns a dictionary mapping angles to intensities. With the keyword argument "get_mults=True", it will also return a dictionary mapping angles to multiplicities

__spectrumify__ turns the dictionary of scattering data into a fake scattering spectrum that can easily understood by humans.
