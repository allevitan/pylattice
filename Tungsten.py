from matplotlib import pyplot as p
from lattice import *

# Set up the crystal structure
lattice = BCC(3.155)
basis = Basis([('W',[0,0,0])],l_const=3.155)
crystal = lattice + basis

# Plot a simulated XRD with copper radiation
angles, values = crystal.powder_XRD(1.5405)
p.plot(angles, values)

# Add some more info to the plot
p.title(r'Simulated Powder XRD of Tungsten, $\lambda = 1.5405$')
p.xlabel(r'$2\theta$')
p.ylabel(r'Scattering Intensity per Cubic Angstrom')
p.show()
