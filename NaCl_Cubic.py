from matplotlib import pyplot as p
from lattice import *

# Set up the crystal structure
lattice = Cubic(5.64)
basis = Basis([('Cl',[0,0,0]),('Cl',[0.5,0.5,0]),
               ('Cl',[0.5,0,0.5]),('Cl',[0,0.5,0.5]),
               ('Na',[0.5,0,0]),('Na',[0,0.5,0]),
               ('Na',[0,0,0.5]),('Na',[0.5,0.5,0.5])],
              l_const=5.64)
crystal = lattice + basis

# Plot a simulated XRD with copper radiation
angles, values = crystal.powder_XRD(1.54)
p.plot(angles, values)

# Add some more info to the plot
p.title(r'Simulated Powder XRD of NaCl, $\lambda = 1.54$')
p.xlabel(r'$2\theta$')
p.ylabel(r'Scattering Intensity per Cubic Angstrom')
p.show()
