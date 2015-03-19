from matplotlib import pyplot as p
from matplotlib import rc
rc('text', usetex=True)
from lattice import *

# Set up the crystal structure
lattice = Cubic(5.64)
basis = Basis([('Cl',[0,0,0]),('Cl',[0.5,0.5,0]),
               ('Cl',[0.5,0,0.5]),('Cl',[0,0.5,0.5]),
               ('N',[0.5,0,0]),('N',[0,0.5,0]),
               ('N',[0,0,0.5]),('N',[0.5,0.5,0.5])],
              l_const=5.64)
crystal = lattice + basis

# Plot a simulated XRD with copper radiation
p.plot(*crystal.powder_XRD(1.54))

# Add some more info to the plot
p.title(r'Simulated Powder XRD of NaCl, $\lambda = 1.54$')
p.xlabel(r'$2\theta$')
p.ylabel(r'Scattering Intensity per Cubic Angstrom')
p.show()
