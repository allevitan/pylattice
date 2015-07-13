from matplotlib import pyplot as p
from lattice import *

# Set up the crystal structure
lattice = Cubic(5.63)
basis = Basis([('Cl',[0,0,0]),('Cl',[0.5,0.5,0]),
               ('Cl',[0.5,0,0.5]),('Cl',[0,0.5,0.5]),
               ('Na',[0.5,0,0]),('Na',[0,0.5,0]),
               ('Na',[0,0,0.5]),('Na',[0.5,0.5,0.5])],
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
