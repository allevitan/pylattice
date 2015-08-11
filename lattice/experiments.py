from __future__ import division, print_function
import numpy as n
import itertools as it
from time import time


def find_accessible_rlvs(crystal, wavelength):
    """
    Generates a list of accessible reciprocal lattice
    vectors. To be accessible, the magnitude of a rlv's
    wavevector must be less than twice that of the input
    radiation's wavenumber."""
    
    # The wavenumber of the input wavelength
    nu = 2*n.pi/wavelength
    
    # Now we find the shortest distance to a wall of a 
    # parallelogram "shell" in the reciprocal lattice
    min_step = min(abs(n.dot(
        (crystal.rlattice[0]+crystal.rlattice[1]
         +crystal.rlattice[2]),
        n.cross(crystal.rlattice[i],crystal.rlattice[j])
        /n.linalg.norm(n.cross(crystal.rlattice[i],crystal.rlattice[j]))))
                   for i,j in [(0,1),(1,2),(2,0)])
    
    # If we look at all the points in this many parallelogram
    # "shells", we can't miss all the accessible wavevectors
    num_shells = int(2*nu / min_step)
    
    # Now we generate these possibilities
    possibilities = [(crystal.rlattice[0]*h + crystal.rlattice[1]*j
                      + crystal.rlattice[2]*k)
                     for h,j,k in it.product(
                             range(-num_shells,num_shells+1),
                             repeat=3)]
    
    # And we filter the possibilities, getting rid of all the
    # rlvs that are too long and the 0 vector
    rlvs = [rlv for rlv in possibilities 
            if n.linalg.norm(rlv) < 2*nu
            and not n.allclose(rlv,0)]
    
    return n.array(rlvs)


def powder_XRD(crystal,wavelength, get_mults=False):
    """
    Generates a powder XRD spectrum for radiation with the
    given wavelength (in angstroms)
    """
    
    # The wavenumber of the input wavelength
    nu = 2*n.pi/wavelength

    # Make a list of the accessible rlvs
    rlvs = find_accessible_rlvs(crystal,wavelength)
    
    # Now we calculate the scattering intensity from each rlv
    intensities = {
        tuple(rlv): n.abs(crystal.structure_factor(rlv))**2
        for rlv in rlvs}
    
    # Now sum up all rlvs with the same magnitude. We also
    # get rid of all the scattering vectors with 0 intensity
    magnitudes = {}
    multiplicities = {}
    for rlv, intensity in intensities.items():
        repeat = False
        mag = n.linalg.norm(rlv)
        for oldmag in magnitudes:
            if n.isclose(mag,oldmag):
                magnitudes[oldmag] += intensity
                multiplicities[oldmag] += 1
                repeat = True
                break
        if not repeat and not n.isclose(mag,0):
            multiplicities[mag] = 1
            magnitudes[mag] = intensity
        
    # Now we reformat the multiplicity data in a nice way
    multiplicities = {2 * n.arcsin(mag / (2 * nu)) * 180 / n.pi:
                      multiplicity
                      for mag, multiplicity in multiplicities.items()
                      if not n.allclose(magnitudes[mag],0)}

    # And now we calculate the scattering intensities
    # (a.u. per steradian) as a function of scattering angle
    intensities = {2 * n.arcsin(mag / (2 * nu)) * 180 / n.pi:
                   intensity * 
                   # This factor corrects for the fact that the same total
                   # power in the debye scherrer rings is more
                   # concentrated when 2\theta is near 0 or 2pi
                   1 / n.sin(2*n.arcsin(mag/(2*nu))) *
                   # This factor corrects for the probability that any
                   # given crystal domain will scatter into the rlv
                   1 / mag *
                   # This factor corrects for polarization effects,
                   # Assuming an unpolarized input beam and no polarization
                   # analysis
                   (1 + n.cos(2*n.arcsin(mag/(2*nu)))**2)/2
                   for mag, intensity in magnitudes.items()
                   if not n.allclose(intensity,0)}
    if get_mults:
        return intensities, multiplicities
    else:
        return intensities



def spectrumify(scattering_data, instr_broadening=0.1):
    """
    This is just a nice function to turn the raw scattering data
    into a human-readable scattering spectrum
    """
    graph_angles = n.linspace(0,180,10000)
    graph_intensities = n.zeros(graph_angles.shape)
    
    for angle, intensity in sorted(scattering_data.items()):
        graph_intensities += intensity * \
                             n.exp(-(graph_angles - angle)**2 / \
                                   (2*(instr_broadening)**2))
        
    return graph_angles, graph_intensities
