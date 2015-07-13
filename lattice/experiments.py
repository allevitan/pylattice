from __future__ import division, print_function
import numpy as n
import itertools as it

def powder_XRD(crystal,wavelength):
    """
    Generates a powder XRD spectrum for radiation with the
    given wavelength (in angstroms)
    """
    # We generate a list of accessible reciprocal lattice
    # vectors. To be accessible, the magnitude of a rlv's
    # wavevector must be less than twice that of the input
    # radiation's wavenumber.
    
    #The input wavenumber.
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
    
    # Now we renormalize the intensities to account for the fact that
    # the same lattice can be described by different unit cells
    unit_vol = n.abs(n.dot(crystal.lattice[0],n.cross(
        crystal.lattice[1],crystal.lattice[2])))
    
    # Now we calculate the scattering intensity from each rlv
    intensities = {
        tuple(rlv): n.abs(crystal.structure_factor(rlv)/unit_vol)**2
        for rlv in rlvs}
    
    # We actually only care about the magnitudes of the rlvs
    magnitudes = {}
    for rlv, intensity in intensities.items():
        repeat = False
        mag = n.linalg.norm(rlv)
        for oldmag in magnitudes:
            if n.isclose(mag,oldmag):
                magnitudes[oldmag] += intensity
                repeat = True
                break
        if not repeat:
            magnitudes[mag] = intensity
        
    # Now we calculate the scattering angles and intensities
    angles = {2 * n.arcsin(mag / (2 * nu)) * 180 / n.pi:
              intensity * 
              # This factor corrects for the fact that the same total
              # power in the debye scherrer rings is more
              # concentrated at 2\theta =  0 and 2pi
              1 / n.sin(2*n.arcsin(mag/(2*nu))) *
              # This factor corrects for the angular dependence of
              # scattering probability given an equal incident
              # scattering wavevector and an equal alowed variance
              # around the scattering vector
              1 * #cos(theta)/cos(theta)
              # This factor corrects for the fact that destructive
              # interference builds up faster at shorter wavelengths,
              # meaning that the allowed variance around the scattering
              # vector is proportional to 1/mag of scattering vector
              1 / mag *
              # This factor corrects for polarization effects,
              # Assuming an unpolarized input beam and no polarization
              # analysis
              (1 + n.cos(2*n.arcsin(mag/(2*nu)))**2)/2
              for mag, intensity in magnitudes.items()
              if not n.allclose(intensity,0)}
    
    return angles
    

def spectrumify(scattering_data):
    """
    This is just a nice function to turn the raw scattering data
    into a human-readable approximation of a scattering spectrum
    """
    graph_angles = n.linspace(0,180,5000)
    graph_intensities = n.zeros(graph_angles.shape)
    
    max_peak = n.max(scattering_data.values())
    
    for angle, intensity in sorted(scattering_data.items()):
        graph_intensities += intensity * \
                             n.exp(-(graph_angles - angle)**2 / (2*(0.1)**2))
        
    return graph_angles, graph_intensities


    
