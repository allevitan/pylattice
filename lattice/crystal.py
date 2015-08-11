from __future__ import print_function, division
import numpy as n
from os import path


def get_form_factors():
    """
    Returns a dictionary containing functions to calculate the
    form factors of common atoms and ions. The constants were
    sourced from the page at
    http://lamp.tu-graz.ac.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    """
    form_factors_file = '/'.join(path.realpath(__file__).split('/')[:-1]) \
                       + '/form_factors.csv'
    constants = list(n.loadtxt(form_factors_file,
        skiprows=1,
        delimiter=',',
        usecols={1,2,3,4,5,6,7,8,9}
    ))
    labels = list(n.loadtxt(form_factors_file,
        skiprows=1,
        delimiter=',',
        usecols={0},
        dtype='a'
    ))
    
    #This lambda function has stupid syntax because python
    #doesn't close over the variables until the function is
    #called, unless we include them as keyword arguments.
    form_factors = [lambda q, a1=a1,b1=b1,a2=a2,b2=b2, \
                    a3=a3,b3=b3,a4=a4,b4=b4,c=c:\
                    a1 * n.exp(-b1*(q/(4*n.pi))**2) + \
                    a2 * n.exp(-b2*(q/(4*n.pi))**2) + \
                    a3 * n.exp(-b3*(q/(4*n.pi))**2) + \
                    a4 * n.exp(-b4*(q/(4*n.pi))**2) + c
                    for a1,b1,a2,b2,a3,b3,a4,b4,c in constants]
    form_factors = {label:form_factor
            for label, form_factor
            in zip(labels, form_factors)}
    return form_factors
        

class Lattice(object):
    """
    Stores a crystal lattice, defined by 3 lattice vectors in
    real space. The vectors should be in units of Angstrom."""
    
    def __init__(self,a1,a2,a3):
        self.lattice = [n.array(a1),n.array(a2),n.array(a3)]
        self.rlattice = [
            2*n.pi*n.cross(a2,a3) /
            n.dot(a1, n.cross(a2,a3)),
            2*n.pi*n.cross(a3,a1) /
            n.dot(a2, n.cross(a3,a1)),
            2*n.pi*n.cross(a1,a2) /
            n.dot(a3, n.cross(a1,a2))
        ]

    def __add__(self,basis):
        if 'Basis' not in str(type(basis)):
            raise TypeError('A Lattice can only be added to a Basis')
        return Crystal(self,basis)



class Basis(object):
    """
    Stores a basis, defined by atomic names and sites like so:
    basis = Basis(('C',[0,0,0]),('C',[0.25,0.25,0.25]),l_const=l_const)
    Atomic sites are specified in cartesian space (not in the lattice
    basis) and in units of l_const, which is 1 Angstrom by default.
    """
    
    def __init__(self,atoms,l_const=1):
        #first_atom is included to ensure something is in the basis
        self.basis = [(atom, l_const*n.array(site)) for atom, site in atoms]


class Crystal(object):
    """
    Stores a lattice and a basis, and automatically calculates
    the reciprocal lattice and the structure factor.
    The calculated reciprocal lattice will include some vectors
    that do not scatter because of symmetry - the structure
    factor at those vectors will be 0.
    """
    
    def __init__(self,lattice,basis):
        self.lattice = lattice.lattice
        self.rlattice = lattice.rlattice
        self.basis = basis.basis
        self.structure_factor = self.gen_structure_factor()
    
    
    def gen_structure_factor(self):
        """
        Returns a function that calculates the structure factor
        per cubic Angstrom"""

        form_factors = get_form_factors()
        if not all(atom in form_factors for atom, site in self.basis):
            raise KeyError('Specified atom has no form factor in database')
        
        unit_vol = n.abs(n.dot(self.lattice[0],n.cross(
            self.lattice[1],self.lattice[2])))

        def structure_factor(q):
            return sum([n.exp(1j*n.dot(q,site)) * \
                            form_factors[atom](n.linalg.norm(q))
                            for atom, site in self.basis]) / unit_vol

        return structure_factor

    
    
#
# This is the section of the code where some nice pre-made lattice
# classes are defined.
#

class FCC(Lattice):
    def __init__(self,l_const):
        super(FCC,self).__init__(
            l_const * n.array([0.5,0.5,0]),
            l_const * n.array([0.5,0,0.5]),
            l_const * n.array([0,0.5,0.5])
        )


class BCC(Lattice):
    def __init__(self,l_const):
        super(BCC,self).__init__(
            l_const * n.array([0.5,0.5,-0.5]),
            l_const * n.array([0.5,-0.5,0.5]),
            l_const * n.array([-0.5,0.5,0.5])
        )


class Cubic(Lattice):
    def __init__(self,l_const):
        super(Cubic,self).__init__(
            l_const * n.array([1,0,0]),
            l_const * n.array([0,1,0]),
            l_const * n.array([0,0,1])
        )
    

class Hexagonal(Lattice):
    def __init__(self, a, c):
        super(Hexagonal,self).__init__(
            a*n.array([1,0,0]),
            a*n.array([-0.5,n.sqrt(3)/2,0]),
            c*n.array([0,0,1])
        )
