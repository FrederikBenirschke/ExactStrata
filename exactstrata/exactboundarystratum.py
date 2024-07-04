import sage

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

from sage.modules.free_module import FreeModule  # pylint: disable=import-error
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
# from exactstrata.profile import *
from exactstrata.exactboundarystratum import *
from exactstrata.iteratedblowup import *
from exactstrata.divtautpolynomial import *
from sage.calculus.var import var

import exactstrata.profile as expf 



class ExactBoundaryStratum(SageObject):
    ''' An ExactBoundary stratum inside a GeneralisedStratum is determined by a boundary stratum (encoded by a profile) and a collection of levels,
    which are the levels where the differential is exact.'''
    

    def __init__(self, profile, levels):
        """
        Initializes an ExactBoundaryStratum object with a given profile and levels.

        Parameters:
            profile (expf.Profile): The profile of the boundary stratum.
            levels (list): A list of levels where the differential is exact.

        Returns:
            None

        Raises:
            AssertionError: If the profile is not of the correct type or if the levels are out of range.
        """
       
        self.profile = profile
        self.levels = tuple(sorted(levels, reverse = True))
        assert isinstance(self.profile,expf.Profile), "Profile is not of the right type"
        assert all(level in self.profile.level_set for level in self.levels), "Levels " + str(self.levels) + " are out of range"
        self.exact_stratum = profile.exact_stratum
        assert self.profile in self.exact_stratum.profiles, str(self.profile) + "is not a valid profile"
        self.X = self.exact_stratum.X
        self.bundle = self.exact_stratum.bundle
        
        
        assert self.profile in self.exact_stratum.profiles
        
        
        
        
        # Only need to compute total Chern classes up to the codimension of the stratum of exact differentials
        self.max_codim = self.exact_stratum.exact_rank
    
        self.ELG = self.X.lookup_graph(self.profile)
        self.ranks = self.exact_stratum.all_ranks(self.ELG)

        
    
    def __repr__(self):
        """
        Return a string representation of the ExactBoundaryStratum object.
        """
        return "Exact stratum: " + self.exact_stratum.__repr__()        + ", Boundary stratum: " + self.profile.__repr__() +        ", Levels: " + str(self.levels)
    
    def __eq__(self, other):
        """
        Check if two ExactBoundaryStratum objects are equal.

        Args:
            other (ExactBoundaryStratum): The other ExactBoundaryStratum object to compare with.

        Returns:
            bool: True if the two ExactBoundaryStratum objects have the same profile and levels, False otherwise.
        """
        return self.profile == other.profile and self.levels == other.levels
    
    def __hash__(self):
        """
        Returns a hash value for the object.

        This method returns a hash value for the object based on the combination of the profile and levels.

        Returns:
            int: The hash value of the object.

        """
        return hash(self.profile.tuple() + self.levels)
    
    
    
    
    def __add__(self, other):
        """
        Adds two ExactBoundaryStratum objects together.

        This method takes another ExactBoundaryStratum object as input and returns a new ExactBoundaryStratum object
        that represents the union of the two original objects. The new profile is the concatenation of the profiles of
        the two original objects, and the new levels are the levels that are present in either of the original objects.

        Parameters:
            other (ExactBoundaryStratum): The other ExactBoundaryStratum object to add.

        Returns:
            ExactBoundaryStratum: A new ExactBoundaryStratum object that represents the union of the two original
            objects.
        """
        
        new_profile = self.profile + other.profile
        new_levels = []
        self_map = new_profile.level_map(self.profile)
        other_map = new_profile.level_map(other.profile)
        for level in new_profile.level_set:
            if self_map[level] in self.levels or other_map[level] in other.levels:
                new_levels.append(level)
        return ExactBoundaryStratum(new_profile, new_levels)
                
        
    def are_disjoint(self, other):
        """
        Checks if the intersection of two ExactBoundaryStratum objects is empty.

        Parameters:
            other (ExactBoundaryStratum): The other ExactBoundaryStratum object to check for intersection.

        Returns:
            bool: True if the intersection is empty, False otherwise. If an exception occurs during the addition of the two
            ExactBoundaryStratum objects, True is returned.
        """
       
        try:
            new = self + other
            if  self.exact_stratum.dim - new.codim < _sage_const_0 :
                return True
            return False
        except:
            return True
#         return self.profile.are_disjoint(other.profile)
    
    def is_contained(self, other):
        """
        Checks if the ExactBoundaryStratum self is contained in other.

        Parameters:
            other (ExactBoundaryStratum): The other ExactBoundaryStratum object to check containment against.

        Returns:
            bool: True if self is contained in other, False otherwise.

        This function checks if the profile of self is a degeneration of the profile of other. If this condition is not met,
        False is returned. Otherwise, it adjusts the levels of self based on the level map between the profiles of self and
        other. It then checks if each adjusted level of self is present in the levels of other and not in the levels of self.
        If any of these conditions are not met, False is returned. Otherwise, True is returned.
        """
       

        # The profile of self needs to be a degeneration of the profile of other
        if not (other.profile > self.profile or self.profile == other.profile):
            return False

        adjusted_levels = [self.profile.level_adjustment(other.profile, level) for level in self.levels]
        # self needs to contain all (adjusted) levels of other
        for level in self.profile.level_set:
            if self.profile.level_adjustment(other.profile, level) in other.levels and level not in self.levels:

                return False
        return True



    
    @property
    def codim(self):
        """
        Returns the codimension of A_P^{I} inside a boundary stratum D_Q. The default is P = Q.

        :return: The sum of the ranks of the levels in self.levels.
        :rtype: int
        """
       
        return sum(self.ranks[-level] for level in self.levels)
    
    @property
    def total_codim(self):
        """
        Returns the total codimension of the exact boundary stratum.

        The total codimension is calculated by adding the codimension of the exact boundary stratum
        (which is the sum of the ranks of the levels in `self.levels`) to the length of the profile.

        :return: The total codimension of the exact boundary stratum.
        :rtype: int
        """
        return self.codim + len(self.profile)
    
    
    
    
    
  
    @cached_method
    def formal_chern(self):
        """
        Returns the formal Chern class of the exact boundary stratum.

        The formal Chern class is computed by multiplying the formal Chern classes of the strata at each level.
        More precisely, it copmutes the total Chern class for A_P^[I] inside D_Q, where Q is an undegeneration of P.
        The default is Q=P.
        Warning: This works with reduced substacks and does not take care of multiplicities.

        :return: The formal Chern class of the exact boundary stratum.
        :rtype: symbolic expression
        """
        
#         if ambient_profile == None:
#             ambient_profile = self.profile
#         assert ambient_profile > self.profile or ambient_profile == self.profile
        return prod(self.profile.formal_chern(level) for level in self.levels) 
#           *\ self.profile.formal_total_chern(ambient_profile = ambient_profile)

    @cached_method
    def symbolic_chern(self):
        """
        Returns the symbolic Chern class of the exact boundary stratum.

        The symbolic Chern class is computed by multiplying the symbolic Chern classes of the strata at each level.
        More precisely, it computes the total Chern class for A_P^[I] inside D_Q, where Q is an undegeneration of P.
        The default is Q=P.
        Warning: This works with reduced substacks and does not take care of multiplicities.

        :return: The symbolic Chern class of the exact boundary stratum.
        :rtype: symbolic expression
        """
        return prod(self.profile.symbolic_chern(level) for level in self.levels) 
        
        
    
    
   
    @cached_method
    def product_decomposition(self, other):
        """
        Computes the product decomposition of the class of a stratum of exact differentials.

        Args:
            other (ExactBoundaryStratum): The other exact boundary stratum to compute the product decomposition with.

        Returns:
            tuple: A tuple containing the intermediate decomposition, the remainder profile, and the remainder levels.
        """
       
        assert isinstance(other, ExactBoundaryStratum) and self.is_contained(other)
          
     
        
        intermediate_boundary = ExactBoundaryStratum(self.profile + other.profile, [])
        
        intermediate  = intermediate_boundary + other
        
        
        # remainder_exact has to have the same profile as self
        assert intermediate.profile == self.profile
        
        remainder_profile = self.profile - other.profile
        remainder_levels = tuple(set(self.levels) - set(intermediate.levels))
        return intermediate, remainder_profile, remainder_levels
        
    
    
    
    def formal_normal_bundle(self, other, mult = False):
        """
        Computes the formal normal bundle of the ExactBoundaryStratum self with respect to another ExactBoundaryStratum other.

        More precisely, it computes the total normal bundle of A_P^{[I]} within A_Q^{[J]}.
        We use the normalization where the last coefficient encodes the fundamental class of the reduced space A_P^[I]
        To obtain the class of a nonreduced space one needs to multiply with the correct multiplicity
        The multiplicity will be determined by the IteratedBlowup and all multiplicities will be computed there.

        Args:
            self (ExactBoundaryStratum): The first ExactBoundaryStratum.
            other (ExactBoundaryStratum): The second ExactBoundaryStratum.
            mult (bool, optional): A flag indicating whether to multiply the normal bundles. Defaults to False.

        Returns:
            tuple: A tuple containing the intermediate decomposition, the remainder profile, and the remainder levels.
        """
       
        intermediate, remainder_profile, remainder_levels = self.product_decomposition(other)
        normal_bundle_A = prod(intermediate.profile.formal_chern(level) for level in remainder_levels)
        normal_bundle_B = remainder_profile.formal_normal_bundle(mult = mult)
        
        return normal_bundle_A * normal_bundle_B 
    

        
    
    
