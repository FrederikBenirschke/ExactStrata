import sage

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

from sage.modules.free_module import FreeModule  # pylint: disable=import-error
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.misc.misc_c import prod
#from exactstrata.exactboundarystratum import *
from exactstrata.iteratedblowup import *
from exactstrata.divtautpolynomial import *
from sage.calculus.var import var

from sage.sets.set import Set



class Profile(SageObject):
    ''' A profile is determined by a GeneralisedStratum and a collection of BICs. It encodes the different boundary strata of a stratum of multi-scale differentials.'''
    def __init__(self, exact_stratum, profile):
        """
        Initializes the Profile object with the given exact_stratum and profile.

        Parameters:
            exact_stratum: The GeneralisedStratum object.
            profile: The collection of BICs representing the profile.

        Returns:
            None
        """
        self.exact_stratum = exact_stratum
        self.X = self.exact_stratum.X
        if type(profile) == int or type(profile) == sage.rings.integer.Integer:
            profile = (profile,)
        self.profile = tuple(profile)
        assert self.profile in self.exact_stratum.profile_list, str(self.profile) + "is not a valid profile"
        
        # Number of levels of profile
        self.L = len(self.profile) +  1 
        self.mult = prod(self.X.bics[bic_index].ell for bic_index in self.profile)
        self.codim = self.L -  1 
        
       
        self.ELG = self.X.lookup_graph(self.profile)
        
        self.ranks = self.exact_stratum.all_ranks(self.ELG)
        self.level_set = Set(range(-self.L+ 1 ,  1 ))
        
        
        
    
        
    
  
    
    def list(self):
        """
        A method that returns a list of the profile attribute.
        """
        return list(self.profile)
    
    def tuple(self):
        """
        Returns the profile attribute as a tuple.

        :return: A tuple representing the profile.
        :rtype: tuple
        """
        return self.profile
    
   
    
    def __eq__(self, other):
        """
        A description of the entire function, its parameters, and its return types.
        """
        if other == None:
            return False
        if not isinstance(other, Profile):
            return False
        return self.profile == other.profile
    
    def __hash__(self):
        """
        A description of the entire function, its parameters, and its return types.
        """
        return hash(self.profile)
        
       
        
    def __repr__(self):
        """
        Returns a string representation of the object.

        :return: A string representation of the object.
        :rtype: str
        """
        return str(self.profile)
    
    def __getitem__(self, key):
        """
        Get the value at the specified key from the profile.

        Parameters:
            key (Any): The key to retrieve the value from the profile.

        Returns:
            Any: The value associated with the specified key.
        """
        return self.profile[key]
    

        
       
    def __add__(self, other):
        """
        Returns the sum of two profiles, corresponding to the profile of the intersection of both boundary strata.

        Parameters:
            other (Profile): The other profile to be added. If not an instance of Profile, it will be converted to a Profile using the same exact_stratum.

        Returns:
            Profile: The profile resulting from the intersection of the two boundary strata.

        Raises:
            AssertionError: If the resulting profile is not in the list of profiles of the exact_stratum.
        """
       
        if not isinstance(other, Profile):
            other = Profile(self.exact_stratum, other)
        intersection = self.list()
        for bic_index in other:
            
            # Only need to add new BICs
            if bic_index not in self:
                
                # Find the first position such that index lies over a BIC in some_profile
                for self_index in self:
                    if self.exact_stratum.X.lies_over(bic_index, self_index):
                        i = intersection.index(self_index)
                        intersection.insert(i, bic_index)
                        break
                # Index has not been inserted yet? Then add in the last position       
                if bic_index not in intersection:
                    intersection.append(bic_index)
                    
        assert Profile(self.exact_stratum, intersection) in self.exact_stratum.profiles, "Profiles do not intersect" 
        return Profile(self.exact_stratum, intersection)
    
    def __radd__(self,other):
        """
        Adds the current object to another object on the right side.

        Parameters:
            other (Any): The object to add to the current object.

        Returns:
            Any: The sum of the current object and the other object.
        """
        return self + other
    
    def __sub__(self, other):
        """
        Subtracts the `other` profile from the current profile.

        Parameters:
            other (Profile): The profile to subtract from the current profile.

        Returns:
            Profile: The resulting profile after the subtraction.

        Raises:
            ValueError: If the profiles are not degenerations of each other.
        """
        if self == other:
            return Profile(self.exact_stratum, ())
        if not isinstance(other, Profile):
            other = Profile(self.exact_stratum, other)
        else:
            if other > self:
                return Profile(self.exact_stratum, [index for index in self.profile if index not in other.profile])
            if self > other:
                return Profile(self.exact_stratum,  [index for index in other.profile if index not in self.profile])
        raise ValueError("Profiles" +str(self)+"," +str(other)+ "are not degenerations of each other")
        
    def __rsub__(self, other):
        return self - other
                
    def __contains__(self, item):
        """
        Check if the given BIC is present in the `profile` attribute.

        Parameters:
            item (int): The BIC index to check.

        Returns:
            bool: True if the BIC index is present in the `profile` attribute, False otherwise.
        """
        return item in self.profile
    
    
    def __len__(self):
        """
        Returns the length of the `profile` attribute.

        :return: An integer representing the length of the `profile` attribute.
        :rtype: int
        """
        return len(self.profile)
    
    def __iter__(self):
        """
        Returns an iterator object that can be used to iterate over the elements of the `profile` attribute.

        Returns:
            An iterator object.
        """
        return iter(self.profile)
    
    
   
    def __gt__(self, other):
        '''  Returns true if other is a (proper) degeneration of self. '''
        assert isinstance(other, Profile), str(other) + ' is not a profile.'
        for bic_index in self:
            if bic_index not in other:
                return False
        return True
    
    def is_empty(self):
        ''' Returns true if the profile corresponds to the ambient stratum.'''
        return len(self) ==  0 
    
    def is_reduced(self, level):
        """
        Determines if the profile is reduced.

        Args:
            level (int): The level to check.

        Returns:
            bool: True if the profile is reduced, False otherwise.
        """
        return self.L <=  2  or (self.L == 3  and level == - 1 )
    
    def are_disjoint(self, other):
        """
        Determines if two profiles are disjoint.

        Args:
            self (Profile): The first profile.
            other (Profile): The second profile.

        Returns:
            bool: True if the profiles are disjoint, False otherwise.
        """
        try:
            self + other
        except:
            return True
        return False
    
    def is_contained(self, other):
        """
        Determines if the profile of self is contained within the profile of other.

        Args:
            other (Profile): The other profile to compare against.

        Returns:
            bool: True if the profile of self is contained within the profile of other, False otherwise.
        """

        # The profile of self needs to be a degeneration of the profile of other
        if not (other.profile > self.profile or self.profile == other.profile):
            return False

        adjusted_levels = [self.level_adjustment(other, level) for level in self.levels]
        # self needs to contain all (adjusted) levels of other
        for level in other.levels:
            if level not in adjusted_levels:
                return False
        return True
    
   
    
    @cached_method
    def codim(self, ambient_profile = None):
        """
        Calculate the codimension of the object with respect to the given ambient profile.
        
        Parameters:
            ambient_profile (Profile, optional): The ambient profile to calculate the codimension with. Defaults to None.
        
        Returns:
            int: The codimension of the object.
        """
        if ambient_profile == None:
            ambient_profile = self.profile
        return self.L - ambient_profile.L
    
    
    
    @cached_method
    def formal_normal_bundle(self, mult = False):
        """
        Computes and returns the formal normal bundle of the exact stratum profile.

        Args:
            mult (bool, optional): A flag indicating whether to multiply the normal bundle by the ell value. Defaults to False.

        Returns:
            DivTautPolynomial: The computed formal normal bundle.
        """
        
        
        
        chern = DivTautPolynomial([1])
        
        if not mult:
            for index in self:
                    chern *= DivTautPolynomial([1, self.exact_stratum.bic_variables[index]])
        else:
            for index in self:
                    chern *= DivTautPolynomial([1, self.X.bics[index].ell*self.exact_stratum.bic_variables[index]])
            
            
        
        return chern
    
 
    

    @cached_method
    def formal_chern(self, level):
        ''' Formal total Chern class of the stratum of exact differentials at level i  
        inside the boundary stratum D_P.'''
        
        
        
        codim = self.ranks[-level]
        if codim ==  1 :
            blowup = self.exact_stratum.blowup(self, level)
            chern, _ = blowup.zero_locus_chern()
            return chern
            
        else:
            return DivTautPolynomial([self.formal_ck(level,k) for k in range( 0 ,codim+ 1 )])
        


    def replace_formal_by_symbolic(self, divTaut):
        """
        Replaces the formal Chern classes in the given `divTaut` polynomial with their symbolic equivalents.

        Parameters:
            divTaut (DivTautPolynomial): The polynomial to replace the formal Chern classes in.

        Returns:
            DivTautPolynomial: The polynomial with the formal Chern classes replaced by their symbolic equivalents.
        """
        for level in self.level_set:
            divTaut = divTaut.substitute( self.formal_chern(level)[ 1 :],                                         self.symbolic_chern(level)[ 1 :])
        return divTaut
    
    
    @cached_method
    def symbolic_chern(self, level, avoiding_blowup = True):
        ''' Recursively computes the total Chern class as a symbolic expression.'''
        blowup = self.exact_stratum.blowup(self, level)
        
        
        # If we are only interested in the fundamental class
        # and every blowup center has the expected codim,
        # then we dont have to blowup
        if self.X.g_H[ 0 ] < 2  and self.profile == () and level == 0  and blowup.has_expected_codim() and avoiding_blowup:
            print('Avoiding blowup')
            # We only compute the class but not the normal bundle
            return DivTautPolynomial([ 0 ]), blowup.class_no_blowup()
        
        chern, _ = blowup.zero_locus_chern()
        if blowup.num_of_blowups ==  0 :
            return chern
        else:
            blowdown = blowup.blowdown(chern)
            for profile in blowup.blownup_profiles:
                
                blowdown = profile.replace_formal_by_symbolic(blowdown)
        return blowdown
    
    
    @cached_method
    def ELG_chern(self, level):
<<<<<<<<<<<<<  ✨ Codeium AI Suggestion  >>>>>>>>>>>>>>
+        """
+        Returns the Chern class of the Exact Stratum at the given level in the ELG of the moduli space of multi-scale differentials.
+
+        Parameters:
+            level (int): The level at which to compute the Chern class.
+
+        Returns:
+            ELGTautClass: The Chern class of the Exact Stratum at the given level in the ELG.
+        """
<<<<<  bot-526980c3-5d0c-44da-a59b-45636325765e  >>>>>
        return self.exact_stratum.to_ELG(self.symbolic_chern(level)[-1])
        

    
    
    def formal_ck(self, level, degree):
        """
        A function that returns the formal Chern class at a given level and degree.

        Parameters:
            level (int): The level at which to compute the Chern class.
            degree (int): The degree of the Chern class.

        Returns:
            The formal Chern class as a symbolic variable.
        """
        if degree <  0  or degree > self.exact_stratum.exact_rank:
            return  0 
        if degree ==  0 :
            return  1 
        
        profile_str = 'c' + ''.join(['_{}'.format(i) for i  in self])
        
        

#         latex_str = 'c_' + ''.join(['{}'.format(i) for i  in self])\
#                           + ''.join('_{}_{}'.format(-level,degree))
        
        
        return var(profile_str + ''.join('__{}_{}'.format(-level,degree)))
        
 
    
    def level_adjustment(self, undegeneration, level):
        ''' Determines the level map to an undegeneration of profile.'''
        # No adjustment needed if both profiles are the same
        if self == undegeneration:
            return level
        
        #Divisorial degeneration?
        if len(undegeneration) + 1  == len(self):
            bic = (undegeneration-self)[ 0 ]
            index = self.profile.index(bic)
            if level >= -index :
                return level
            else:
                return level +  1 

        else:
            
            # If the degeneration is not divisorial we write it as a chain of divisorial ones
            last_bic = (self-undegeneration)[- 1 ]
            intermediate = Profile(self.exact_stratum, [index for index in self.list() if index !=last_bic])
            
            new_level = self.level_adjustment(intermediate, level)
            return intermediate.level_adjustment(undegeneration, new_level)
        
    def level_map(self, other):
        # Need other to be an undegeneration of self
        assert other > self, "Argument is not an undegeneration"
        
        level_map = { i:self.level_adjustment(other, i)                                for i in self.level_set}
        return level_map
        
   
    
    def profile_splits_level(self, degeneration, i):
        ''' Returns the profile that (proper) degeneration of self obtained by splitting level i.'''
        
        if degeneration == self:
            return True
        if not  self > degeneration:
            return False
        remainder = degeneration - self

        for index in remainder:
            if not self.bic_splits_level(index, i):
                return False
        return True



    def bic_splits_level(self, bic_index, i):
        X = self.exact_stratum.X
        if i >  0  or i < -len(self):
            
            raise ValueError("Level not in range")
        if len(self) ==  0 :
            return True
        index = -i
        if  index- 1  >=  0 :
            if not X.lies_over(self[index- 1 ], bic_index):
                return False

        if index <= len(self)- 1 :
            if not X.lies_over(bic_index, self[index]):
                return False

        return True
    
    
    def splitting_profiles(self, level):
        """
        Returns a list of profiles that can be obtained by splitting the current profile at the given level.

        Parameters:
            level (int): The level at which the profile is to be split.

        Returns:
            list: A list of profiles that can be obtained by splitting the current profile at the given level.
        """
        return [profile for profile in self.exact_stratum.profiles if self.profile_splits_level(profile, level) ]

    @cached_method
    def splitting_bics(self, level):
        """
        Returns a list of boundary intersection cycles (BICs) from the `self.exact_stratum.bics` list
        that can be split at the given `level`.

        Parameters:
            level (int): The level at which the BICs are to be split.

        Returns:
            list: A list of BICs that can be split at the given `level`.
        """
        bics = self.exact_stratum.bics
        return [bic for bic in  bics if self.bic_splits_level(bic, level)]
    
    
    @cached_method
    def building_set(self, level):
        """
        This function calculates the building set and building divisors based on the given level.
        
        Parameters:
            level (int): The level at which the building set is to be calculated.
        
        Returns:
            tuple: A tuple containing two dictionaries - building_set and building_divs.
        """
        building_set = {}
        building_divs = {}
        for bic_index in self.splitting_bics(level):
            bd = self.exact_stratum.basic_exact_stratum(self+bic_index, level)
            
            #Ignore strata of negative codimension, since they have necessarily zero class
            if self.exact_stratum.dim - bd.total_codim >= 0 :
            
                # A_{P+bic}^{[i]} is a divisor in D_P exactly if A_{P+bic}^{[i]}= D_{P+bic}
                if bd.codim ==  0 : # Reminder: codim of exact boundary stratum A_P^{i} is within D_P     
                    building_divs[bic_index] = bd
                else:    
                    # Exact Boundary Strata of codim inside D_P > rank of the globally defined bundle have to be empty   
                    if bd.codim + 1  <= self.ranks[-level]:
                        building_set[bic_index] = bd
        return building_set, building_divs

