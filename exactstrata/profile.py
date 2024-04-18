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
        return list(self.profile)
    
    def tuple(self):
        return self.profile
    
   
    
    def __eq__(self, other):
        if other == None:
            return False
        if not isinstance(other, Profile):
            return False
        return self.profile == other.profile
    
    def __hash__(self):
        return hash(self.profile)
        
       
        
    def __repr__(self):
        return str(self.profile)
    
    def __getitem__(self, key):
        return self.profile[key]
    

        
       
    def __add__(self, other):
        ''' Returns the sum of two profiles, corresponding to the profile of the intersection of both boundary strata.'''
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
        return self + other
    
    def __sub__(self, other):
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
        return item in self.profile
    
    
    def __len__(self):
        return len(self.profile)
    
    def __iter__(self):
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
        return self.L <=  2  or (self.L == 3  and level == - 1 )
    
    def are_disjoint(self, other):
        try:
            self + other
        except:
            return True
        return False
    
    def is_contained(self, other):

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
        if ambient_profile == None:
            ambient_profile = self.profile
        return self.L - ambient_profile.L
    
    
    
    @cached_method
    def formal_normal_bundle(self, mult = False):
        
        
        
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
        return self.exact_stratum.to_ELG(self.symbolic_chern(level)[-1])
        

    
    
    def formal_ck(self, level, degree):
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
        return [profile for profile in self.exact_stratum.profiles if self.profile_splits_level(profile, level) ]

    @cached_method
    def splitting_bics(self, level):
        bics = self.exact_stratum.bics
        return [bic for bic in  bics if self.bic_splits_level(bic, level)]
    
    
    @cached_method
    def building_set(self, level):
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

