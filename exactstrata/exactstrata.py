
from sage.all_cmdline import *   # import sage library

#  1 = Integer(1);  _2 = Integer(2);  3 = Integer(3)#######################################################
# A (generalised) stratum of exact differentials.
# An ExactStratum is uniquely identified by the following information:
#
#   * X: a GeneralisedStratum, corresponding to the ambient space 
#   * bundle : list of relative homology cycles, i.e. [R_1,...,R_n] where each R_l is
#      a list of integers, corresponding to marked zeros of the differentials. The entries in R_l need 
#      to add to zero and determine which relative periods are zero
# 
#  The main functionalities are
# exact_stratum_class(): Returns the fundamental class of the locus of lambda exact differentials 
#               as a SymbolicExpression representing a divisorial tautological class
# 
# to_ELG(symbolic_class): Takes  a SymbolicExpression representing a divisorial tautological class
# and returns the corresponding ELGTautclass
#
# pushforward(ELGTautClass): Takes an ELGTautClass and computes the pushforward to the moduli space Mg,n

class ExactStratum(SageObject):
    
    
    
    def __init__(self, X, bundle):
        self.X = X
        self._dim = None
        self.bundle = bundle
        self.bundle = self.bundle_basis()
        self._total_rank =  None
        
        # The rank of the residue local system
        self.residue_rank = sum(max(0 , len(signature.pole_ind)- 1 ) for signature in self.X.sig_list())
        # The rank of the relative period local system
        self.torsion_rank = len(self.bundle)
        self._profiles = None
        
        self._bics = None
        
        
        

        self._boundary_strata = None
       
    # Shortcut for the profile () corresponding to the ambient GeneralisedStratum   
    @property    
    def empty_profile(self):
        return Profile(self,())
 
        
   
    # The dimension of the ambient GeneralisedStratum
    @property
    def dim(self):
        if self._dim == None:
           
            self._dim = self.X.dim()
        return self._dim
    
    
    
    # The rank of the globally defined vector bundle or equivalently the codimension of the stratum
    # of exact differentials within the GeneralisedStratum
    @property
    def total_rank(self):
        if self._total_rank == None:
            
            # The rank is the dimension of absolute cohomology of the punctured(!) surface
            # Plus the number of linearly independent relative equations
            
            self._total_rank = self.exact_rank +self.residue_rank + self.torsion_rank
        return self._total_rank
    
    
    # The rank of the absolute cohomology (not inlcuding residues)
    @property
    def exact_rank(self):
        return  2 *sum(self.X.g) 
    
    
    # Lists all BICs such that the restriction of the residue subspace containing i many marked poles
    # to the two level graph is the same as the restriction containing only i-1 marked poles    
    def res_bics(self, index):
        assert index >= 1  and index<= self.residue_rank
        return [bic for bic in self.bics_unordered() if self.residue_ranks_bic(bic)[index]== self.residue_ranks_bic(bic)[index- 1 ]]
    
    
    # Same as res_bics but instead for the bundle of relative periods
    def torsion_bics(self, index):
        assert index >= 1  and index<= self.torsion_rank
        return [bic for bic in self.bics_unordered() if self.torsion_ranks_bic(bic)[index]== self.torsion_ranks_bic(bic)[index- 1 ]]
        
    # The class of the locus of residueless differentials
    # as SymbolicExpression
    def res_class(self):
        result =  1 
        for index in range(1,self.residue_rank+1):
            result *= var('xi')+sum(self.X.bics[bic].ell * self.bic_variables[bic] for bic in self.res_bics(index))
        
        if self.to_ELG(result) == self.X.ZERO:
            return 0 
        return (- 1 )**(self.residue_rank) * result
    
    # The class of the locus of  lambda-exact differentials inside the stratum of exact differentials
    # as SymbolicExpression
    def torsion_class(self):
        result =  1 
        for index in range(1,self.torsion_rank+1):
            result *= var('xi')+sum(self.X.bics[bic].ell * self.bic_variables[bic] for bic in self.torsion_bics(index))
        return (- 1 )**(self.torsion_rank) * result
    
    
    # Product of class of the locus of residueless differentials and the class 
    # of the locus of  lambda-exact differentials inside the stratum of exact differentials
    # as SymbolicExpression
    # In general this class does not have geometric meaning
    # in g=0 this is the class of the stratum of lambda-exact differentials
    def res_torsion_class(self):
        return self.res_class() * self.torsion_class()
    
    
    # Returns the ranks of a total flag for the local system of relative homology to a two level graphs 
    def torsion_ranks_bic(self, bic_index):
        torsion_ranks = []
        
        for index in range(len(self.bundle)+ 1 ):
            sliced_bundle = self.bundle[:index]
            
            exc = ExactStratum(self.X, sliced_bundle)
            
            bund = len(exc.bundle_to_top(self.X.bics[bic_index]).bundle)
            
            torsion_ranks.append(bund)
        return torsion_ranks
    
    # Returns the ranks of a total flag for the local system of residues to a two level graphs
    def residue_ranks_bic(self, bic_index):
        
        X = self.X
        all_poles = X._polelist
#         show('All poles: ',all_poles)
        BIC = self.X.bics[bic_index]
        bic_poles = [BIC.dmp_inv[pole] for pole in all_poles]
#         show('BIC poles:', bic_poles)
        
        LG = BIC.LG
        res_list = []
        for index in range(len(all_poles)+ 1 ):
            index_poles = bic_poles[:index]
#             show('Index: ', index, 'Index poles:', index_poles)
#             bic_poles = [leg for leg in bic_poles if BIC.LG.levelofleg(leg)==0 ]
                
            residue_conditions = 0 
            for v in range(len(LG.genera)):
                if LG.levelofvertex(v) == 0 :



                    poles_on_vertex = len([l for l in LG.legsatvertex(v) if LG.orderatleg(l) < 0 ])
                    poles_from_list = len([pole for pole in index_poles if LG.vertex(pole) == v])
#                     show(v, [l for l in LG.legsatvertex(v) if LG.orderatleg(l) < 0], [pole for pole in index_poles if LG.vertex(pole) == v])
                    # Need to take  into account that residues add up to zero
                    if poles_on_vertex == poles_from_list:
                        residue_conditions += max(0 ,poles_from_list -1 )
                    else:
                        residue_conditions += poles_from_list
            res_list.append(residue_conditions)

        return res_list

       
    
    
    
    # Returns a list of integers
    # The i-th entry is the the number of independent residue cycles when using i marked poles
    # Which poles are chosen depends on the ordering of the signature of the underlying GeneralisedStratum
    def _residue_ranks(self,pole_list):
        X = self.X
        all_poles = X._polelist
        

        def residue_rank(index):


            assert index <= len(all_poles)
            pole_list = all_poles[:index]
           # show(pole_list)
            residue_conditions = 0 
            for i, sig in enumerate(X.sig_list()):
                poles_on_vertex = len(sig.pole_ind)
                poles_from_list = len([pole for pole in pole_list if (pole[0 ]==i and pole[ 1 ] in sig.pole_ind)])

                # Need to take  into account that residues add up to zero
                if(poles_on_vertex==poles_from_list):
                    residue_conditions += max(0 ,poles_from_list - 1 )
                else:
                    residue_conditions += poles_from_list

            return residue_conditions

        return [residue_rank(index) for index in range(len(all_poles)+ 1 )]

    
    
                
    # Creates a profile from a tuple
    @cached_method
    def profile(self, profile_tuple):
        return Profile(self, profile_tuple)
    
    # Returns all profiles (as Profiles) of the underlying GeneralisedStratum, ordered by ascending dimension
    @property
    def profiles(self):
        if self._profiles == None:
            self._profiles = sorted([Profile(self, pf) for pf in self.profile_list],key=len, reverse=True)
        return self._profiles
    
    # All list of all profiles (as tuples)
    # Used for checking that a Profile is valid
    @property
    def profile_list(self):
        return [pf for pf_list in self.X.lookup_list for pf in pf_list]
        
    
    
    
    #Returns a list of all indices of BICs, respecting the partial order on two-level graphs
    @property
    def bics(self):
        if self._bics == None:
            old_list = list(range(len(self.X.bics.copy())))
            new_list = []
            
            # We need to reorder the list, to ensure that it preserves the partial order on level graphs
           
            while len(old_list) > 0 :
                index = old_list.pop()
                
                for pos, other_index in enumerate(new_list):
                    
                    
                    # LG.index > LG.other_index ?
                    if self.X.lies_over(index, other_index):
                        new_list.insert(pos, index)
                        break
                if index not in new_list:
                    new_list.append(index)
            self._bics = new_list
        return self._bics


    
    #Returns a list of all indices of BICs, not necessarily respecting the partial order
    @cached_method
    def bics_unordered(self):
        return list(range(0,len(self.X.bics)))

    
    
#     @cached_method
#     def mult(self, profile, amb_profile):
#         return prod(self.X.bics[bic_index].ell for bic_index in remainder(amb_profile, profile))
    
    
    
    
    # Returns the boundary stratum A_P^{[i]}        
    def basic_exact_stratum(self, profile, i):
#         assert type(i) is sage.rings.integer.Integer, str(i) + " is not an integer but has type "+ str(type(i))
        return self.exact_boundary_stratum(profile, (i,))
    
    # Creates an ExactBoundaryStratum from a profile (Profile or tuple)
    # and a level (int), corresponding to subscheme A_P^{[i]} of exact differentials 
    # exact on the top level
    @cached_method
    def exact_boundary_stratum(self, profile, levels):
        profile = Profile(self, profile)
        return ExactBoundaryStratum(profile, levels)
    
#     # Creates an ExactBoundaryStratum without levels from a profile (Profile or tuple)
#     # Corresponds to a boundary stratum D_P 
#     @cached_method
#     def boundary_stratum(self, profile):
#         return self.exact_boundary_stratum(profile, ())
    
    
    
    
#     @cached_method
#     def xi_at_level(self, profile, i):
#         assert i in profile.level_set
#         if len(profile) == 0 or (len(profile)==1 and i==0):
#             return self.X.xi
#         else:
#             #Remaining case: two levels and i=-1, here we need to use Prop. 3.8
#             if len(profile) == 1:
#                 bic= profile[0]

#                 xi_bot = self.X.xi + profile.mult* self.D_lg(bic)
#                 xi_bot += sum([self.X.bics[nbic].ell * self.D_lg(nbic) for nbic in self.bics if self.X.lies_over(nbic, bic)]) 
#                 return xi_bot
#             if i == 0:
#                 return self.xi_level(Profile(self, (profile[0],)), 0)
#             else:
#                 return self.xi_level(Profile(self, (profile[-i-1],)), -1)
            
    
    # Returns the ELGTautClass of a BIC (as integer).
    @cached_method
    def D_lg(self, bic):
        return self.X.additive_generator(((bic,),0 ))
    
    
    # Returns dictionary convering from a SymbolicExpression 'D_{bic}' to the ELGTautClass
    # of the corresponding boundary divisor
    @property
    def divisor_classes(self):
        return {var('D_{}'.format(bic)) : self.D_lg(bic) for bic in self.bics}
    
    
    # Returns the xi class at a level of a boundary stratum D_P as a SymbolicExpression.
    # The result is a sum of 'xi' and 'D_{bic}' for different BICs
    @cached_method    
    def formal_xi_at_level(self, profile, level):
        assert level in profile.level_set, show(level, profile.level_set)
        if len(profile) == 0  or  level==0 :
            return var('xi')
        if len(profile) ==  1 :
            bic_index = profile[0 ]
            xi_bot = var('xi') + self.X.bics[bic_index].ell* self.bic_variables[bic_index] 
            xi_bot += sum([self.X.bics[nbic].ell * self.bic_variables[nbic]                            for nbic in self.bics if self.X.lies_over(nbic, bic_index)]) 
            return xi_bot
            
           
        else:
            return self.formal_xi_at_level(Profile(self, (profile[-level- 1 ],)), - 1 )
        
       
        

    # Dictionary converting formal xi classes and boundary classes 'D_{bic}' to ELGTautClasses
    @cached_method
    def formal_to_AG_dic(self):
        return self.divisor_classes | {var('xi'):self.X.xi }
                   
    @property
    def bic_variables_inv(self):
        return { var('D_{}'.format(i)): i for i,_ in enumerate(self.bics)} 
    
    @property 
    def bic_variables(self):
        return {value: key for key, value in self.bic_variables_inv.items()}
 
    def __repr__(self):
        return str([[sig.sig for sig in self.X.sig_list()], self.bundle])  
    
    
      ######################################################### 
    # -------- Blow up algorithm --------
    #########################################################
    #### Input: amb_profile - Profile
    #####.      i - level
    #####
    
    
    # Constructs an iterated blowup. 
    # Given a boundary stratum D_P and a level i
    # this is an iterated blow up along A_{P,Q}^{[i]} where P+Q is a degeneration of P
    # obtained by splitting the i-th level (except for boundary strata such that A_{P,Q}^{[i]})
    # is a divisor in D_P)
    # 
    def blowup(self, amb_profile , i):

        amb_profile = self.profile(amb_profile)

        assert i in amb_profile.level_set
        # Step 1) Find the building set, i.e. all two-level graphs that split level i

        building_set, building_divs = amb_profile.building_set(i)

        blowup_centers = list(building_set.values())
        blowup = IteratedBlowup(self, amb_profile, i, self.exact_boundary_stratum(amb_profile,()),                                 blowup_centers)

        return blowup
    
    #--------------------------------------------------------------------------
    # Chern polynomial logic
    #--------------------------------------------------------------------------
    
    # Returns the fundamental class of a stratum of lambda-exact differentials.
    # The result is a SymbolicExpression in terms of 'xi' and 'D_{bic}'.
    # 
    def exact_stratum_class(self, avoiding_blowup=True):
        # We first reduce to the case without relative periods
        if len(self.bundle) > 0 :
            new_exact_stratum = ExactStratum(self.X, [])
            return new_exact_stratum.exact_stratum_class(avoiding_blowup = avoiding_blowup) *  self.torsion_class()
        
        # Now we deal with all residue equations
        else:
            res_cl = self.res_class()
            if res_cl == 0 :
                return 0 
            else:
                # Finally we need to compute the class corresponding to the remaining ansolute periods
                return self.res_class() * self.symbolic_class(self.empty_profile, 0, avoiding_blowup= avoiding_blowup)
        
    
    # Returns the fundamental class of A_P^{[i]} inside D_P.
    # We are really working inside the stratum of residueless differentials.
    # If one wants the class inside the GenerlisedStratum this needs to be multiplied by
    # self.res_class()
    def symbolic_class(self, profile, level, avoiding_blowup = True):
        profile = Profile(self, profile)
        
        # The class of the exact stratum is the last coefficient in the Chern polynomial
        return profile.symbolic_chern(level, avoiding_blowup=avoiding_blowup)[-1]
        
               
    

    # Converts a SymbolicExpression to an ELGTautClass.
    # Requires the symbolic expression to be a polynomial in 'xi' and 'D_{bic}'.
    def to_ELG(self, symbolic_class):
        X = self.X

        if  type(symbolic_class) is int:
            return symbolic_class * X.ONE
        
        result = []
        xi = var("xi")

        symbolic_class = symbolic_class.expand()
        for product in symbolic_class.iterator():
            tmp = X.ONE
            for D in product.variables():
                if D == xi:
                    tmp *= X.xi_pow(product.degree(D))
                else:
                    tmp *= self.formal_to_AG_dic()[D]**product.degree(D)
            result.append(tmp)
        return X.ELGsum(result)
    
    
    # Computes the pushforward of a ELGTautClass on the underlying GeneralisedStratum
    # to a TautClass on Mgnbar.
    def pushforward(self, fund_class):
        return fund_class.to_prodtautclass().pushforward()

    
  
    #--------------------------------------------------------------------------
    # Bundle logic
    #--------------------------------------------------------------------------
    
  
        
    # The bundle records the relative equations and  is given by a list of vectors. 
    # Returns a basis for the space of relative cycles.
    def bundle_basis(self):
        
        num_marked_points =  sum([len(sig.sig) for sig in self.X.sig_list()]) #number of marked points

        point_module = FreeModule(QQ, num_marked_points)
        for vector in self.bundle:
            if(sum(vector)!=0 ):
                raise ValueError("Not a valid relative homology class")
        return list(point_module.submodule(self.bundle).basis())

  
    
    # Restricts the ExactBoundaryStratum to the top level of a BIC.
    # The result is a new ExactBoundaryStratum supported on the underlying GenerlisedStratum of the top level
    # and computes the images of absolute and relative cycles under the specialization morphism.
    def bundle_to_top(self, BIC):
        
        if not BIC.is_bic():
            raise ValueError("Level graph is not a BIC")
        marked_topl ,num_marked_topl, rel_module, topl_dict = self.top_level_homology(BIC)
        total_relations = self.relative_bundle_basis(BIC) + self.absolute_bundle_basis(BIC)
        basis = rel_module.submodule(total_relations).basis()
        
        return ExactStratum(BIC.top,list(basis))
    
        
    # Helper function to organize the data associated to the top level of a bic.
    # Input:
    # 'ELG' - an EmbeddedLevelGraph, should only be used for BICs
    # Output:
    # 'marked_topl' - list of integers, the entries correspond to the marked points  on the top level of the BIC
    # 'num_marked_topl' - integer, the number of marked points on top level
    # 'topl_dict' - dictionary (marked point on top level of BIC: index_marked point ) enumerating the marked points
    def top_level_homology(self, ELG):
        marked_topl = [l for i,_ in enumerate(ELG.LG.genera)                       if ELG.LG.levelofvertex(i)==0  for l in ELG.LG.legsatvertex(i) ]
        num_marked_topl = len(marked_topl)
        rel_module = FreeModule(QQ, num_marked_topl)
        topl_dict = {point : index for index, point in enumerate(marked_topl)}
        return marked_topl, num_marked_topl, rel_module, topl_dict
    
    
    
    # For a BIC finds all relative cycles in the image of the absolute homology under the specialization map 
    # to the top level
    def absolute_bundle_basis(self, ELG):
        
        basis= ELG.LG.stgraph.cycle_basis()
        marked_topl , num_marked_topl, _,topl_dict = self.top_level_homology(ELG)
        
        cycle_basis =[]
        for path in basis:
            vec = [0 ] * num_marked_topl
            for i,e in  enumerate(path):
                leg = ELG.LG.edges[i][0]
                vec[topl_dict[leg]]= e
            cycle_basis.append(vec)
        return cycle_basis
    
        
    # Computes a basis for the image of the relative cycles (here we mean purely relative)
    # under the specialization morphism to the top level.
    def relative_bundle_basis(self, ELG):
        

        #only works for BICs:
        if not ELG.is_bic():
            raise ValueError("ELG ist not a BIC")
            
        marked_topl, num_marked_topl, _, topl_dict = self.top_level_homology( ELG)

        #First step: need to write the element in relative homology as sum of paths
        # an element i nrelative homology has closed boundary so that is always possible
        new_bundle =[]
        dic = ELG.dmp_inv # dictionary sending marked points in ambient generalsed stratum to marked points in ELG
        for vector in self.bundle:
            if sum(vector)!=0 :
                raise ValueError("Not a relative homology class")
            new_vector = [0 ] * num_marked_topl
            for i, x in enumerate(vector):
                #we choose x * path from i-th marked point to some marked point
                #ensure that the marked point is on the same connected component
                start = ELG.dmp_inv[marked_to_component(self.X)[i]]
                #end = ELG.dmp_inv[(0,0)]
                end = ELG.dmp_inv[(marked_to_component(self.X)[i][0 ],0 )]
                path = path_between_markings(ELG, start, end)
                #now go through the path and add x*(j+1-th point - j-th point) but only for the points in pts,
                #since the other points correspond to marked points in lower levels
                for j in range(len(path)- 1 ):

                    a,b = path[j],path[j+ 1 ]
                    if ELG.LG.levelofleg(a)==0  and ELG.LG.levelofleg(b)==0 :
                        new_vector[topl_dict[a]]+= -x  
                        new_vector[topl_dict[b]]+= x
            new_bundle.append(new_vector)
        return new_bundle
    
    
    

    

    # Given an EmbeddedLevelGraph computes the exact rank of the restriction of the bundle of relative homology
    # Note that this does not take into account residues or relative cycles
    # but only absolute homology hmology cycles
    
    def all_ranks(self, ELG):
        

        # If the input ELG is a profile instead of an EmbeddeedLevelGraph,
        # we can use any graph with this profile.
        # The result is independent of the choice
        # and such a graph can be obtained using self.X.lookup_graph(profile)
        if type(ELG) == list or type(ELG) == tuple:
            ELG = self.X.lookup_graph(ELG)
            
        # Needs to be an EmbeddedLevelGraph embedded into the stratum self.X
        if ELG.X!= self.X:
            raise ValueError("ELG is defined for the wrong stratum")
        ranks =[]
        
        # By squishing a level graph we obtain a two level, where the rank of the restriction to the top level
        # is the sum of the ranks of V restricted to all levels that are being squished to the top level
        # This allows to compute the rank of the restriction to every level recursively
        for i in range( 1 ,len(ELG.dlevels)):
            squish = ELG.delta(i)
#             squish_index = self.X.bics.index(squish)
            rank_i = self.bundle_to_top(squish).total_rank - self.bundle_to_top(squish).residue_rank
            ranks.append(rank_i - sum(ranks))
            
        # The rank on the bottom level is the rank of V minus the rank on all higher levels
        ranks.append(self.exact_rank - sum(ranks))
        return ranks
                 

class ExactBoundaryStratum(SageObject):
    

    def __init__(self, profile, levels):
        self.profile = profile
        self.levels = tuple(sorted(levels, reverse = True))
        assert isinstance(self.profile,Profile), "Profile is not of the right type"
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
        return "Exact stratum: " + self.exact_stratum.__repr__()        + ", Boundary stratum: " + self.profile.__repr__() +        ", Levels: " + str(self.levels)
    
    def __eq__(self, other):
        return self.profile == other.profile and self.levels == other.levels
    
    def __hash__(self):
        return hash(self.profile.tuple() + self.levels)
    
    
    
    # Intersects two exact boundary strata
    # The underlying profile is the sume of the profiles
    # The level set is the union of level sets, adjusted by the level maps of the undegeneration.
    def __add__(self, other):
        new_profile = self.profile + other.profile
        new_levels = []
        self_map = new_profile.level_map(self.profile)
        other_map = new_profile.level_map(other.profile)
        for level in new_profile.level_set:
            if self_map[level] in self.levels or other_map[level] in other.levels:
                new_levels.append(level)
        return ExactBoundaryStratum(new_profile, new_levels)
                
        
    def are_disjoint(self, other):
        try:
            new = self + other
            if  self.exact_stratum.dim - new.codim < 0 :
                return True
            return False
        except:
            return True
#         return self.profile.are_disjoint(other.profile)
    
    def is_contained(self, other):

        # The profile of self needs to be a degeneration of the profile of other
        if not (other.profile > self.profile or self.profile == other.profile):
            return False

        adjusted_levels = [self.profile.level_adjustment(other.profile, level) for level in self.levels]
        # self needs to contain all (adjusted) levels of other
        for level in self.profile.level_set:
            if self.profile.level_adjustment(other.profile, level) in other.levels and level not in self.levels:

                return False
        return True



    # Computes the codimension of A_P^{I} inside a  boundary stratum D_Q
    # Default is P = Q
    @property
    def codim(self):
        return sum(self.ranks[-level] for level in self.levels)
    
    @property
    def total_codim(self):
        return self.codim + len(self.profile)
    
    
    
    
    
    # Formal expression for the total Chern class for A_P^[I] inside D_Q, where Q is an undegenration of P
    # By default Q=P
    # Warning: This works with reduced substacks and does not take care of multiplicities
    @cached_method
    def formal_chern(self):
#         if ambient_profile == None:
#             ambient_profile = self.profile
#         assert ambient_profile > self.profile or ambient_profile == self.profile
        return prod(self.profile.formal_chern(level) for level in self.levels) 
#           *\ self.profile.formal_total_chern(ambient_profile = ambient_profile)

    @cached_method
    def symbolic_chern(self):
        return prod(self.profile.symbolic_chern(level) for level in self.levels) 
        
        
    
    
    # Decomposes
    @cached_method
    def product_decomposition(self, other):
        assert isinstance(other, ExactBoundaryStratum) and self.is_contained(other)
          
     
        
        intermediate_boundary = ExactBoundaryStratum(self.profile + other.profile, [])
        
        intermediate  = intermediate_boundary + other
        
        
        # remainder_exact has to have the same profile as self
        assert intermediate.profile == self.profile
        
        remainder_profile = self.profile - other.profile
        remainder_levels = tuple(set(self.levels) - set(intermediate.levels))
        return intermediate, remainder_profile, remainder_levels
        
    
    
    # Computes the total normal bundle of A_P^{[I]} within A_Q^{[J]} 
    # We use the normalization so that the 
    # The last coefficient encodes the fundamental class of the reduced space A_P^[I]
    # To obtain the class of a nonreduced space one needs to multiply with the correct multiplicity
    # The multiplicity will be determined by the IteratedBlowup and all multiplicities will be computed there
    def formal_normal_bundle(self, other, mult = False):
        intermediate, remainder_profile, remainder_levels = self.product_decomposition(other)
        normal_bundle_A = prod(intermediate.profile.formal_chern(level) for level in remainder_levels)
        normal_bundle_B = remainder_profile.formal_normal_bundle(mult = mult)
        
        return normal_bundle_A * normal_bundle_B 
    

        
    
    

#Helper functions


# Dictionary listing the marked points of a generalised stratum     
def marked_to_component(X):
    return { i+j : (i,j) for i,sig in enumerate(X.sig_list()) for j,_ in enumerate(sig.sig)}


            

    
    
def path_between_markings(ELG, start, end):
    G=ELG.LG
    return bfs(G,start,end)

# Breadth-first search to find a path between two vertices
# Only works if both vertices are in the same connected component of the graph  
def bfs(LG, start_leg,goal_leg):
    visited = [] # List of (vertex,leg) to keep track of visited nodes.
    path_dict= {}
    queue = []     # Initialize a queue
    visited.append([LG.vertex(start_leg),start_leg,0 ])

    queue.append(LG.vertex(start_leg))

    while queue:
        s = queue.pop(0 ) 

    #print (s, end = " ") 
        neighbors = []
        for i,j in  LG.edgesatvertex(s):
            if LG.vertex(i)== s:
                neighbors.append([LG.vertex(j),j,i])
            if LG.vertex(j)==s:
                neighbors.append([LG.vertex(i),i,j])

        for vertex,leg_i,leg_o in neighbors:
            if vertex not in [visit_vertex for visit_vertex,_,_ in visited]:
                visited.append([vertex,leg_i,leg_o])
                queue.append(vertex)
            if vertex==LG.vertex(goal_leg):

                # Now need to compute the path from goal_leg to start_leg
                position=vertex

                path=[goal_leg]
                while(position!=LG.vertex(start_leg)):


                    for V in visited:
                        if V[0 ]==position:
                            path+=[V[ 1 ],V[ 2 ]]
                            break
                    position = LG.vertex(path[- 1 ])

                path.append(start_leg)
                path.reverse()
                return path


class Profile(SageObject):
    

    

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
    

        
     #Sum of profiles   
    def __add__(self, other):
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
    
    
    # Returns true if other is a (proper) degeneration of self
    def __gt__(self, other):
        assert isinstance(other, Profile), str(other) + ' is not a profile.'
        for bic_index in self:
            if bic_index not in other:
                return False
        return True
    
    def is_empty(self):
        return len(self) == 0 
    
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
        
        
        
        chern = DivTautPolynomial([ 1 ])
        
        if not mult:
            for index in self:
                    chern *= DivTautPolynomial([ 1 , self.exact_stratum.bic_variables[index]])
        else:
            for index in self:
                    chern *= DivTautPolynomial([ 1 , self.X.bics[index].ell*self.exact_stratum.bic_variables[index]])
            
            
        
        return chern
    
 
    
#     # Formal total Chern class of the stratum of exact differentials at level i  
#     # inside the boundary stratum D_P
    @cached_method
    def formal_chern(self, level):
        
        
        
        codim = self.ranks[-level]
        if codim ==  1 :
            blowup = self.exact_stratum.blowup(self, level)
            chern, _ = blowup.zero_locus_chern()
            return chern
            
        else:
            return DivTautPolynomial([self.formal_ck(level,k) for k in range(0 ,codim+ 1 )])
        


    def replace_formal_by_symbolic(self, divTaut):
        for level in self.level_set:
            divTaut = divTaut.substitute( self.formal_chern(level)[ 1 :],                                         self.symbolic_chern(level)[ 1 :])
        return divTaut
    
    # Recursively computes the total Chern class as a symbolic expression
    @cached_method
    def symbolic_chern(self, level, avoiding_blowup = True):
        blowup = self.exact_stratum.blowup(self, level)
        
        
        # If we are only interested in the fundamental class
        # and every blowup center has the expected codim,
        # then we dont have to blowup
        if self.X.g[0] < 2  and self.profile == () and level ==0  and blowup.has_expected_codim() and avoiding_blowup:
            print('Avoiding blowup')
            # We only compute the class but not the normal bundle
            return DivTautPolynomial([0 ]), blowup.class_no_blowup()
        
        chern, _ = blowup.zero_locus_chern()
        if blowup.num_of_blowups == 0 :
            return chern
        else:
            blowdown = blowup.blowdown(chern)
            for profile in blowup.blownup_profiles:
                
                blowdown = profile.replace_formal_by_symbolic(blowdown)
        return blowdown
    
    
    @cached_method
    def ELG_chern(self, level):
        return self.exact_stratum.to_ELG(self.symbolic_chern(level)[- 1 ])
        

    
    
    def formal_ck(self, level, degree):
        if degree < 0  or degree > self.exact_stratum.exact_rank:
            return 0 
        if degree == 0 :
            return  1 
        
        profile_str = 'c' + ''.join(['_{}'.format(i) for i  in self])
        
        

#         latex_str = 'c_' + ''.join(['{}'.format(i) for i  in self])\
#                           + ''.join('_{}_{}'.format(-level,degree))
        
        
        return var(profile_str + ''.join('__{}_{}'.format(-level,degree)))
        
 
    # Determines the level map to an undegeneration of profile
    def level_adjustment(self, undegeneration, level):
        # No adjustment needed if both profiles are the same
        if self == undegeneration:
            return level
        
        #Divisorial degeneration?
        if len(undegeneration) + 1  == len(self):
            bic = (undegeneration-self)[0 ]
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
        
   
    # Returns the profile that (proper) degeneration of self obtained by splitting level i    
    def profile_splits_level(self, degeneration, i):
        
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
        if i > 0  or i < -len(self):
            
            raise ValueError("Level not in range")
        if len(self) == 0 :
            return True
        index = -i
        if  index- 1  >= 0 :
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
            if self.exact_stratum.dim - bd.total_codim >=0 :
            
                # A_{P+bic}^{[i]} is a divisor in D_P exactly if A_{P+bic}^{[i]}= D_{P+bic}
                if bd.codim == 0 : # Reminder: codim of exact boundary stratum A_P^{i} is within D_P     
                    building_divs[bic_index] = bd
                else:    
                    # Exact Boundary Strata of codim inside D_P > rank of the globally defined bundle have to be empty   
                    if bd.codim + 1  <= self.ranks[-level]:
                        building_set[bic_index] = bd
        return building_set, building_divs


       
# Tuple of symbolic expressions in the variables xi and D_i where i is the index of a BIC
# Representing a polynomial in divisorial classes
# We want DivTautPolynomial to be cachable, so need to work with tuples
class DivTautPolynomial(SageObject):
    
    def __init__(self, div_taut_poly):
        self._div_taut_poly = div_taut_poly
        if (type(self._div_taut_poly)) != tuple:
            try:
                self._div_taut_poly = tuple(self._div_taut_poly)
            except:
                print('Cannot convert to tuple')
                
        if type(self._div_taut_poly) != tuple:
            raise ValueError('DivTautPolynomial needs to work with tuples')
#             self._div_taut_poly = (div_taut_poly)
        if len(self.div_taut_poly) == 0 :
            self._div_taut_poly = (0 ,)
            
            
    
    @property
    def div_taut_poly(self):
        return self._div_taut_poly
    
    def tuple(self):
        return self.div_taut_poly
    
    def list(self):
        return list(self.tuple())
    
    def __eq__(self, other):
        if other == None:
            return False
        return self.div_taut_poly == other.div_taut_poly
    
    def __hash__(self):
        return hash(self._div_taut_poly)
        
       
        
    def __repr__(self):
        return str(self.div_taut_poly)
    
    def __getitem__(self, key):
        return self.div_taut_poly[key]
    

    def __add__(self, other):
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) + S(other.div_taut_poly)).list())
        
    
    def __mul__(self, other):
        if isinstance(other, sage.rings.integer.Integer):
            return DivTautPolynomial([coeff*other for coeff in self.div_taut_poly ])
        if not isinstance(other, DivTautPolynomial):
            
            raise ValueError(str(other)+'is of type' +str(type(other))+' and not a DivTautPolynomial')
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) * S(other.div_taut_poly)).list())
    
    def __rmul__(self, other):
        if isinstance(other, sage.rings.integer.Integer):
            return DivTautPolynomial([coeff*other for coeff in self.div_taut_poly ])
        
    def __sub__(self, other):
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) - S(other.div_taut_poly)).list())
        
        
    
    def __pow__(self, exponent):
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial(((S(self.div_taut_poly))**exponent).list())
        
        
    @property
    def codim(self):
        return len(self.list())
    
    
    def reciprocal(self):
        return DivTautPolynomial(tuple(reversed(self.div_taut_poly)))
    
    
    def substitute(self, old_variables, new_variables, verbose = False):
        if len(old_variables) != len(new_variables):
                show(old_variables, new_variables)
                raise ValueError("Number of variables is not the same")
        substitution_dict = {old_variables[i]: new_variables[i]  for i, _ in enumerate(old_variables)}
        if verbose:
            show(substitution_dict)
        substituted_list = []
        for term in self.div_taut_poly:
            substituted_list.append(term.subs(substitution_dict))
         
        return DivTautPolynomial(substituted_list)
    
    
    
    def proper_transform(self, var):
        S = SR['t']; (t,) = S._first_ngens(1)
        poly =  S(self.div_taut_poly)
        return DivTautPolynomial(poly(t= S(t+var)).list())
    
    def tensor(self, var):
        return self.reciprocal().proper_transform(var).reciprocal()
    
   
    
    
    # Inverts the DivTautPolynomial in the Chow ring CH^*(X)
    # dim - dimension of ambient space = maximal length of the inverse - 1
    def invert(self, dim):
        # Leading term must be non-zero
        assert self[0 ]!=0 
        Z = PowerSeriesRing(SR, 'z', default_prec=dim+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        assert Z(self.div_taut_poly) is not 0 , self.div_taut_poly
        inverted = ( 1 /Z(self.div_taut_poly)).list()
        return DivTautPolynomial(inverted)
    
    def divide(self, other, length):
        Z = PowerSeriesRing(SR, 'z', default_prec=length+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        inverted = ( 1 /Z(self.div_taut_poly)).list()
        return (self * other.invert(length)).trim(length)
        
    
    # Ensure that the polynomial has degree (i.e. codim) at most k
    def trim(self, length):
        # Codimension is greather than length?
        if length >=0  and length <= len(self.list())- 1 :
            #Then cut off the remainder
            return DivTautPolynomial(self.list()[0 :length+ 1 ])
        #Otherwise dont need to change anything
        else:
            return self
        
    
    
    # Makes sure the tuple has the right length
    # Since DivTautPolynomial should be hashable we are returning a new DivTautPolynomial
    def extend_to_length(self, target_len):
        new_poly =  self.tuple()[:target_len] + (0 ,)*(target_len - len(self.div_taut_poly))
        return DivTautPolynomial(new_poly)
     

class IteratedBlowup(SageObject):
    
    def __init__(self, exact_stratum, global_profile, level, base, blowup_centers):
        
        
      
        self.exact_stratum = exact_stratum
        self.X = exact_stratum.X
        self.base = base
        
        
        # List of ExactBoundaryStrata that are being blow up
        self.blowup_centers = blowup_centers
        self.center_index = { center : i for i, center in enumerate(self.blowup_centers)}
        
        # Ambient boundary stratum D_P
        self.global_profile = global_profile
        self.level = level
        
        self.num_of_blowups = len(self.blowup_centers)
        self.xi_at_level = {i: var('xi_'+ ''.join('{}'.format(-i))) for i in self.global_profile.level_set}
        
        self.blownup_profiles = [bd.profile for bd in self.blowup_centers]
        self.basic_exact_stratum = self.exact_stratum.basic_exact_stratum(self.global_profile, self.level)
        
        self.divisors = self.global_profile.building_set(self.level)[- 1 ]
        
        
        # Exceptional divisors
        self.E = {center : var("E_{}".format(i)) for i, center in enumerate(self.blowup_centers)}
 
        
    def __repr__(self):
        return "Blowup of ExactStratum " + self.exact_stratum.__repr__()         + "\n"+ "Global Profile: " + self.global_profile.__repr__()        + "\n"+ "Base:" + str(self.base.profile)+','+str( self.base.levels)        + "\n"+ "Blowup centers: " + str([(bd.profile, bd.levels) for  bd in self.blowup_centers])        + '\n'+ 'Divisors: ' + str(list(self.divisors.keys()))
    
    def codim(self, deg, undeg):
        assert deg.is_contained(undeg)
        return deg.total_codim - undeg.total_codim
   
    # Creates the iterated blow up
    def sub_blowup(self, exc_bd):
        assert exc_bd.is_contained(self.base), "ExactBoundaryStratum not a subspace of base"
        
        
        
        step_adjustment = {- 1 :- 1 }
        new_centers = []
        for i, center in enumerate(self.blowup_centers):
            
            if not center.are_disjoint(exc_bd):
                
                # Since we might remove some centers,
                # we need to keep track which step in the old blow up corresponds to which step
                # in the new blow up
                step_adjustment[i] = len(new_centers)
                new_centers.append(exc_bd + center)
            else:
                #We are omitting this step
                step_adjustment[i] = len(new_centers)- 1 
                
        
        
        # When replacing formal by symbolic chern classes we need to iterate through all possible profiles
        # that we have been blowing up
        self.blownup_profiles = self.blownup_profiles + [bd.profile for bd in new_centers]
        
        return IteratedBlowup(self.exact_stratum, self.global_profile, self.level, exc_bd, new_centers),                step_adjustment
    
#     # Computes the normal bundle of A_P,Q^{[I]} inside A_{P,S}^{[K]}
#     @cached_method
#     def normal_bundle(self, deg, undeg, mult = False):
      
#         nb = deg.formal_normal_bundle(undeg, mult = mult)
        
        
#         return nb, nb[-1]
    
    
    def has_expected_codim(self):
        for bd in self.blowup_centers:
            # All boundary strata have expected codim?
            if bd.total_codim == self.basic_exact_stratum.total_codim:
                return True
    
   
  
    
    
    
    # Decomposes the blowup of Y=undeg inside X=deg as
    # X -> I -> Y, 
    # Here I has the same underlying profie as X but potentially less levels,
    # where the normal bundle N(I/Y) is the normal bundle of a boundary stratum in the stratum
    # The boundary stratum intersects the base transversely, so the proper transform is simply a pullback
    # and the normal bundle N(X/I) is the product of normal bundles of proper transforms 
    # N(A_{P,Q}^{[i]}/ A_{P,Q}^{[J]})
    def aluffi_decomposition(self, deg, undeg, step):
        intermediate, remainder_profile, remainder_levels = deg.product_decomposition(undeg)
        remainder_nb = remainder_profile.formal_normal_bundle(mult = True)
        
        if intermediate != deg:
            
            new_blowup, new_steps = self.sub_blowup(intermediate)
            
            recursive_nb, recursive_fc = new_blowup.proper_transform(deg, new_steps[step])
            return remainder_nb * recursive_nb, remainder_nb[- 1 ]* recursive_fc
        else: 
            return remainder_nb, remainder_nb[- 1 ]
        
        
        
        
    
    
    
    
    # For the i-th building in the blowup, compute the proper transform
    # on the last step before the i-th building is blown up
    # First output is the total Chern class, second is the fundamental class with multiplicity
    @cached_method
    def building_transform(self, step):
        
        bd = self.blowup_centers[step]
        return self.proper_transform(bd, step- 1 )
    
    
    # Computes the Segre class of the proper transform of the i-th building
    # on the i-1-th blowup
    @cached_method
    def building_segre(self, step):
        #The dimension of the blow-up center, gives upper bound for the Length of the total Segre class
        dim = self.exact_stratum.dim - self.blowup_centers[step].total_codim
        total_chern, fc = self.building_transform(step)
       
        return [term* fc for term in total_chern.invert(dim)]
    
   
    
    
    # Returns the difference between total_transform of pullback and proper transform
    # Only implemented for the case that exc_bd is not contained in ANY blowup-center
    def final_blowup_difference(self, exc_bd):
        return self.proper_transform(exc_bd,len(step)- 1 )
        
    
    
    
    @cached_method
    # Returns the total Chern class and fundamental class of proper transform
    def proper_transform(self, exc_bd, step):
        if step == - 1 :
            nb = exc_bd.formal_normal_bundle(self.base, mult = True)
            return nb, nb[- 1 ] 
        else:
            nb_old, fc_old = self.proper_transform(exc_bd,step- 1 )
            nb_diff, fc_diff = self.proper_transform_diff(exc_bd, step)
            return  nb_old+nb_diff, fc_old+fc_diff 
        
        
        
    
    
    # Returns the difference between the total Chern class of the proper transform and the pullback
    # on the i-th step of the blow up
    @cached_method
    def proper_transform_diff(self, exc_bd, step):
        
        if exc_bd in self.blowup_centers:
            assert step < self.blowup_centers.index(exc_bd)
        
        
        # No difference in the first step
        if step == - 1 :
            return DivTautPolynomial([0 ]),0 
        
        
        assert step >=0  and step< self.num_of_blowups, str(step) + 'is no in range ' + str(self.num_of_blowups)
        # Locate blowup_center
        center = self.blowup_centers[step]
        base = self.base
        
        # exc_bd is disjoint from blowup center?
        if center.are_disjoint(exc_bd):
            # No difference
            return DivTautPolynomial([0 ]),0 
        else:
            intersection = exc_bd + center
            
             # Codimensions of X,Y,W in Z, in that order
            codims = [self.codim(intersection, base), self.codim(exc_bd, base), self.codim(center, base)]
            codimX, codimY, codimW = codims
            # Check if the intersection is transversal:
            if codimX == codimY+codimW:
                # No difference
                return DivTautPolynomial([0 ]),0 
            
            #We are now setting up aluffis formula
            # Y - variety whose proper transform we are trying to compute
            # W - center of blowup
            # X - intersection of variety and blowup center
            # Z - ambient space
            cXW, fXW = self.aluffi_decomposition(intersection, center, step- 1 )
            cXY, _ = self.aluffi_decomposition(intersection, exc_bd, step- 1 )
            cWZ, _ = self.aluffi_decomposition(center, base, step- 1 )
            cherns = [cXY, cXW, cWZ]
            E = self.E[self.blowup_centers[step]]
            

            # It remains to change the multiplicity of the last term
            # We use the conventation that 
            
            
            total_chern, fundamental_class = self.aluffi_formula(cherns, codims, E, fXW) 
            
            return total_chern, fundamental_class
     
    
    # If all blow up centers have exptected codim,
    # we can obtain the class of A_P^{[i]} simply via the top chern Class of the globall defined bundle
    # - boundary divisors - classes of A_{P,BIC}^{[i]}
    def class_no_blowup(self):
#         print("Computing the class of " + str(self.basic_exact_stratum) +" without blowing up.")
        _, building_divs = self.global_profile.building_set(self.level)
        div_class = -self.exact_stratum.formal_xi_at_level(self.global_profile, self.level)
        
        
        
        for bic_index in building_divs.keys():
            bic = self.exact_stratum.profile((bic_index,))
            div_class -= bic.formal_normal_bundle(mult = True)[- 1 ]
            
#         print("Line bundle class: ", div_class)
        rank = self.global_profile.ranks[-self.level]
        
        # Total Chern of the proper transform of A_P^{[0]}
        chern = DivTautPolynomial([ 1 ,div_class])**rank
        fc = chern[- 1 ]
#         print("Fundamental class: ", fc)

        diff = sum(exc_bd.formal_normal_bundle(self.base, mult  = True)[- 1 ] for exc_bd in self.blowup_centers)
#         for exc_bd in self.blowup_centers:
            
#             print("Difference for BIC ",exc_bd.profile[0], ' is :',\
#                   exc_bd.formal_normal_bundle(self.base, mult  = True))
#             print('Mult: ',self.X.bics[exc_bd.profile[0]].ell )
        return fc - diff
        
        
        
        
    
    @cached_method
    def zero_locus_chern(self):
        _, building_divs = self.global_profile.building_set(self.level)
        div_class = -self.exact_stratum.formal_xi_at_level(self.global_profile, self.level)
        
        
        
        for bic_index in building_divs.keys():
            bic = self.exact_stratum.profile((bic_index,))
            div_class -= bic.formal_normal_bundle(mult = True)[- 1 ]
        for exceptional_index in self.blowup_centers:
            div_class -= self.E[exceptional_index]

        rank = self.global_profile.ranks[-self.level]
        
        # Total Chern of the proper transform of A_P^{[0]}
        chern = DivTautPolynomial([ 1 ,div_class])**rank
        fc = chern[-1]
       
        # In order to obtain the total Chern class of the Pullback instead
        # we neeed to compute the difference using Aluffis formula
        exc_bd = self.basic_exact_stratum
        for step in range(self.num_of_blowups):
            difference_ch, difference_fc =  self.proper_transform_diff(exc_bd,step)
            chern -= difference_ch
            
            fc -= difference_fc
                                       
        
        
        
        
        return chern, fc
    
    def zero_locus_blowdown(self):
            return self.blowdown(self.zero_locus_chern()[0]), self.blowdown(self.zero_locus_chern()[-1])
    
    
    
    
    def blowdown(self, divTaut):
        result = divTaut
        reverse_centers = self.blowup_centers[::-1]
        for center in reverse_centers:
            step = self.blowup_centers.index(center)
            result = self._blowdown_step(step,result)
        return result
        
        
        

    # Computes the pushforward of a DivTautpolynomial on the i-th blowup
    # Should not be called by it self
    def _blowdown_step(self, step, divTaut):
        
        if type(divTaut) is DivTautPolynomial:
            return DivTautPolynomial([self._blowdown_step(step, term) for term in divTaut.list()])
        
        
        center = self.blowup_centers[step]
        E =  self.E[center]
        
        #Codim within base
        codim = center.total_codim - self.base.total_codim
        dim = self.exact_stratum.dim - center.total_codim
        segre = []
        
        def segre_k(k):
            if k<0  or k > dim:
                return 0 
            else:
                
                #Only compute the Segre class if necessary
                segre = self.building_segre(step)
                return segre[k]
        
 

        # Use formula (5.2) from the paper for the pushforward

        coeff = [ (-1)**(k-1 )* aj* segre_k(k - codim) for aj, k in divTaut.coefficients(E)]

        # The constant term stays the same
        coeff[0] = divTaut.coefficients(E)[0][0]
        return sum(coeff)
    
 
    ######
    # Aluffi logic
    ########
    
    
    def aluffi_formula(self, cherns, codims, E, fXW):
        cXY, cXW, cWZ = cherns
        total_chern = self.aluffi_normal_bundle(cherns, codims, E)
        fundamental_class = self.aluffi_proper_transform(cherns, codims, E, fXW)
       
        return total_chern, fundamental_class
    
    def aluffi_normal_bundle(self, cherns, codims, E):
        codimX, codimY, codimW = codims
        cXY, cXW, cWZ = cherns
        d = codimY 
        e = codimX - codimW
        excess = cWZ.divide(cXY, d-e)
       
        return self.aluffi_formal_tensor(cXW, excess, E).trim(codimY)
    
    def aluffi_proper_transform(self, cherns, codims, E, fXW):
        codimX, codimY, codimW = codims
        cXY, cXW, cWZ = cherns
        
        # To obtain the fundamental class of Xbar we need to consider the correct multiplicity
        Xbar = fXW
        d = codimY 
        e = codimX - codimW
        
        excess = cWZ.divide(cXY, d-e)
        segre = [(-E)**i for i in range(1,d-e+1)]
        
        return - sum(excess[d-e-k]*Xbar*segre[k-1] for k in range(1,d-e+1))



    def aluffi_formal_tensor(self, A,B, E):
        transform = B.reciprocal().proper_transform(-E).reciprocal()

        result =  (A * transform).list()
        #Now need remove all all terms without E coefficient
        for i, term in enumerate(result):
            result[i] = sum( coeff* E**power for coeff,power in term.coefficients(E)[ 1 :])

       
        return DivTautPolynomial(result)

    
   

       

  
            
            
    
            
    
        
    
    

