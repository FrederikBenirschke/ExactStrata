import sage

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

from sage.modules.free_module import FreeModule  # pylint: disable=import-error
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from exactstrata.profile import *
from exactstrata.exactboundarystratum import *
from exactstrata.iteratedblowup import *
from exactstrata.divtautpolynomial import *
from sage.calculus.var import var

from exactstrata.helper import *





#######################################################
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
        self.residue_rank = sum(max(0, len(signature.pole_ind)-1) for signature in self.X._sig_list_H)
        # The rank of the relative period local system
        self.torsion_rank = len(self.bundle)
        self._profiles = None
        
        self._bics = None
        
        
        

        self._boundary_strata = None
       
    # Shortcut for the profile () corresponding to the ambient GeneralisedStratum   
    @property    
    def empty_profile(self):
        return Profile(self,())
 
        
   
    
    @property
    def dim(self):
        '''Returns the dimension of the ambient GeneralisedStratum'''
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
        return 2*sum(self.X.g_H) 
    
    
    # Lists all BICs such that the restriction of the residue subspace containing i many marked poles
    # to the two level graph is the same as the restriction containing only i-1 marked poles    
    def res_bics(self, index):
        assert index >=1 and index<= self.residue_rank
        return [bic for bic in self.bics if self.residue_ranks_bic(bic)[index]== self.residue_ranks_bic(bic)[index-1]]
    
    
    # Same as res_bics but instead for the bundle of relative periods
    def torsion_bics(self, index):
        assert index >=1 and index<= self.torsion_rank
        return [bic for bic in self.bics if self.torsion_ranks_bic(bic)[index]== self.torsion_ranks_bic(bic)[index-1]]
        
    # The class of the locus of residueless differentials
    # as SymbolicExpression
    def res_class(self):
        result = 1
        for index in range(1, self.residue_rank+1):
            result *= var('xi')+sum(self.X.bics[bic].ell * self.bic_variables[bic] for bic in self.res_bics(index))
        
        if self.to_ELG(result) == self.X.ZERO:
            return 0
        return (-1)**(self.residue_rank) * result
    
    # The class of the locus of  lambda-exact differentials inside the stratum of exact differentials
    # as SymbolicExpression
    def torsion_class(self):
        result = 1
        for index in range(1,self.torsion_rank+1):
            result *= var('xi')+sum(self.X.bics[bic].ell * self.bic_variables[bic] for bic in self.torsion_bics(index))
        return (-1)**(self.torsion_rank) * result
    
    
    
    def res_torsion_class(self):
        '''Product of class of the locus of residueless differentials and the class 
        of the locus of  lambda-exact differentials inside the stratum of exact differentials
        as SymbolicExpression
        In general, this class does not have geometric meaning
        but for g=0 this is the class of the stratum of lambda-exact differentials.'''
        return self.res_class() * self.torsion_class()
    
    
    
    def torsion_ranks_bic(self, bic_index):
        '''Returns the ranks of a total flag for the local system of relative homology to a two level graphs.'''
        torsion_ranks = []
        
        for index in range(len(self.bundle)+1):
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
        for index in range(len(all_poles)+1):
            index_poles = bic_poles[:index]
#             show('Index: ', index, 'Index poles:', index_poles)
#             bic_poles = [leg for leg in bic_poles if BIC.LG.levelofleg(leg)==0 ]
                
            residue_conditions = 0
            for v, _ in enumerate(LG.genera):
                if LG.levelofvertex(v) == 0:



                    poles_on_vertex = len([l for l in LG.legsatvertex(v) if LG.orderatleg(l) < 0])
                    poles_from_list = len([pole for pole in index_poles if LG.vertex(pole) == v])

                    # Need to take  into account that residues add up to zero
                    if poles_on_vertex == poles_from_list:
                        residue_conditions += max(0,poles_from_list -1)
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
            show(pole_list)
            residue_conditions = 0
            for i, sig in enumerate(X._sig_list_H):
                poles_on_vertex = len(sig.pole_ind)
                poles_from_list = len([pole for pole in pole_list if (pole[0]==i and pole[1] in sig.pole_ind)])

                # Need to take  into account that residues add up to zero
                if(poles_on_vertex==poles_from_list):
                    residue_conditions += max(0,poles_from_list -1)
                else:
                    residue_conditions += poles_from_list

            return residue_conditions

        return [residue_rank(index) for index in range(len(all_poles)+1)]

    
    
                
    
    @cached_method
    def profile(self, profile_tuple):
        ''' Creates a profile from a tuple. The entries of the tuple are the indices of the BICs.'''
        return Profile(self, profile_tuple)
    
    # Returns all profiles (as Profiles) of the underlying GeneralisedStratum, ordered by ascending dimension
    @property
    def profiles(self):
        if self._profiles == None:
            self._profiles = sorted([Profile(self, pf) for pf in self.profile_list],key=len, reverse=True)
        return self._profiles
    
    
    @property
    def profile_list(self):
        '''Returns a list of all profiles (as tuples).
        It is used for checking that a profile is valid.'''
        return [pf for pf_list in self.X.lookup_list for pf in pf_list]
        
    
    
    
   
    @property
    def bics(self):
        '''  Returns a list of all indices of BICs, respecting the partial order on two-level graphs.'''
        if self._bics == None:
            old_list = list(range(len(self.X.bics.copy())))
            new_list = []
            
            # We need to reorder the list, to ensure that it preserves the partial order on level graphs
           
            while len(old_list) > 0:
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
    

  
            
    def basic_exact_stratum(self, profile, i):
        ''' Returns the boundary stratum A_P^{[i]}.
        Here P is a profile and i(int) denotes the level where the differential is exact.'''
#         assert type(i) is sage.rings.integer.Integer, str(i) + " is not an integer but has type "+ str(type(i))
        return self.exact_boundary_stratum(profile, (i,))
    
   
    @cached_method
    def exact_boundary_stratum(self, profile, levels):
        ''' Creates an ExactBoundaryStratum from a profile (Profile or tuple)
        and a  list of levels I,  corresponding to subscheme A_P^{[I]} of exact differentials 
        exact on each level in I.'''
        profile = Profile(self, profile)
        return ExactBoundaryStratum(profile, levels)
    

            
    
   
    @cached_method
    def D_lg(self, bic):
        ''' Returns the ELGTautClass of a BIC (as integer).'''
        return self.X.additive_generator(((bic,),0))
    
    
    # Returns dictionary converting from a SymbolicExpression 'D_{bic}' to the ELGTautClass
    # of the corresponding boundary divisor
    @property
    def divisor_classes(self):
        return {var('D_{}'.format(bic)) : self.D_lg(bic) for bic in self.bics}
    
    
    # Returns the xi class at a level of a boundary stratum D_P as a SymbolicExpression.
    # The result is a sum of 'xi' and 'D_{bic}' for different BICs
    @cached_method    
    def formal_xi_at_level(self, profile, level):
        assert level in profile.level_set, show(level, profile.level_set)
        if len(profile) == 0 or  level==0:
            return var('xi')
        if len(profile) == 1:
            bic_index = profile[0]
            xi_bot = var('xi') + self.X.bics[bic_index].ell* self.bic_variables[bic_index] 
            xi_bot += sum([self.X.bics[nbic].ell * self.bic_variables[nbic] \
                           for nbic in self.bics if self.X.lies_over(nbic, bic_index)]) 
            return xi_bot
            
           
        else:
            return self.formal_xi_at_level(Profile(self, (profile[-level-1],)), -1)
        
       
        

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
        return str([[sig.sig for sig in self.X._sig_list_H], self.bundle])  
    
    
      ######################################################### 
    # -------- Blow up algorithm --------
    #########################################################
   
    
    
    # Constructs an iterated blowup. 
    # Given a boundary stratum D_P and a level i
    # this is an iterated blow up along A_{P,Q}^{[i]} where P+Q is a degeneration of P
    # obtained by splitting the i-th level.
    # (Except for boundary strata such that A_{P,Q}^{[i]})
    # is a divisor in D_P)
    # 
    def blowup(self, amb_profile , i):

        amb_profile = self.profile(amb_profile)

        assert i in amb_profile.level_set
        # Step 1) Find the building set, i.e. all two-level graphs that split level i

        building_set, building_divs = amb_profile.building_set(i)

        blowup_centers = list(building_set.values())
        blowup = IteratedBlowup(self, amb_profile, i, self.exact_boundary_stratum(amb_profile,()), \
                                blowup_centers)

        return blowup
    
    #--------------------------------------------------------------------------
    # Chern polynomial logic
    #--------------------------------------------------------------------------
    
    # Returns the fundamental class of a stratum of lambda-exact differentials.
    # The result is a SymbolicExpression in terms of 'xi' and 'D_{bic}'.
    def exact_stratum_class(self, avoiding_blowup=True):
        # We first reduce to the case without relative periods
        if len(self.bundle) > 0:
            new_exact_stratum = ExactStratum(self.X, [])
            return new_exact_stratum.exact_stratum_class(avoiding_blowup = avoiding_blowup) * \
        self.torsion_class()
        
        # Now we deal with all residue equations
        else:
            res_cl = self.res_class()
            if res_cl == 0:
                return self.X.ZERO
            else:
                # Finally we need to compute the class corresponding to the remaining ansolute periods
                return self.res_class() * \
            self.symbolic_class(self.empty_profile, 0, avoiding_blowup= avoiding_blowup)
        
    
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
        
        result =  0
        
        if type(symbolic_class) is sage.rings.integer.Integer:
            return symbolic_class *  X.ONE
        
        if type(symbolic_class) is not sage.symbolic.expression.Expression:
            try:
                symbolic_class = SR(symbolic_class)
            except:
                raise ValueError('Cannot convert {}'.format(type(symbolic_class)) + 'to ELG.')
        if len(symbolic_class.variables())==0:
            return symbolic_class * X.ONE
        else:
            # Convert the first variable that appears in the symbolic expression.
            variable = symbolic_class.variables()[0]

            AG = self.formal_to_AG_dic()[variable]

            for coefficient, power in symbolic_class.coefficients(variable):
                # Each coefficient is a SymbolicClass with one variable less
                # So we can convert recursively
                result += AG.__pow__(power) * self.to_ELG(coefficient)


        return result
    
    
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
        
        num_marked_points =  sum([len(sig.sig) for sig in self.X._sig_list_H]) #number of marked points

        point_module = FreeModule(QQ, num_marked_points)
        for vector in self.bundle:
            if(sum(vector)!=0):
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
        marked_topl = [l for i,_ in enumerate(ELG.LG.genera)\
                       if ELG.LG.levelofvertex(i)==0 for l in ELG.LG.legsatvertex(i) ]
        num_marked_topl = len(marked_topl)
        rel_module = FreeModule(QQ, num_marked_topl)
        topl_dict = {point : index for index, point in enumerate(marked_topl)}
        return marked_topl, num_marked_topl, rel_module, topl_dict
    
    
    
    # # For a BIC finds all relative cycles in the image of the absolute homology under the specialization map 
    # # to the top level
    # def absolute_bundle_basis(self, ELG):
        
    #     basis= ELG.LG.stgraph.cycle_basis()
    #     marked_topl , num_marked_topl, _,topl_dict = self.top_level_homology(ELG)
        
    #     cycle_basis =[]
    #     for path in basis:
    #         vec = [0] * num_marked_topl
    #         for i,e in  enumerate(path):
    #             leg = sort_edge(ELG,ELG.LG.edges[i])[0]
    #             vec[topl_dict[leg]]= e
    #         cycle_basis.append(vec)
    #     return cycle_basis



    # For a BIC finds all relative cycles in the image of the absolute homology under the specialization map 
    # to the top level
    def absolute_bundle_basis(self, ELG):
        
        basis= ELG.LG.stgraph.cycle_basis()
        marked_topl , num_marked_topl, _,topl_dict = self.top_level_homology(ELG)
        
        cycle_basis =[]
        for path in basis:
            vec = [0] * num_marked_topl
            for i,e in  enumerate(path):
                leg = ELG.LG.edges[i][0]
                vec[topl_dict[leg]] = e
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
            if sum(vector)!=0:
                raise ValueError("Not a relative homology class")
            new_vector = [0] * num_marked_topl
            for i, x in enumerate(vector):
                #we choose x * path from i-th marked point to some marked point
                #ensure that the marked point is on the same connected component
                start = ELG.dmp_inv[marked_to_component(self.X)[i]]
                #end = ELG.dmp_inv[(0,0)]
                end = ELG.dmp_inv[(marked_to_component(self.X)[i][0],0)]
                path = path_between_markings(ELG, start, end)
                #now go through the path and add x*(j+1-th point - j-th point) but only for the points in pts,
                #since the other points correspond to marked points in lower levels
                for j in range(len(path)-1):

                    a,b = path[j],path[j+1]
                    if ELG.LG.levelofleg(a)==0 and ELG.LG.levelofleg(b)==0:
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
        for i in range(1,len(ELG.dlevels)):
            squish = ELG.delta(i)
#             squish_index = self.X.bics.index(squish)
            rank_i = self.bundle_to_top(squish).total_rank - self.bundle_to_top(squish).residue_rank
            ranks.append(rank_i - sum(ranks))
            
        # The rank on the bottom level is the rank of V minus the rank on all higher levels
        ranks.append(self.exact_rank - sum(ranks))
        return ranks
                 
