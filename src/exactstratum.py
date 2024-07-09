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







class ExactStratum(SageObject):
    ''' #######################################################
        # A (generalized) stratum of exact differentials.
        # An ExactStratum is uniquely identified by the following information:
        #
        #   * X: a GeneralisedStratum, corresponding to the ambient space 
        #   * bundle: list of relative homology cycles, i.e. [R_1,...,R_n] where each R_l is
        #      a list of integers, corresponding to marked zeros of the differentials. The entries in R_l need 
        #      to add to zero and determine which relative periods are zero
        # 
        #  The main functionalities are
        #  exact_stratum_class(): Returns the fundamental class of the locus of lambda exact differentials 
        #                         as a SymbolicExpression representing a divisorial tautological class
        # 
        #  to_ELG(symbolic_class): Takes  a SymbolicExpression representing a divisorial tautological class
        #                          and returns the corresponding ELGTautclass.
        #
        # pushforward(ELGTautClass): Takes an ELGTautClass and computes the push forward to the moduli space Mg,n.'''
    
    
    
    def __init__(self, X, bundle):
        """
        Initializes an instance of the class with the given arguments.

        Args:
            X (GeneralisedStratum): The ambient space corresponding to the ExactStratum.
            bundle (list): A list of relative homology cycles, where each cycle is represented as a list of integers.
                Each integer corresponds to the marked zeros of the differentials. The entries in each cycle need to add to zero
                and determine which relative periods are zero.

        Initializes the following instance variables:
            X (GeneralisedStratum): The ambient space corresponding to the ExactStratum.
            _dim (None): The dimension of the ExactStratum.
            bundle (list): The bundle of relative homology cycles.
            _total_rank (None): The total rank of the ExactStratum.
            residue_rank (int): The rank of the residue local system.
            torsion_rank (int): The rank of the relative period local system.
            _profiles (None): The profiles of the ExactStratum.
            _bics (None): The boundary intersection cycles of the ExactStratum.
            _boundary_strata (None): The boundary strata of the ExactStratum.
        """
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
       
    
    @property    
    def empty_profile(self):
        """
        Returns the empty profile of the ExactStratum.

        :return: An instance of the Profile class with an empty profile.
        :rtype: Profile
        """
       
        return Profile(self,())
 
        
   
    
    @property
    def dim(self):
        """
        Returns the dimension of the ambient GeneralisedStratum.

        :return: The dimension of the GeneralisedStratum.
        :rtype: int
        """
       
        if self._dim == None:
           
            self._dim = self.X.dim()
        return self._dim
    
    
    @property
    def total_rank(self):
        """
        Returns the total rank of the ExactStratum.

        The total rank is the sum of the exact rank, residue rank, and torsion rank.
        The total rank is equal to the codimension of the stratum of exact differentials within the GeneralisedStratum.

        :return: The total rank of the ExactStratum.
        :rtype: int
        """
       
        if self._total_rank == None:
            
            # The rank is the dimension of absolute cohomology of the punctured(!) surface
            # Plus the number of linearly independent relative equations
            self._total_rank = self.exact_rank +self.residue_rank + self.torsion_rank
        return self._total_rank
    
    
    
    @property
    def exact_rank(self):
        """
        Returns the exact rank of the ExactStratum.

        The exact rank is the sum of twice the number of elements in the list `self.X.g_H`. 
        This represents the rank of the bundle of absolute cohomology (not including residues).

        :return: The exact rank of the ExactStratum.
        :rtype: int
        """
      
        return 2*sum(self.X.g_H) 
    
    
    def res_bics(self, index):
        """
        Returns a list of boundary intersection cycles (BICs) from the `self.bics` list
        with the following property:    
        The restriction of the residue subspace containing the first i marked poles
        to the two level graph with index BIC is the same as the restriction containing only the first i-1 marked poles.

        Parameters:
            index (int): The index of the residue rank to compare.

        Returns:
            list: A list of BICs that have the same residue rank at the specified index.
        """
          
        assert index >=1 and index<= self.residue_rank
        return [bic for bic in self.bics if self.residue_ranks_bic(bic)[index]== self.residue_ranks_bic(bic)[index-1]]
    
    
   
    def torsion_bics(self, index):
        """
        Lists all Boundary Intersection Cycles (BICs) that satisfy the property where 
        the restriction of the subspace of relative periods containing the first `index` relative periods
        to the two-level graph with index `BIC` is equal to the restriction containing only the first `index` relative periods.
        
        This function is equivalent to `self.res_bics(index)`, but specifically for the bundle of relative periods.

        Parameters:
            index (int): The index of the relative period to compare.

        Returns:
            list: A list of BICs that meet the specified property.
        """
        assert index >=1 and index<= self.torsion_rank
        return [bic for bic in self.bics if self.torsion_ranks_bic(bic)[index]== self.torsion_ranks_bic(bic)[index-1]]
        
   
    def res_class(self):
        """
        Calculates the class of the locus of residueless differentials.

        Returns:
            int or SymbolicExpression: The class of the locus of residueless differentials.
                Returns 0 if the class is zero.
                Returns a SymbolicExpression representing the class otherwise.

        Description:
            This function calculates the class of the locus of residueless differentials.
            It iterates over the residue indices and calculates the class for each index.
            The class is represented as a product of the variable 'xi' and the sum of the product
            of the 'ell' attribute of the BICs and the corresponding bic_variable.

            If the class is zero, the function returns 0.
            Otherwise, it returns the negative of the class multiplied by the residue rank.

      
        """
        
        result = 1
        for index in range(1, self.residue_rank+1):
            result *= var('xi')+sum(self.X.bics[bic].ell * self.bic_variables[bic] for bic in self.res_bics(index))
        
        if self.to_ELG(result) == self.X.ZERO:
            return 0
        return (-1)**(self.residue_rank) * result
    
    
    def torsion_class(self):
        """
        Calculates the class of the locus of lambda-exact differentials inside the stratum of exact differentials
        as a SymbolicExpression.

        Returns:
            SymbolicExpression: The class of the locus of lambda-exact differentials.

        Description:
            This function calculates the class of the locus of lambda-exact differentials by iterating over the torsion rank.
            For each index, it multiplies the variable 'xi' by the sum of the product of the 'ell' attribute of the BICs
            and the corresponding bic_variables. The result is then multiplied by (-1) raised to the power of the torsion rank.

        """
      
        result = 1
        for index in range(1,self.torsion_rank+1):
            result *= var('xi')+sum(self.X.bics[bic].ell * self.bic_variables[bic] for bic in self.torsion_bics(index))
        return (-1)**(self.torsion_rank) * result
    
    
    
    def res_torsion_class(self):
        '''Returns the product of class of the locus of residueless differentials and the class 
        of the locus of  lambda-exact differentials inside the stratum of exact differentials
        as SymbolicExpression.
        In general, this class does not have geometric meaning
        but for g=0 this is the class of the stratum of lambda-exact differentials.'''
        return self.res_class() * self.torsion_class()
    
    
    
    def torsion_ranks_bic(self, bic_index):
        '''Returns the ranks of a total flag for the local system of relative homology on a two level graph.'''
        torsion_ranks = []
        
        for index in range(len(self.bundle)+1):
            sliced_bundle = self.bundle[:index]
            
            exc = ExactStratum(self.X, sliced_bundle)
            
            bund = len(exc.bundle_to_top(self.X.bics[bic_index]).bundle)
            
            torsion_ranks.append(bund)
        return torsion_ranks
    
    
    def residue_ranks_bic(self, bic_index):
        """
        Calculates the ranks of a total flag for the local system of residues on a two-level graph.

        Parameters:
            bic_index (int): The index of the BIC in the list of BICs.

        Returns:
            list: A list of integers representing the ranks of the total flag.
        """
       
        
        X = self.X
        all_poles = X._polelist
        BIC = self.X.bics[bic_index]
        bic_poles = [BIC.dmp_inv[pole] for pole in all_poles]

        
        LG = BIC.LG
        res_list = []
        for index in range(len(all_poles)+1):
            index_poles = bic_poles[:index]           
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

       
    
    
    
    
    def _residue_ranks(self,pole_list):
        ''' Returns a list of integers, where the i-th entry is the number of independent residue cycles for the first i marked poles.
        Which poles are chosen depends on the ordering of the signature of the underlying GeneralisedStratum.'''
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
    
    
    @property
    def profiles(self):
        """
        Returns a sorted list of all the profiles of the underlying GeneralisedStratum, ordered by ascending dimension.

        Returns:
            list: A list of Profile objects, sorted by the length of the profile in descending order.
        """
      
        if self._profiles == None:
            self._profiles = sorted([Profile(self, pf) for pf in self.profile_list],key=len, reverse=True)
        return self._profiles
    
    
    @property
    def profile_list(self):
        """
        Returns a list of all profiles (as tuples) of the underlying GeneralisedStratum.
        
        This property generates a list of all profiles by iterating over the lookup_list of the GeneralisedStratum X.
        Each element in the lookup_list is a list of tuples, where each tuple represents a profile.
        The method then flattens the list of tuples into a single list and returns it.
        
        Returns:
            list: A list of all profiles (as tuples) of the underlying GeneralisedStratum.
        """
       
        return [pf for pf_list in self.X.lookup_list for pf in pf_list]
        
    
    
    
   
    @property
    def bics(self):
        """
        Returns a list of all indices of BICs, respecting the partial order on two-level graphs.
        
        This property generates a list of all indices of BICs by reordering the list of indices based on the partial order on two-level graphs. The partial order is determined by the `lies_over` method of the `X` object.
        
        Returns:
            list: A list of all indices of BICs, respecting the partial order on two-level graphs.
        """
        
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
        Here P is a profile (tuple) and i (int) denotes the level where the differential is exact.'''
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
    
    
    
    @property
    def divisor_classes(self):
        ''' Returns a dictionary converting a SymbolicExpression 'D_{bic}' to the ELGTautClass
        of the corresponding boundary divisor.'''
        return {var('D_{}'.format(bic)) : self.D_lg(bic) for bic in self.bics}
    
    
    
    @cached_method    
    def formal_xi_at_level(self, profile, level):
        '''  Returns the xi - class at a given level of a boundary stratum D_P as a SymbolicExpression.
        The result is a sum of 'xi' and 'D_{bic}' for different BICs.'''
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
        
       
        

   
    @cached_method
    def formal_to_AG_dic(self):
         ''' Returns a dictionary converting formal xi classes and boundary classes 'D_{bic}' to ELGTautClasses.'''
        return self.divisor_classes | {var('xi'):self.X.xi }
                   
    @property
    def bic_variables_inv(self):
        ''' Returns a dictionary converting a SymbolicVariable to the index in self.bics of the corresponding BIC.'''
        return { var('D_{}'.format(i)): i for i,_ in enumerate(self.bics)} 
    
    @property 
    def bic_variables(self):
        """
        Returns a dictionary that maps the values of `bic_variables_inv` to their corresponding keys.

        Returns:
            dict: A dictionary mapping the values of `bic_variables_inv` to their corresponding keys.
        """
        return {value: key for key, value in self.bic_variables_inv.items()}
 
    def __repr__(self):
        """
        Returns a string representation of the object.

        Returns:
            str: A string representation of the object.
        """
        return str([[sig.sig for sig in self.X._sig_list_H], self.bundle])  
    
    
      ######################################################### 
    # -------- Blow up algorithm --------
    #########################################################
   
    
    
    
    def blowup(self, amb_profile , i):
        """ 
        Constructs an iterated blowup.
        
        Args:
            D_P (object): The boundary stratum.
            i (int): The level to split.
        
        Returns:
            IteratedBlowup: The iterated blowup object.
        
        Description:
            Given a boundary stratum D_P and a level i, 
            this function constructs an iterated blow up along all subvarieties A_{P,Q}^{[i]},
            where P+Q is a degeneration of P obtained by splitting the i-th level.
            (Ignores boundary strata such that A_{P,Q}^{[i]} is a divisor in D_P since here the blow-up does not do anything).
        
        is a divisor in D_P since here the blow-up does not do anything).
        """

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
    # In the next section we implement basic functionality for dealing with Chern classes in the Chow ring
    # of the moduli space of multi-scale differentials.
    
    
    def exact_stratum_class(self, avoiding_blowup=True):
        '''
        Returns the fundamental class of a stratum of lambda-exact differentials.

        Args:
            avoiding_blowup (bool, optional): Flag to indicate if avoiding blowup is enabled. Defaults to True.

        Returns:
            SymbolicExpression: The fundamental class of the stratum as a SymbolicExpression in terms of 'xi' and 'D_{bic}'.

        Description:
            The function returns the fundamental class of a stratum of lambda-exact differentials as a SymbolicExpression in terms of 'xi' and 'D_{bic}'.
            First it reduces to the case without relative periods.
            Then, it computes the class corresponding to the remaining absolute periods and returns the result.
        '''
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
        
    
    
    def symbolic_class(self, profile, level, avoiding_blowup = True):
        """
        Returns the fundamental class of A_P^{[i]} inside D_P.
        The result is computed inside the stratum of residueless differentials.
        If one wants the class inside the GenerlisedStratum, this needs to be multiplied by
        self.res_class().

        Args:
            profile (Profile): The profile of the class.
            level (int): The level of the class.
            avoiding_blowup (bool, optional): Whether to avoid blowup in the computation. Defaults to True.

        Returns:
            SymbolicClass: The fundamental class of A_P^{[i]} inside D_P.

        """
        profile = Profile(self, profile)
        
        # The class of the exact stratum is the last coefficient in the Chern polynomial
        return profile.symbolic_chern(level, avoiding_blowup=avoiding_blowup)[-1]
        


    
    def to_ELG(self, symbolic_class):
        """
        Converts a SymbolicExpression to an ELGTautClass.

        Args:
            symbolic_class (sage.symbolic.expression.Expression or int): The symbolic expression to convert.

        Returns:
            ELGTautClass: The converted ELGTautClass.

        Raises:
            ValueError: If the input is not a symbolic expression or an integer.

        Note:
            The symbolic expression should be a polynomial in 'xi' and 'D_{bic}'.

        """
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
    
    
    .
    def pushforward(self, fund_class):
              """
              Pushes forward a  class to a product taut class and returns the resulting taut class.
      
              Args:
                  fund_class (DivTautClass): A  DivTautClass.
      
              Returns:
                  TautClass: A taut class.
              """
        return fund_class.to_prodtautclass().pushforward()

    
  
    #--------------------------------------------------------------------------
    # Vector bundle logic
    #--------------------------------------------------------------------------
    
  
    def bundle_basis(self):
        """
        Compute the basis for the space of relative cycles.

        Returns:
            list: A list of vectors representing the basis for the space of relative cycles.

        Raises:
            ValueError: If the bundle contains a non-zero vector, indicating an invalid relative homology class.
        """
        
        num_marked_points =  sum([len(sig.sig) for sig in self.X._sig_list_H]) # Number of marked points

        point_module = FreeModule(QQ, num_marked_points)
        for vector in self.bundle:
            if(sum(vector)!=0):
                raise ValueError("Not a valid relative homology class")
        return list(point_module.submodule(self.bundle).basis())

  
    
    def bundle_to_top(self, BIC):
        """
        Compute the top level homology of an EmbeddedLevelGraph (ELG) and return an ExactStratum object supported on the
        underlying GenerlisedStratum of the top level.

        Parameters:
            BIC (EmbeddedLevelGraph): The EmbeddedLevelGraph representing a BIC.

        Returns:
            ExactStratum: An ExactStratum object representing the top level homology of the BIC.

        Raises:
            ValueError: If the input BIC is not a valid BIC.
        """
        
        if not BIC.is_bic():
            raise ValueError("Level graph is not a BIC")
        marked_topl ,num_marked_topl, rel_module, topl_dict = self.top_level_homology(BIC)
        total_relations = self.relative_bundle_basis(BIC) + self.absolute_bundle_basis(BIC)
        basis = rel_module.submodule(total_relations).basis()
        
        return ExactStratum(BIC.top,list(basis))
    
        
    def top_level_homology(self, ELG):
        """
        Compute the top level homology of an EmbeddedLevelGraph (ELG).

        Parameters:
            ELG (EmbeddedLevelGraph): The input EmbeddedLevelGraph.

        Returns:
            marked_topl (list): A list of integers representing the marked points on the top level of the BIC.
            num_marked_topl (int): The number of marked points on the top level.
            rel_module (FreeModule): A free module with coefficients in QQ and basis of size num_marked_topl.
            topl_dict (dict): A dictionary mapping marked points on the top level to their indices.
        """
        marked_topl = [l for i,_ in enumerate(ELG.LG.genera)\
                       if ELG.LG.levelofvertex(i)==0 for l in ELG.LG.legsatvertex(i) ]
        num_marked_topl = len(marked_topl)
        rel_module = FreeModule(QQ, num_marked_topl)
        topl_dict = {point : index for index, point in enumerate(marked_topl)}
        return marked_topl, num_marked_topl, rel_module, topl_dict
    
    
    
    
    def absolute_bundle_basis(self, ELG):
        """
        Computes the absolute bundle basis for a given ELG.

        Args:
            ELG (ExactLG): The ExactLG object representing the ELG.

        Returns:
            list: A list of vectors representing the absolute bundle basis. 
            Each vector corresponds to a cycle in the cycle basis of the stgraph of ELG. 
            Each element in the vector corresponds to the value of the cycle at a marked top level vertex. 
            The value is determined by the leg of the corresponding edge in ELG.

        """
        
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
        """
        Computes a basis for the image of the relative cycles (here we mean purely relative)
        under the specialization morphism to the top level.

        Parameters:
            ELG (ExactStratum): The exact stratum to compute the relative bundle basis for.

        Returns:
            list: A list of vectors representing the basis for the image of the relative cycles.

        Raises:
            ValueError: If ELG is not a BIC (Branched Ideal Complex).

        """
        

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
    
    
    

    

    
    def all_ranks(self, ELG):
        """
        Calculates the ranks of the restriction of the bundle of relative homology to each level of the EmbeddedLevelGraph.
        
        Args:
            ELG (EmbeddedLevelGraph or tuple/list): The EmbeddedLevelGraph or the profile of the EmbeddedLevelGraph.
        
        Returns:
            list: A list of integers representing the ranks of the restriction of the bundle of relative homology to each level.
        
        Raises:
            ValueError: If the input ELG is defined for the wrong stratum.
        
        Note:
            If the input ELG is a profile instead of an EmbeddedLevelGraph, we can use any graph with this profile.
            The result is independent of the choice and such a graph can be obtained using self.X.lookup_graph(profile).
        """
        

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
                 
