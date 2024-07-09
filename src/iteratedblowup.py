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


     

class IteratedBlowup(SageObject):
    
    def __init__(self, exact_stratum, global_profile, level, base, blowup_centers):
        """
        Initializes the IteratedBlowup object with the given parameters.
        
        Parameters:
            exact_stratum (ExactStratum): The exact stratum object.
            global_profile (GlobalProfile): The global profile object.
            level (int): The level of the blow-up.
            base (str): The base.
            blowup_centers (list): List of ExactBoundaryStrata objects to be blown up.
        
        Returns:
            None
        """
        
        
      
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
        """
        Returns a string representing the blowup of ExactStratum, the global profile, the base, the blowup centers, and the divisors.
        """
        return "Blowup of ExactStratum " + self.exact_stratum.__repr__()         + "\n"+ "Global Profile: " + self.global_profile.__repr__()        + "\n"+ "Base:" + str(self.base.profile)+','+str( self.base.levels)        + "\n"+ "Blowup centers: " + str([(bd.profile, bd.levels) for  bd in self.blowup_centers])        + '\n'+ 'Divisors: ' + str(list(self.divisors.keys()))
    
    def codim(self, deg, undeg):
        """
        A function that calculates the codimension of a differential object.

        Parameters:
            self: The IteratedBlowup object.
            deg: The degree of the differential object.
            undeg: The undegenerate part of the differential object.

        Returns:
            The codimension of the differential object.
        """
        assert deg.is_contained(undeg)
        return deg.total_codim - undeg.total_codim
   
    # Creates the iterated blow up
    def sub_blowup(self, exc_bd):
        """
        A function that performs a sub blowup operation.
        
        Parameters:
            self: The IteratedBlowup object.
            exc_bd: The exact boundary stratum object.
        
        Returns:
            A tuple containing the new IteratedBlowup object and the step adjustment dictionary.
        """
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
        """
        Check if all boundary strata in the `blowup_centers` list have the expected codimension.

        Returns:
            bool: True if all boundary strata have the expected codimension, False otherwise.
        """
        for bd in self.blowup_centers:
            # All boundary strata have expected codim?
            if bd.total_codim == self.basic_exact_stratum.total_codim:
                return True
    
   
  
    
    
    
  
    def aluffi_decomposition(self, deg, undeg, step):
        """
        Decomposes the blowup of Y=undeg inside X=deg as X -> I -> Y, where I has the same underlying profile as X but potentially less levels.
        The normal bundle N(I/Y) is the normal bundle of a boundary stratum in the stratum, which intersects the base transversely.
        The proper transform is simply a pullback and the normal bundle N(X/I) is the product of normal bundles of proper transforms.
        
        Parameters:
            deg (ExactStratum): The exact stratum object representing X.
            undeg (ExactStratum): The exact stratum object representing Y.
            step (int): The step number.
        
        Returns:
            tuple: A tuple containing the product of normal bundles and the last normal bundle multiplied by the proper transform.
        """
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
        """
        Caches the result of the `building_transform` method.
        
        Calculates the building transform for a given step.
        For the i-th building in the blowup, compute the proper transform
        on the last step before the i-th building is blown up
        First output is the total Chern class, second is the fundamental class with multiplicity
        
        Parameters:
            step (int): The step number.
        
        Returns:
            tuple: A tuple containing the product of normal bundles and the last normal bundle multiplied by the proper transform.
        """
        
        bd = self.blowup_centers[step]
        return self.proper_transform(bd, step- 1 )
    
    
   
    @cached_method
    def building_segre(self, step):
        """
        Caches the result of the `building_segre` method.
        
        Computes the Segre class for a given step in the blowup.
        The dimension of the blow-up center gives an upper bound for the Length of the total Segre class.
        Returns the Segre class by multiplying each term with the fundamental class and inverting the total Chern class.
        """
        #The dimension of the blow-up center, gives upper bound for the Length of the total Segre class
        dim = self.exact_stratum.dim - self.blowup_centers[step].total_codim
        total_chern, fc = self.building_transform(step)
       
        return [term* fc for term in total_chern.invert(dim)]
    
   
    
    
    
    def final_blowup_difference(self, exc_bd):
        """
        Calculate the final blowup difference between the given exceptional boundary stratum (exc_bd) and the last step of the blowup.

        Parameters:
            exc_bd (ExactBoundaryStratum): The exceptional boundary stratum to calculate the final blowup difference for.

        Returns:
            tuple: A tuple containing the product of normal bundles and the last normal bundle multiplied by the proper transform.

       Note: Only implemented for the case that exc_bd is not contained in ANY blowup-center
        """
        return self.proper_transform(exc_bd,len(step)- 1 )
        
    
    
    
    @cached_method
    def proper_transform(self, exc_bd, step):
        """
        Returns the total Chern class and fundamental class of proper transform

        Parameters:
            self: The current instance of the class
            exc_bd: The exceptional boundary stratum
            step: The step in the transformation process
        
        Returns:
            Tuple containing the total Chern class and fundamental class of the proper transform
        """
        if step == - 1 :
            nb = exc_bd.formal_normal_bundle(self.base, mult = True)
            return nb, nb[- 1 ] 
        else:
            nb_old, fc_old = self.proper_transform(exc_bd,step - 1 )
            nb_diff, fc_diff = self.proper_transform_diff(exc_bd, step)
            return  nb_old+nb_diff, fc_old+fc_diff 
        
        
        
    
    
    @cached_method
    def proper_transform_diff(self, exc_bd, step):
        """
        Returns the difference between the total Chern class of the proper transform and the pullback
        on the i-th step of the blow up
        
        Parameters:
            self: The current instance of the class
            exc_bd: The exceptional boundary stratum
            step: The step in the transformation process
        
        Returns:
            Tuple containing the total Chern class and fundamental class of the proper transform
        """
        
        if exc_bd in self.blowup_centers:
            assert step < self.blowup_centers.index(exc_bd)
        
        
        # No difference in the first step
        if step == - 1 :
            return DivTautPolynomial([ 0 ]), 0 
        
        
        assert step >= 0  and step< self.num_of_blowups, str(step) + 'is no in range ' + str(self.num_of_blowups)
        # Locate blowup_center
        center = self.blowup_centers[step]
        base = self.base
        
        # exc_bd is disjoint from blowup center?
        if center.are_disjoint(exc_bd):
            # No difference
            return DivTautPolynomial([ 0 ]), 0 
        else:
            intersection = exc_bd + center
            
             # Codimensions of X,Y,W in Z, in that order
            codims = [self.codim(intersection, base), self.codim(exc_bd, base), self.codim(center, base)]
            codimX, codimY, codimW = codims
            # Check if the intersection is transversal:
            if codimX == codimY+codimW:
                # No difference
                return DivTautPolynomial([ 0 ]), 0 
            
            #We are now setting up aluffis formula
            # Y - variety whose proper transform we are trying to compute
            # W - center of blowup
            # X - intersection of variety and blowup center
            # Z - ambient space
            cXW, fXW = self.aluffi_decomposition(intersection, center, step- 1 )
            cXY, _ = self.aluffi_decomposition(intersection, exc_bd, step- 1 )
            cWZ, _ = self.aluffi_decomposition(center, base, step - 1 )
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
        """
        Computes the class of the given exact stratum without blowing up.

        This function computes the class of the basic exact stratum without blowing up by using the formal xi at the given level.
        It iterates over the building divisors and subtracts the formal normal bundle of each bic index. 
        The rank is obtained from the global profile ranks at the given level. The total Chern of the proper transform of A_P^{[0]} is computed using the div_class and rank. The difference is then computed by summing the formal normal bundle of each blowup center. The function returns the fundamental class minus the difference.

        Returns:
           SymbolicExpression: The fundamental class minus the difference.
        """
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

        return fc - diff
        
        
        
        
    
    @cached_method
    def zero_locus_chern(self):
        """
        Calculates the Chern class and fundamental class of the zero locus of the given exact stratum.

        This method calculates the Chern class and fundamental class of the zero locus of the given exact stratum. 
        It does this by computing the formal xi at the given level, subtracting the formal normal bundle of each building divisor,
         and subtracting the formal normal bundle of each blowup center. 
         The rank is obtained from the global profile ranks at the given level. 
         The total Chern of the proper transform of A_P^{[0]} is computed using the div_class and rank. 
         The difference is then computed by summing the formal normal bundle of each blowup center. 
         The function returns the Chern class and fundamental class of the zero locus.

        Returns:
            Tuple[DivTautPolynomial, SymbolicExpression]: A tuple containing the Chern class and fundamental class of the zero locus.
        """
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
        fc = chern[- 1 ]
       
        # In order to obtain the total Chern class of the Pullback instead
        # we neeed to compute the difference using Aluffis formula
        exc_bd = self.basic_exact_stratum
        for step in range(self.num_of_blowups):
            difference_ch, difference_fc =  self.proper_transform_diff(exc_bd,step)
            chern -= difference_ch
            
            fc -= difference_fc
                                       
        
        
        
        
        return chern, fc
    
    def zero_locus_blowdown(self):
        """
        Returns the blowdown of the zero locus.

        This function calculates the blowdown of the zero locus by calling the `blowdown` function
        twice, once with the Chern class of the zero locus and once with the fundamental class.

        Returns:
            A tuple containing the blowdown of the Chern class and the blowdown of the fundamental class.
        """
            
            return self.blowdown(self.zero_locus_chern()[ 0 ]), self.blowdown(self.zero_locus_chern()[- 1 ])
    
    
    
    
    def blowdown(self, divTaut):
        """
        Apply the blowdown operation to the given divisorial tautological polynomial.

        Parameters:
            divTaut (DivTautPolynomial): The divisorial tautological polynomial to apply the blowdown operation to.

        Returns:
            DivTautPolynomial: The result of applying the blowdown operation to the given divisorial tautological polynomial.
        """
        result = divTaut
        reverse_centers = self.blowup_centers[::- 1 ]
        for center in reverse_centers:
            step = self.blowup_centers.index(center)
            result = self._blowdown_step(step,result)
        return result
        
        
        

       def _blowdown_step(self, step, divTaut):
        """
        Applies the blowdown operation to the given divisorial tautological polynomial on the i-th blowup step.

        Parameters:
            step (int): The index of the blowup step.
            divTaut (DivTautPolynomial): The divisorial tautological polynomial to apply the blowdown operation to.

        Returns:
            DivTautPolynomial: The result of applying the blowdown operation to the given divisorial tautological polynomial.

        Raises:
            None

        Algorithm:
            1. If the input is a DivTautPolynomial, recursively apply the blowdown operation to each term.
            2. Retrieve the center of the blowup step and the corresponding E-scheme.
            3. Compute the codimension within the base and the dimension of the exact stratum.
            4. Define a helper function `segre_k` to compute the Segre class.
            5. Use formula (5.2) from the paper to compute the pushforward of the divisorial tautological polynomial.
            6. Set the constant term of the pushforward to the constant term of the input.
            7. Return the sum of the pushforward coefficients.

        Note:
            The blowdown operation is applied iteratively to the input divisorial tautological polynomial, with each step representing a blowup step.
        """
        
        if type(divTaut) is DivTautPolynomial:
            return DivTautPolynomial([self._blowdown_step(step, term) for term in divTaut.list()])
        
        
        center = self.blowup_centers[step]
        E =  self.E[center]
        
        #Codim within base
        codim = center.total_codim - self.base.total_codim
        dim = self.exact_stratum.dim - center.total_codim
        segre = []
        
        def segre_k(k):
            if k< 0  or k > dim:
                return  0 
            else:
                
                #Only compute the Segre class if necessary
                segre = self.building_segre(step)
                return segre[k]
        
 

        # Use formula (5.2) from the paper for the pushforward

        coeff = [ (- 1 )**(k- 1 )* aj* segre_k(k - codim)                  for aj, k in divTaut.coefficients(E)]

        # The constant term stays the same
        coeff[ 0 ] = divTaut.coefficients(E)[ 0 ][ 0 ]
        return sum(coeff)
    
 
    ######
    # Aluffi logic
    ########
    
    
    def aluffi_formula(self, cherns, codims, E, fXW):
        """
        Calculate the total chern class and fundamental class for the given cherns, codims, E, and fXW.

        Parameters:
            cherns: Tuple containing cXY, cXW, and cWZ
            codims: List containing the codimensions of X, Y, and W in Z
            E: Value associated with the blowup center
            fXW: Value representing the proper transform

        Returns:
            Tuple containing the total chern class and fundamental class
        """
        cXY, cXW, cWZ = cherns
        total_chern = self.aluffi_normal_bundle(cherns, codims, E)
        fundamental_class = self.aluffi_proper_transform(cherns, codims, E, fXW)
       
        return total_chern, fundamental_class
    
    def aluffi_normal_bundle(self, cherns, codims, E):
        """
        Calculate the normal bundle of a blowup given the Chern classes, codimensions, and a specific value E.

        Parameters:
            cherns: Tuple containing the Chern classes cXY, cXW, and cWZ
            codims: List containing the codimensions of X, Y, and W in Z
            E: Value associated with the blowup center

        Returns:
            The normal bundle of the blowup trimmed to the codimension of Y.
        """
        codimX, codimY, codimW = codims
        cXY, cXW, cWZ = cherns
        d = codimY 
        e = codimX - codimW
        excess = cWZ.divide(cXY, d-e)
       
        return self.aluffi_formal_tensor(cXW, excess, E).trim(codimY)
    
    def aluffi_proper_transform(self, cherns, codims, E, fXW):
        """
        Calculate the proper transform of a blowup using Aluffi's formula.

        Parameters:
            cherns (tuple): A tuple containing the Chern classes cXY, cXW, and cWZ.
            codims (list): A list containing the codimensions of X, Y, and W in Z.
            E (int): The value associated with the blowup center.
            fXW (object): The fundamental class of Xbar.

        Returns:
            object: The proper transform of the blowup.
        """
        codimX, codimY, codimW = codims
        cXY, cXW, cWZ = cherns
        
        # To obtain the fundamental class of Xbar we need to consider the correct multiplicity
        Xbar = fXW
        d = codimY 
        e = codimX - codimW
        
        excess = cWZ.divide(cXY, d-e)
        segre = [(-E)**i for i in (ellipsis_range( 1 ,Ellipsis,d-e))]
        
        return - sum(excess[d-e-k]*Xbar*segre[k- 1 ] for k in (ellipsis_range( 1 ,Ellipsis,d-e)))



    def aluffi_formal_tensor(self, A,B, E):
        """
        Calculate the Aluffi formal tensor of two DivTautPolynomials.

        Parameters:
            A (DivTautPolynomial): The first DivTautPolynomial.
            B (DivTautPolynomial): The second DivTautPolynomial.
            E (int): The value associated with the blowup center.

        Returns:
            DivTautPolynomial: The Aluffi formal tensor of A and B.

        This function calculates the Aluffi formal tensor of two DivTautPolynomials A and B.
        It first calculates the proper transform of B using the reciprocal of E.
        Then it multiplies A with the proper transform and converts the result to a list.
        Finally, it removes all terms without the E coefficient and returns the resulting DivTautPolynomial.
        """
        transform = B.reciprocal().proper_transform(-E).reciprocal()

        result =  (A * transform).list()
        #Now need remove all all terms without E coefficient
        for i, term in enumerate(result):
            result[i] = sum( coeff* E**power for coeff,power in term.coefficients(E)[ 1 :])

       
        return DivTautPolynomial(result)

    
   

       

  
            
