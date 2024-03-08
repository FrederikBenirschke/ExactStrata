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
            nb_old, fc_old = self.proper_transform(exc_bd,step - 1 )
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
            return self.blowdown(self.zero_locus_chern()[ 0 ]), self.blowdown(self.zero_locus_chern()[- 1 ])
    
    
    
    
    def blowdown(self, divTaut):
        result = divTaut
        reverse_centers = self.blowup_centers[::- 1 ]
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
        segre = [(-E)**i for i in (ellipsis_range( 1 ,Ellipsis,d-e))]
        
        return - sum(excess[d-e-k]*Xbar*segre[k- 1 ] for k in (ellipsis_range( 1 ,Ellipsis,d-e)))



    def aluffi_formal_tensor(self, A,B, E):
        transform = B.reciprocal().proper_transform(-E).reciprocal()

        result =  (A * transform).list()
        #Now need remove all all terms without E coefficient
        for i, term in enumerate(result):
            result[i] = sum( coeff* E**power for coeff,power in term.coefficients(E)[ 1 :])

       
        return DivTautPolynomial(result)

    
   

       

  
            
