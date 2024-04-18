import sage

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

from sage.modules.free_module import FreeModule  # pylint: disable=import-error
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.rings.power_series_ring import PowerSeriesRing
# from exactstrata.profile import *
# from exactstrata.exactboundarystratum import *
# from exactstrata.iteratedblowup import *
# from exactstrata.divtautpolynomial import *
from sage.calculus.var import var






class DivTautPolynomial(SageObject):
    ''' A DivTautPolynomial is a tuple where each entry is a SymbolicExpression, i.e. a polynomial  in the variables xi and D_i, where i is the index of a BIC.
    A DivTautPolynomial formally represents a tautological class in the Chow ring of the moduli space of multi-scale differentials.
    DivTautPolynomial needs to be cachable and its thue implemented using tuples.'''
    
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
        if len(self.div_taut_poly) ==  0 :
            self._div_taut_poly = (0,)
            
            
    
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
        ''' Returns the sum of two DivTautPolynomials (as DivTautPolynomial).'''
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) + S(other.div_taut_poly)).list())
        
    
    def __mul__(self, other):
        ''' Returns the product of two DivTautPolynomials (or of a DivTautPolynomial with an integer). The result is a DivTautPolynomial.'''
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
        ''' Returns the difference of two DivTautPolynomials. '''
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) - S(other.div_taut_poly)).list())
        
        
    
    def __pow__(self, exponent):
        ''' Raises a DivTautPolynomial to the power exponent.'''
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial(((S(self.div_taut_poly))**exponent).list())
        
        
    @property
    def codim(self):
        ''' The degree of the DivTautPolynomial. In the case that DivTautPolynomials represent the normal bundle of a subvariety
        this agree with the codimension. '''
        return len(self.list())
    
    
    def reciprocal(self):
        ''' Reverses the list of coefficients of the polynomial.'''
        return DivTautPolynomial(tuple(reversed(self.div_taut_poly)))
    
    
    def substitute(self, old_variables, new_variables, verbose = False):
        ''' Replaces the variables of DivTautPolynomial by a different set of variables.'''
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
        ''' Implements Fulton's formula for the class of a proper transform under a blow-up.'''
        S = SR['t']; (t,) = S._first_ngens(1)
        poly =  S(self.div_taut_poly)
        return DivTautPolynomial(poly(t= S(t+var)).list())
    
    def tensor(self, var):
        return self.reciprocal().proper_transform(var).reciprocal()
    
   
    
    
    
    def invert(self, dim):
        ''' Inverts the DivTautPolynomial in the Chow ring CH^*(X).
        dim - dimension of ambient space = maximal length of the inverse - 1.'''
        # Leading term must be non-zero
        assert self[ 0 ]!= 0 
        Z = PowerSeriesRing(SR, 'z', default_prec=dim+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        assert Z(self.div_taut_poly) is not  0 , self.div_taut_poly
        inverted = ( 1 /Z(self.div_taut_poly)).list()
        return DivTautPolynomial(inverted)
    
    def divide(self, other, length):
        ''' Retursn the division of DivTautPolynomial by another one.
        ''' Requires the expected length of the result as additional parameter.'''
        Z = PowerSeriesRing(SR, 'z', default_prec=length+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        inverted = ( 1 / Z(self.div_taut_poly)).list()
        return (self * other.invert(length)).trim(length)
        
    
    
    def trim(self, length):
        ''' Ensures that the polynomial has degree (i.e. codim) at most k.'''
        # Codimension is greather than length?
        if length >= 0  and length <= len(self.list()) - 1 :
            #Then cut off the remainder
            return DivTautPolynomial(self.list()[ 0 :length + 1 ])
        #Otherwise dont need to change anything
        else:
            return self
        
    
    
    
    def extend_to_length(self, target_len):
        ''' Returns a new DivTautPolynomial by adding zeros until the tuple encoding the DivTautPolynomial has the desired length.
        Since DivTautPolynomial needs to be hashable, the result is a new DivTautPolynomial. '''
        new_poly =  self.tuple()[:target_len] + ( 0 ,)*(target_len - len(self.div_taut_poly))
        return DivTautPolynomial(new_poly)
     
