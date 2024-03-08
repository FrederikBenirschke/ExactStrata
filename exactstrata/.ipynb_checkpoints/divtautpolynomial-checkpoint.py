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
        if len(self.div_taut_poly) ==  0 :
            self._div_taut_poly = ( 0 ,)
            
            
    
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
        assert self[ 0 ]!= 0 
        Z = PowerSeriesRing(SR, 'z', default_prec=dim+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        assert Z(self.div_taut_poly) is not  0 , self.div_taut_poly
        inverted = ( 1 /Z(self.div_taut_poly)).list()
        return DivTautPolynomial(inverted)
    
    def divide(self, other, length):
        Z = PowerSeriesRing(SR, 'z', default_prec=length+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        inverted = ( 1 /Z(self.div_taut_poly)).list()
        return (self * other.invert(length)).trim(length)
        
    
    # Ensure that the polynomial has degree (i.e. codim) at most k
    def trim(self, length):
        # Codimension is greather than length?
        if length >= 0  and length <= len(self.list())- 1 :
            #Then cut off the remainder
            return DivTautPolynomial(self.list()[ 0 :length+ 1 ])
        #Otherwise dont need to change anything
        else:
            return self
        
    
    
    # Makes sure the tuple has the right length
    # Since DivTautPolynomial should be hashable we are returning a new DivTautPolynomial
    def extend_to_length(self, target_len):
        new_poly =  self.tuple()[:target_len] + ( 0 ,)*(target_len - len(self.div_taut_poly))
        return DivTautPolynomial(new_poly)
     
