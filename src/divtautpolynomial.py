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
        """
        Initializes the DivTautPolynomial object with the given div_taut_poly.
        Checks if the input div_taut_poly is a tuple, converts it to a tuple if needed, and handles empty tuples.
        """
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
        """
        Returns the value of the `_div_taut_poly` attribute.

        :return: The value of the `_div_taut_poly` attribute.
        :rtype: tuple
        """
        return self._div_taut_poly
    
    def tuple(self):
        """
        Returns the div_taut_poly attribute of the DivTautPolynomial object.
        """
        return self.div_taut_poly
    
    def list(self):
        """
        Returns a list representation of the `div_taut_poly` attribute of the `DivTautPolynomial` object.

        :return: A list representation of the `div_taut_poly` attribute.
        :rtype: list
        """
        return list(self.tuple())
    
    def __eq__(self, other):
        """
        Check if the current object is equal to another object.

        Parameters:
            other (object): The object to compare with.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if other == None:
            return False
        return self.div_taut_poly == other.div_taut_poly
    
    def __hash__(self):
        """
        Return the hash value of the `_div_taut_poly` attribute.
        """
        return hash(self._div_taut_poly)
        
       
        
    def __repr__(self):
        """
        Returns a string representation of the `div_taut_poly` attribute of the `DivTautPolynomial` object.

        :return: A string representation of the `div_taut_poly` attribute.
        :rtype: str
        """
        return str(self.div_taut_poly)
    
    def __getitem__(self, key):
        """
        Returns the value associated with the given key in the `div_taut_poly` dictionary.

        Parameters:
            key (Any): The key to access the value in the dictionary.

        Returns:
            Any: The value associated with the given key.
        """
        return self.div_taut_poly[key]
    

    def __add__(self, other):
        """
        Adds two `DivTautPolynomial` objects together.

        Parameters:
            other (DivTautPolynomial): The other `DivTautPolynomial` object to add.

        Returns:
            DivTautPolynomial: The sum of the two `DivTautPolynomial` objects.

        Raises:
            TypeError: If `other` is not an instance of `DivTautPolynomial`.
        """
       
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) + S(other.div_taut_poly)).list())
        
    
    def __mul__(self, other):
        """
        Multiplies two `DivTautPolynomial` objects together.

        Parameters:
            other (DivTautPolynomial or int): The other `DivTautPolynomial` object to multiply with, or an integer to scale the polynomial.

        Returns:
            DivTautPolynomial: The product of the two `DivTautPolynomial` objects, or the scaled `DivTautPolynomial` if `other` is an integer.

        Raises:
            ValueError: If `other` is not an instance of `DivTautPolynomial` or an integer.
        """
       
        if isinstance(other, sage.rings.integer.Integer):
            return DivTautPolynomial([coeff*other for coeff in self.div_taut_poly ])
        if not isinstance(other, DivTautPolynomial):
            
            raise ValueError(str(other)+'is of type' +str(type(other))+' and not a DivTautPolynomial')
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) * S(other.div_taut_poly)).list())
    
    def __rmul__(self, other):
        """
        Right multiplication method to multiply a `DivTautPolynomial` object by an integer.
        
        Parameters:
            other (sage.rings.integer.Integer): The integer to multiply the `DivTautPolynomial` object by.
        
        Returns:
            DivTautPolynomial: The product of the `DivTautPolynomial` object and the integer.
        """
        if isinstance(other, sage.rings.integer.Integer):
            return DivTautPolynomial([coeff*other for coeff in self.div_taut_poly ])
        
    def __sub__(self, other):
        """
        Returns the difference of two DivTautPolynomials.

        Parameters:
            other (DivTautPolynomial): The other DivTautPolynomial to subtract.

        Returns:
            DivTautPolynomial: The difference of the two DivTautPolynomials.
        """
        
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial((S(self.div_taut_poly) - S(other.div_taut_poly)).list())
        
        
    
    def __pow__(self, exponent):
        """
        Raises the DivTautPolynomial to the power of the given exponent.

        Parameters:
            exponent (int): The exponent to raise the DivTautPolynomial to.

        Returns:
            DivTautPolynomial: The result of raising the DivTautPolynomial to the power of the exponent.
        """
       
        S = SR['t']; (t,) = S._first_ngens(1)
        return DivTautPolynomial(((S(self.div_taut_poly))**exponent).list())
        
        
    @property
    def codim(self):
        """
        Returns the codimension of the DivTautPolynomial. 

        :return: An integer representing the codimension of the DivTautPolynomial.
        :rtype: int
        """
        
        return len(self.list())
    
    
    def reciprocal(self):
        """
        Returns a new `DivTautPolynomial` object with the coefficients of the original polynomial reversed.

        Returns:
            DivTautPolynomial: A new `DivTautPolynomial` object with the coefficients of the original polynomial reversed.
        """
        return DivTautPolynomial(tuple(reversed(self.div_taut_poly)))
    
    
    def substitute(self, old_variables, new_variables, verbose = False):
        """
        Replaces the variables of DivTautPolynomial by a different set of variables.

        Args:
            old_variables (list): A list of variables to be replaced.
            new_variables (list): A list of variables to replace the old variables with.
            verbose (bool, optional): If True, prints the substitution dictionary. Defaults to False.

        Raises:
            ValueError: If the number of old variables is not equal to the number of new variables.

        Returns:
            DivTautPolynomial: A new DivTautPolynomial object with the variables replaced.
        """
       
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
        """
        Performs a proper transformation on the DivTautPolynomial object by substituting the variable 't' with 't + var'.

        Parameters:
            var (int or float): The value to add to the variable 't'.

        Returns:
            DivTautPolynomial: A new DivTautPolynomial object with the proper transformation applied.
        """
       
        S = SR['t']; (t,) = S._first_ngens(1)
        poly =  S(self.div_taut_poly)
        return DivTautPolynomial(poly(t= S(t+var)).list())
    
    def tensor(self, var):
        """
        Calculates the tensor product of the current `DivTautPolynomial` object with a given variable.

        Parameters:
            var (int or float): The variable to perform the tensor product with.

        Returns:
            DivTautPolynomial: The tensor product of the current object with the given variable.

        Note:
            The tensor product is calculated by first taking the reciprocal of the current object,
            then performing a proper transformation with the given variable, and finally taking the
            reciprocal again.
        """
        return self.reciprocal().proper_transform(var).reciprocal()
    
   
    
    
    
    def invert(self, dim):
        """
        Inverts the DivTautPolynomial in the Chow ring CH^*(X).

        Args:
            dim (int): The dimension of the ambient space, which is the maximal length of the inverse minus 1.

        Returns:
            DivTautPolynomial: The inverted DivTautPolynomial.

        Raises:
            AssertionError: If the leading term of the current DivTautPolynomial is zero.

        Notes:
            This function uses the leading term of the DivTautPolynomial to check if it is non-zero.
            It then creates a PowerSeriesRing with the given dimension and a single variable 'z'.
            The current DivTautPolynomial is then converted to a polynomial in this ring.
            The inverse of this polynomial is computed and returned as a new DivTautPolynomial.
        """
       
        # Leading term must be non-zero
        assert self[ 0 ]!= 0 
        Z = PowerSeriesRing(SR, 'z', default_prec=dim+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        assert Z(self.div_taut_poly) is not  0 , self.div_taut_poly
        inverted = ( 1 /Z(self.div_taut_poly)).list()
        return DivTautPolynomial(inverted)
    
    def divide(self, other, length):
        """
        Divides the current DivTautPolynomial by another DivTautPolynomial and returns the result as a new DivTautPolynomial.
        
        Args:
            other (DivTautPolynomial): The DivTautPolynomial to divide by.
            length (int): The expected length of the result.
        
        Returns:
            DivTautPolynomial: The result of the division, trimmed to the specified length.
        
        Raises:
            AssertionError: If the leading term of the current DivTautPolynomial is zero.
        
        """
       
        Z = PowerSeriesRing(SR, 'z', default_prec=length+ 1 , names=('z',)); (z,) = Z._first_ngens(1)
        inverted = ( 1 / Z(self.div_taut_poly)).list()
        return (self * other.invert(length)).trim(length)
        
    
    
    def trim(self, length):
        """
        Ensures that the polynomial has degree (i.e. codim) at most k.
        
        Args:
            length (int): The maximum degree allowed for the polynomial.
        
        Returns:
            DivTautPolynomial: The trimmed polynomial based on the specified length.
        
        """
       
        # Codimension is greather than length?
        if length >= 0  and length <= len(self.list()) - 1 :
            #Then cut off the remainder
            return DivTautPolynomial(self.list()[ 0 :length + 1 ])
        #Otherwise dont need to change anything
        else:
            return self
        
    
    
    
    def extend_to_length(self, target_len):
        """
        Extends the DivTautPolynomial to a specified target length by padding with zeros.
        
        Args:
            target_len (int): The desired length of the DivTautPolynomial.
        
        Returns:
            DivTautPolynomial: The extended DivTautPolynomial to the specified length.
        """
        
        new_poly =  self.tuple()[:target_len] + ( 0 ,)*(target_len - len(self.div_taut_poly))
        return DivTautPolynomial(new_poly)
     
