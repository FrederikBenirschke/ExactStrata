# ExactStrata







<p align="right">(<a href="#readme-top">back to top</a>)</p>

This ```Sage``` package implements the algorithm developed in the paper [`Tautological classes of strata of exact differentials`](https://arxiv.org/abs/2304.04064). Exact differentials are meromorphic $1$-forms on a compact Riemann surface with zero periods. The goal of (loc. cit.) is to compute the cohomology class of the locus of exact differentials in the moduli space of all differential forms. On the  other hand, exact differentials are closely related to meromorphic functions from a Riemann surface to $\mathbb{P}^1$ and the cohomology classes of strata of exact differentials can be used to compute the cohomology classes of Riemann surfaces together with meromorphic functions of fixed ramification, so called admissible cover cycles. For admissible cover cycles the `Sage` package [`admcycles`](https://gitlab.com/modulispaces/admcycles) allows efficient computations and we can compare our algorithm to the computations performed with admcycles.


<!-- GETTING STARTED -->
## Getting Started


The basic usage of exact strata is 
```sh

# Create an ExactStratum for genus one Riemann surfaces
# and differentials with 3 simple zeros and a pole of order 3
A= ExactStratum(Stratum([-3,1,1,1]),[])

# Compute the class of the stratum of exact differentials in the Chow ring
# The formula is a symbolic formula only
divtaut =  A.exact_stratum_class()
print("Divisorial tautological class: ", divtaut)

# A.to_ELG(divtaut) converts the symbolic class into an ELGTautClass, a class provided by the package diffstrata
# to implement tautological classes in strata of multi-scale differentials
# Finally, A.pushforward(A.to_ELG(divtaut)) takes the ELGTautClass and pushes it forward
# to the moduli space of curves (where computations are handled by the package admcycles).
cl= A.pushforward(A.to_ELG(divtaut))
 ```


_For more examples, see the  notebook [Examples/exactstrata_notebook.ipynb](Examples/exactstrata_notebook.ipynb) ._ 


### Prerequisites



Requires the package ```admcycles```. 
[`AdmCycles`](https://gitlab.com/modulispaces/admcycles)






<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the MIT License. 
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Frederik Benirschke - benirschke.math@gmail.com
Project Link: [https://github.com/FrederikBenirschke/ExactStrata](https://github.com/FrederikBenirschke/ExactStrata)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

