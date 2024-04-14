# ExactStrata







<p align="right">(<a href="#readme-top">back to top</a>)</p>

This ```Sage``` package implements the algorithm developed in the paper [`Tautological classes of strata of exact differentials`](https://arxiv.org/abs/2304.04064).


<!-- GETTING STARTED -->
## Getting Started


The basic usage of exact strata is 
```sh
A= ExactStratum(Stratum([-3,1,1,1]),[])
divtaut =  A.exact_stratum_class()
print("Divisorial tautological class: ", divtaut)
cl= A.pushforward(A.to_ELG(divtaut))
 ```


_For more examples, see the notebook Examples/exactstrata_notebook.ipynb._


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

