# ExactStrata







<p align="right">(<a href="#readme-top">back to top</a>)</p>

This ```Sage``` package implements the algorithm developed in the paper [`Tautological classes of strata of exact differentials`](https://arxiv.org/abs/2304.04064). Exact differentials are meromorphic $1$-forms on a compact Riemann surface with zero periods. The goal of (loc. cit.) is to compute the cohomology class of the locus of exact differentials in the moduli space of all differential forms. 

On the  other hand, exact differentials are closely related to meromorphic functions from a Riemann surface to $\mathbb{P}^1$ and the cohomology classes of strata of exact differentials can be used to compute the cohomology classes of Riemann surfaces together with meromorphic functions of fixed ramification, so called admissible cover cycles. For admissible cover cycles the `Sage` package [`admcycles`](https://gitlab.com/modulispaces/admcycles) allows efficient computations and we can compare our algorithm to the computations performed with admcycles.


## Table of contents
- [Features](#features)
- [Installation](#installation)
- [Outline of the algorithm](#outline-of-the-algorithm)
- [Getting Started](#getting-started)
- [License](#license)
- [Contact](#contact)


## Features

- Efficient computation of the cohomology class of the locus of exact differentials.
- Integration with the `admcycles` package for comparative analysis.
- Symbolic computation of divisorial tautological classes.
- Conversion to `ELGTautClass` for pushing forward to the moduli space of curves.


<!-- GETTING STARTED -->
## Installation

To get started with `ExactStrata`, you need to have the following prerequisites:

- `Sage` (version 9.7 or later)
- `admcycles` package (available at [AdmCycles](https://gitlab.com/modulispaces/admcycles)) You can install `admcycles` by running the following command in your `Sage` terminal:

    ```bash 
    sage -pip install git+https://gitlab.com/modulispaces/admcycles
    ```

To install `ExactStrata follows these steps:

1. Clone the repository.

    ```bash
    git clone https://github.com/FrederikBenirschke/ExactStrata.git
    ```

2. Navigate to the root directory.
    ```bash
    cd ExactStrata
    ```

3. Install the package.
    ```bash
    sage -pip install .
    ```

## Outline of the algorithm
The general algorithm presented in the paper is intricate and highly technical, so we will provide a simplified overview of the main idea.

Let $\mathcal{H}(\mu)$ be a stratum of meromorphic differentials and let $\text{Exc}$ be the locus of exact differentials in $\mathcal{H}(\mu)$. Our focus is on the multi-scale compactification $\overline{\mathcal{H}}(\mu)$ , with the objective of determining the class of $\text{Exc}$ in the cohomology ring (Chow ring) of $\overline{\mathcal{H}}(\mu)$.

Exact differentials are described by the vanishing of all periods, which allows us to represent $\text{Exc}$ as the zero locus of a section of a vector bundle $E$ on $\mathcal{H}(\mu)$.
The class of $\text{Exc}$ in $\mathcal{H}(\mu)$ can be computed in terms of the Chern classes of $E$, which are all tautological classes.

The goal is to extend the computation to $\overline{\mathcal{H}}(\mu)$. Both the vector bundle $E$ and the section vanishing on $\text{Exc}$ can be extended to $\overline{\mathcal{H}}(\mu)$.
The extended section vanishes on $\overline{\text{Exc}}\subseteq \overline{\mathcal{H}}(\mu)$, but it vanishes on additional components supported on the boundary. Our approach involves the following steps:
1. **Identify Additional Components**: These components are related to exact differentials supported on boundary components, and their classes can be computed inductively.
2.  **Subtract Additional Contributions:**  After identifying the additional components, we subtract their contributions. The remaining class corresponds to the class of $\overline{\text{Exc}}$.

This straightforward approach works directly when an additional component is a divisor. The main technical challenge arises from the presence of many components of higher codimension. To handle these higher codimension components, we convert them into divisors through a series of blow-ups. This process must be done carefully and in the correct order, as the components can intersect each other. After all blow-ups have been performed we can compute the class of   $\overline{\text{Exc}}$ using a simple Chern class computation. Finally, we have to pushforward the class from the blow-up down to $\overline{\mathcal{H}}(\mu)$. At this stage, we need to know the classes (and all the Chern classes of the normal bundles) of all additional components, which can be determined recursively.




## Getting started

To use `ExactStrata`, you can import the necessary modules in your Sage code:

```python
    from admcycles import * # Provides computations in
    # the tautological ring of the moduli space of curves
    from exactstrata import *
```

Here's an example of how to create an `ExactStratum` object and compute the class of the stratum of exact differentials:

```python
    from admcycles import * 
    from exactstrata import * 

    # Create an ExactStratum for genus one Riemann surfaces
    # and differentials with 3 simple zeros and a pole of order 3
    A= ExactStratum(Stratum([-3,1,1,1]),[])

    # Compute the class of the stratum of exact differentials in the Chow ring
    # The formula is a symbolic formula only
    divtaut =  A.exact_stratum_class()
    print("Divisorial tautological class: ", divtaut)

    # A.to_ELG(divtaut) converts the symbolic class into an ELGTautClass, 
    # a class provided by the package diffstrata
    # to implement tautological classes in strata of multi-scale differentials
    # Finally, A.pushforward(A.to_ELG(divtaut)) takes the ELGTautClass and pushes it forward
    # to the moduli space of curves (where computations are handled by the package admcycles).
    cl= A.pushforward(A.to_ELG(divtaut))
   ```


_For more examples, you can refer to the [Examples notebook](Examples/exactstrata_notebook.ipynb)._ 






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

