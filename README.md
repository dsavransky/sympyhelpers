# sympyhelpers

Helper routines and definitions for sympy.  In particular, these methods wrap existing `sympy` functionality to do Newton-Euler and Lagrangian mechanics using Kane-style notation and format the outputs to match (some of) the conventions used in various textbooks including Kasdin & Paley (2011).


## Installation 

`sympyhelpers` is pip installable.  Just run:

```
pip install sympyhelpers
```

## Notational Conventions

The code in this repository adopts the following notational conventions:
* All vectors are represented by lowercase bolded letters ($\mathbf{v}$)
* Unit vectors are represented by lowercase bolded letters with hats ($\mathbf{\hat{v}}$ is the unit vector of $\mathbf v$)
