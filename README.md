# BendersOptim
Benders algorithm to solve mixed integer linear programming, especially stochastic programming in seconds!

Benders decomposition (or Benders' decomposition) is a technique in mathematical programming that allows the solution of very large linear programming problems that have a special block structure. This block structure often occurs in applications such as stochastic programming as the uncertainty is usually represented with scenarios. The technique is named after Jacques F. Benders.

https://en.wikipedia.org/wiki/Benders_decomposition

# How to use

## 1. Import 

Put the file `Benders.jl`, `Benders_milp.jl` and `Benders_lshaped.jl` in your working folder. Write the following code in your file:  
```Julia
using Benders
```

## 2. Generic Benders Decomposition for Mixed Integer Linear Programming

![Standard MILP]（../images/standard_milp.png）

```Julia
Bender.milp(n_x, n_y, vec_min_y, vec_max_y, vec_c, vec_f, 
  vec_b, mat_a, mat_b, epsilon, timesIterationMax)
```

## 3. L-Shaped Benders Decomposition for Stochastic Programming without Integer Variables in Second Stage

![Stochastic Programming without Integer Variables in Second Stage]（../images/l-shaped.png）

```Julia
Bender.lshaped(n_x, n_y, vec_min_y, vec_max_y, vec_f, 
  vec_pi, mat_c, mat_h, mat3_t, mat3_w, epsilon, timesIterationMax)
```

***

Edward J. Xu  
edxu96@outlook.com  
Version: 2.0  
Date: April 5th, 2019

