# BendersOptim

Benders algorithm to solve mixed integer linear programming, especially stochastic programming in seconds!

This repo has been archived in 190728. New update will be made to edxu96/MatrixOptim, which is the aggregation of robust optimization and matrix optimization. It's in matrix form as well, and there is a new function tool to convert the model to matrix form.

***

According to wikipedia:
> Benders decomposition (or Benders' decomposition) is a technique in mathematical programming that allows the solution
> of very large linear programming problems that have a special block structure. This block structure often occurs in
> applications such as stochastic programming as the uncertainty is usually represented with scenarios. The technique is
> named after Jacques F. Benders.

https://en.wikipedia.org/wiki/Benders_decomposition

There are two algorithms, with one for standard MILP, and the other one specifically for stochastic programming without
integer variables in second stage.  

For detailed explanation, refer
to [Cookbook for Benders Decomposition, EDXU](files/cookbook.pdf).  

# How to use

## 1. Import

Put the file `Benders.jl`, `Benders_milp.jl` and `Benders_lshaped.jl` in your working folder. Write the following code
in your file:  
```Julia
using Benders
```

## 2. Generic Benders Decomposition for Mixed Integer Linear Programming

![Standard MILP](images/standard_milp.png)

```Julia
Bender.milp(n_x, n_y, vec_min_y, vec_max_y, vec_c, vec_f,
  vec_b, mat_a, mat_b, epsilon, timesIterationMax)
```

## 3. L-Shaped Benders Decomposition for Stochastic Programming without Integer Variables in Second Stage

![Stochastic Programming without Integer Variables in Second Stage](images/l-shaped.png)

```Julia
Bender.lshaped(n_x, n_y, vec_min_y, vec_max_y, vec_f,
  vec_pi, mat_c, mat_h, mat3_t, mat3_w, epsilon, timesIterationMax)
```

***

Edward J. Xu (edxu96@outlook.com) (edxu96.github.io)  
Version: 2.1  
Date: July 28th, 2019
