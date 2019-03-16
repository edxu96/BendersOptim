# Test Benders Algorithm for Mixed Integer Linear Programming
# Edward J. Xu
# March 16th, 2019
# ---------------------------------------------------------------------------
push!(LOAD_PATH, "$(homedir())/Desktop/BendersOptim_Julia_EDXU")
cd("$(homedir())/Desktop/BendersOptim_Julia_EDXU")
using BendersMILP_EDXU
using JuMP
using GLPKMathProgInterface
# ---------------------------------------------------------------------------
n_x = 1
n_y = 1
vec_max_y = [10]
vec_c = [5]
vec_f = [-3]
vec_b = [4; 0; -13]
mat_a = [1; 2; 1]
mat_b = rand(3,1)
mat_b[1] = 2
mat_b[2] = -1
mat_b[3] = -3
# ---------------------------------------------------------------------------
BendersMILP(n_x = n_x,
            n_y = n_y,
            vec_max_y = vec_max_y,
            vec_c = vec_c,
            vec_f = vec_f,
            vec_b = vec_b,
            mat_a = mat_a,
            mat_b = mat_b)
