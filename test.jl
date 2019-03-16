# Test Benders Algorithm for Mixed Integer Linear Programming
# Edward J. Xu
# March 16th, 2019
# ---------------------------------------------------------------------------
push!(LOAD_PATH, "$(homedir())/Desktop/BendersOptim_Julia_EDXU")
cd("$(homedir())/Desktop/BendersOptim_Julia_EDXU")
using BendersMilp_EDXU
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
BendersMilp(n_x = n_x,
            n_y = n_y,
            vec_max_y = vec_max_y,
            vec_c = vec_c,
            vec_f = vec_f,
            vec_b = vec_b,
            mat_a = mat_a,
            mat_b = mat_b,
            epsilon = 0,
            timesIterationMax = 5)
# ---------------------------------------------------------------------------
n_x_2 = 2
n_y_2 = 2
vec_max_y_2 = [10; 10]
vec_c_2 = [5; 3]
vec_f_2 = [-3; 1]
vec_b_2 = [4; 0; -13]
mat_a_2 = [1 3; 2 1; 1 -5]
mat_b_2 = [2 -4; -1 2; -3 1]
# ---------------------------------------------------------------------------
BendersMilp(n_x = n_x_2,
            n_y = n_y_2,
            vec_max_y = vec_max_y_2,
            vec_c = vec_c_2,
            vec_f = vec_f_2,
            vec_b = vec_b_2,
            mat_a = mat_a_2,
            mat_b = mat_b_2,
            epsilon = 0,
            timesIterationMax = 10)
