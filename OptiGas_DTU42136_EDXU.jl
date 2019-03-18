# OptiGas: Optimize Gas Network using Benders Algorithm
# Edward J. Xu, edxu96@outlook.com
# March 18th, 2019
push!(LOAD_PATH, "$(homedir())/Desktop/OptiGas, DTU42136")
cd("$(homedir())/Desktop/OptiGas, DTU42136")
using BendersMilp_EDXU
#-------------------------------------------------------------------
# 1. Parameters
vec_nameNodes = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "j"]
numNodes = length(vec_nameNodes)
vec_xNode = [1    3    8    2    4    5    6    7    9    9]
vec_yNode = [1    2    1    7    9    5    3    7    9    4]
mat_arcTwoNodes = [0    1    0    0    0    0    0    0    0    0;
                   0    0    0    0    0    1    0    0    0    0;
                   0    0    0    0    0    0    0    0    0    1;
                   1    0    0    0    1    0    0    0    0    0;
                   0    0    0    0    0    1    0    0    0    0;
                   1    0    0    1    0    0    0    1    0    0;
                   0    1    1    0    0    0    0    0    0    0;
                   0    0    0    0    0    0    0    0    1    0;
                   0    0    0    0    1    0    0    0    0    0;
                   0    0    0    0    0    1    1    0    1    0]
vec_netEject = [0    0   0    0    0   90    0    0   0   80 ]
vec_netInject = [10   30  20   10   40    0   25   20  15    0 ]
mat_distance = zeros(Float64, numNodes, numNodes)
for m = 1: numNodes
    for n = 1: numNodes
        mat_distance[m, n] = floor(sqrt( (vec_xNode[m] - vec_xNode[n]) * (vec_xNode[m] - vec_xNode[n]) +
                                 (vec_yNode[m] - vec_yNode[n]) * (vec_yNode[m] - vec_yNode[n]) ) )
    end
end
mat_fixedCost = zeros(Float64, numNodes, numNodes)
for m = 1: numNodes
    for n = 1: numNodes
        mat_fixedCost[m, n] = 10 * mat_distance[m, n]#*(1-mat_arcTwoNodes[m,n])
    end
end
# Fomulation of Matrix A
mat_a_1 = zeros(numNodes^2, numNodes^2)
for mn = 1: numNodes^2
    mat_a_1[mn, mn] = -1
end
mat_a_2 = zeros(numNodes^2, numNodes^2)
mat_a_3 = zeros(numNodes, numNodes^2)
seq_zeroTenNinety = collect(0:10:90)
for n = 1: numNodes
    mat_a_3[n, (10*n-9): (10*n)] = repeat([1], 10)
    mat_a_3[n, (seq_zeroTenNinety + repeat([n],10))] = repeat([-1], 10)
end
mat_a = vcat(mat_a_1, mat_a_2, mat_a_3)
# Fomulation of Matrix B
mat_b_1 = zeros(numNodes^2, numNodes^2)
for mn = 1: numNodes^2
    mat_b_1[mn, mn] = 170
end
mat_b_2 = zeros(numNodes^2, numNodes^2)
for mn = 1: numNodes^2
    mat_b_2[mn, mn] = 1
end
mat_b_3 = zeros(numNodes, numNodes^2)
mat_b = vcat(mat_b_1, mat_b_2, mat_b_3)
# Fomulation of Vector b
vec_b_1 = zeros(numNodes^2)
vec_b_2 = zeros(numNodes^2)
for m = 1: numNodes
    for n = 1: numNodes
        vec_b_2[10 * (m - 1) + n] = mat_arcTwoNodes[m, n]
    end
end
vec_b_3 = zeros(numNodes)
for n = 1: numNodes
    vec_b_3[n] = vec_netInject[n] - vec_netEject[n]
end
vec_b = vcat(vec_b_1, vec_b_2, vec_b_3)
vec_b = hcat(vec_b)
#
vec_max_y = repeat([1], numNodes^2)
vec_max_y = hcat(vec_max_y)
#
vec_c = zeros(numNodes^2)
for m = 1: numNodes
    for n = 1: numNodes
        vec_c[10 * (m - 1) + n] = mat_distance[m, n]
    end
end
vec_c = hcat(vec_c)
#
vec_f = zeros(numNodes^2)
for m = 1: numNodes
    for n = 1: numNodes
        vec_f[10 * (m - 1) + n] = mat_fixedCost[m, n]
    end
end
vec_f = hcat(vec_f)
#
BendersMilp(n_x = numNodes^2,
            n_y = numNodes^2,
            vec_max_y = vec_max_y,
            vec_c = vec_c,
            vec_f = vec_f,
            vec_b = vec_b,
            mat_a = mat_a,
            mat_b = mat_b,
            epsilon = 0.0001,
            timesIterationMax = 1000)
# obj - sum(mat_arcTwoNodes .* mat_fixedCost)
