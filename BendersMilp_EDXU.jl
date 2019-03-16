# Benders Algorithm to Solve of Mixed Integer Linear Programming more Quickly
# Edward J. Xu, edxu96@outlook.com
# March 16th, 2019
module BendersMilp_EDXU
    export BendersMilp
    using JuMP
    using GLPKMathProgInterface
    function BendersMILP(; n_x, n_y, vec_max_y, vec_c, vec_f, vec_b, mat_a, mat_b, epsilon, timesIterationMax)
        """
        Benders Algorithm to Solve of Mixed Integer Linear Programming more Quickly

            Example:
            push!(LOAD_PATH, "$(homedir())/Desktop/...")
            cd("$(homedir())/Desktop/...")
            using BendersMILP_EDXU
            using JuMP
            using GLPKMathProgInterface
            n_x = 1
            n_y = 1
            vec_max_y = [10]
            vec_c = [5]
            vec_f = [-3]
            vec_b = [4; 0; -13]
            mat_a = [1; 2; 1]
            # mat_b must be defined strictly as a column vector.
            mat_b = rand(3,1)
            mat_b[1] = 2
            mat_b[2] = -1
            mat_b[3] = -3
            BendersMILP(n_x = n_x,
                        n_y = n_y,
                        vec_max_y = vec_max_y,
                        vec_c = vec_c,
                        vec_f = vec_f,
                        vec_b = vec_b,
                        mat_a = mat_a,
                        mat_b = mat_b,
                        epsilon = 0,
                        timesIterationMax = 5)

            Version: 1.0
            Date: 190316
            Edward J. Xu
        """
        n_constraint = length(mat_a[:, 1])
        let
            model_mas = Model(solver = GLPKSolverMIP())
            @variable(model_mas, q)
            @variable(model_mas, vec_y[1: n_y] >= 0, Int)
            @constraint(model_mas, vec_y .<= vec_max_y)
            # add master-problem objective
            @objective(model_mas, Min, transpose(vec_f) * vec_y + q)
            boundUp = Inf
            boundLow = - Inf
            # definition of model_sub variables
            vec_uBar = zeros(n_constraint, 1)
            # initial value of master variables
            vec_yBar = zeros(n_y, 1)
            #
            obj_mas = 0
            obj_sub = 0
            timesIteration = 1
            while (boundUp - boundLow > epsilon && timesIteration < timesIterationMax)
                # model_sub-problem
                let
                    model_sub = Model(solver = GLPKSolverLP())
                    # create full model_sub-problem model, using the value of master variables ybar
                    @variable(model_sub, vec_u[1: n_constraint] >= 0)
                    @objective(model_sub, Max, (transpose(vec_b - mat_b * vec_yBar) * vec_u)[1])
                    @constraint(model_sub, transpose(mat_a) * vec_u .<= vec_c)
                    solve(model_sub)
                    obj_sub = getobjectivevalue(model_sub)
                    for i = 1: n_constraint
                        vec_uBar[i] = getvalue(vec_u[i])
                    end
                end
                boundUp = min(boundUp, obj_sub + (transpose(vec_f) * vec_yBar)[1])
                # master-problem only add the feasibility constraint
                # If mat_b is not defined strictly as a column vector, there will be a problem here.
                @constraint(model_mas, (transpose(vec_uBar) * (vec_b - mat_b * vec_y))[1] <= q)
                solve(model_mas)
                obj_mas = getobjectivevalue(model_mas)
                for i = 1: n_y
                    vec_yBar[i] = getvalue(vec_y[i])
                end
                boundLow = max(obj_mas, boundLow)
                println("timesIteration: $(timesIteration)\t boundUp: $(boundUp)\t boundLow: $(boundLow)")
                timesIteration += 1
            end
        end
    end
    println("Correct Ending")
end
