# Benders template with rays
# Edward J. Xu, edxu96@outlook.com
# March 16th, 2019
module BendersMilp_EDXU
    export BendersMilp
    using JuMP
    using GLPKMathProgInterface
    function BendersMilp(; n_x, n_y, vec_max_y, vec_c, vec_f, vec_b, mat_a, mat_b, epsilon, timesIterationMax)
        # Define Master problem
        n_constraint = length(mat_a[:, 1])
        model_mas = Model(solver = GLPKSolverMIP())
        @variable(model_mas, q)
        @variable(model_mas, vec_y[1: n_y] >= 0, Int)
        @objective(model_mas, Min, (transpose(vec_f) * vec_y + q)[1])
        @constraint(model_mas, vec_y[1: n_y] .<= vec_max_y)


        function solve_master(vec_uBar, opt_cut::Bool)
            if opt_cut
                @constraint(model_mas, (transpose(vec_uBar) * (vec_b - mat_b * vec_y))[1] <= q)
            else  # Add feasible cut Constraints
                @constraint(model_mas, (transpose(vec_uBar) * (vec_b - mat_b * vec_y))[1] <= 0)
            end
            @constraint(model_mas, (transpose(vec_uBar) * (vec_b - mat_b * vec_y))[1] <= q)
            solve(model_mas)
            return getobjectivevalue(model_mas)
        end


        function solve_sub(vec_uBar, vec_yBar, n_constraint, vec_b, mat_b, mat_a, vec_c)
            model_sub = Model(solver = GLPKSolverLP())
            @variable(model_sub, vec_u[1: n_constraint] >= 0)
            @objective(model_sub, Max, (transpose(vec_b - mat_b * vec_yBar) * vec_u)[1])
            constraintsForDual = @constraint(model_sub, transpose(mat_a) * vec_u .<= vec_c)
            solution_sub = solve(model_sub)
            print("------------------------------ Sub Problem ------------------------------\n")  # , model_sub)
            vec_uBar = getvalue(vec_u)
            if solution_sub == :Optimal
                vec_result_x = zeros(length(vec_c))
                vec_result_x = getdual(constraintsForDual)
                return (true, getobjectivevalue(model_sub), vec_uBar, vec_result_x)
            end
            if solution_sub == :Unbounded
                return (false, getobjectivevalue(model_sub), vec_uBar, repeat([NaN], length(vec_c)))
            end
        end


        function solve_ray(vec_uBar, vec_yBar, n_constraint, vec_b, mat_b, mat_a)
            # model_ray = Model(solver = GurobiSolver())
            model_ray = Model(solver = GLPKSolverLP())
            @variable(model_ray, vec_u[1: n_constraint] >= 0)
            @objective(model_ray, Max, 1)
            @constraint(model_ray, (transpose(vec_b - mat_b * vec_yBar) * vec_u)[1] == 1)
            @constraint(model_ray, transpose(mat_a) * vec_u .<= 0)
            solve(model_ray)
            print("------------------------------ Ray Problem ------------------------------\n")  # , model_ray)
            vec_uBar = getvalue(vec_u)
            obj_ray = getobjectivevalue(model_ray)
            return (obj_ray, vec_uBar)
        end


        # Begin Calculation
        let
            boundUp = Inf
            boundLow = - Inf
            epsilon = 0
            # initial value of master variables
            vec_uBar = zeros(n_constraint, 1)
            vec_yBar = zeros(n_y, 1)
            vec_result_x = length(n_x)
            obj_sub = 0
            timesIteration = 1
            while (boundUp - boundLow > epsilon)
                (bool_solutionSubModel, obj_sub, vec_uBar, vec_result_x) = solve_sub(vec_uBar, vec_yBar, n_constraint,
                                                                                     vec_b, mat_b, mat_a, vec_c)
                if bool_solutionSubModel
                    boundUp = min(boundUp, obj_sub + (transpose(vec_f) * vec_yBar)[1])
                else
                    (obj_ray, vec_uBar) = solve_ray(vec_uBar, vec_yBar, n_constraint, vec_b, mat_b, mat_a)
                end
                obj_mas = solve_master(vec_uBar, bool_solutionSubModel)
                vec_yBar = getvalue(vec_y)
                boundLow = max(obj_mas, boundLow)
                if bool_solutionSubModel
                    println("------------------- Result in $(timesIteration)-th Iteration with Sub",
                            "-------------------\n", "boundUp: $(boundUp), boundLow: $(boundLow), Sub: $(obj_sub).")
                else
                    println("------------------- Result in $(timesIteration)-th Iteration with Ray",
                            "-------------------\n", "boundUp: $(boundUp), boundLow: $(boundLow), Sub: $(obj_ray).")
                end
                timesIteration += 1
            end
            println("----------------------------- Master Problem ----------------------------\n", model_mas)
            print("-------------------------------------------------------------------------\n",
                  "--------------------------------- Result --------------------------------\n",
                  "-------------------------------------------------------------------------\n")
            println("boundUp: $(boundUp), boundLow: $(boundLow), difference: $(boundUp - boundLow)")
            println("(final) obj_sub: $(obj_sub)")
            vec_result_y = getvalue(vec_y)
            println("vec_y: ", vec_result_y)
            println("vec_x: ", vec_result_x)
        end
        println("-------------------------------------------------------------------------\n",
                "----------------------------- Correct Ending ----------------------------\n",
                "-------------------------------------------------------------------------\n")
    end
end
