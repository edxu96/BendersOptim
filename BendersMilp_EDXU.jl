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


        function solve_master(vec_uBar, opt_cut::Bool )
            if opt_cut
                @constraint(model_mas, (transpose(vec_uBar) * (vec_b - mat_b * vec_y))[1] <= q)
            else  # Add feasible cut Constraints
                @constraint(model_mas, (transpose(vec_uBar) * (vec_b - mat_b * vec_y))[1] <= 0)
            end
            @constraint(model_mas, (transpose(vec_uBar) * (vec_b - mat_b * vec_y))[1] <= q)
            solve(model_mas)
            return getobjectivevalue(model_mas)
        end


        function solve_sub(vec_uBar, vec_yBar, n_constraint, vec_b, mat_b, mat_a)
            model_sub = Model(solver = GLPKSolverLP())
            @variable(model_sub, vec_u[1: n_constraint] >= 0)
            @objective(model_sub, Max, (transpose(vec_b - mat_b * vec_yBar) * vec_u)[1])
            @constraint(model_sub, transpose(mat_a) * vec_u .<= vec_c)
            solution_sub = solve(model_sub)
            print("------------------------------ Sub Problem ------------------------------\n",
                  model_sub)
            for i = 1: n_constraint
                vec_uBar[i] = getvalue(vec_u[i])
            end
            if solution_sub == :Optimal
                return (true, getobjectivevalue(model_sub), vec_uBar)
            end
            if solution_sub == :Unbounded
                return (false, getobjectivevalue(model_sub), vec_uBar)
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
            print("------------------------------ Ray Problem ------------------------------\n", model_ray)
            for i = 1: n_constraint
                vec_uBar[i] = getvalue(vec_u[i])
            end
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
            timesIteration = 1
            while (boundUp - boundLow > epsilon)
                (bool_solutionSubModel, sub_obj, vec_uBar) = solve_sub(vec_uBar, vec_yBar, n_constraint,
                                                                       vec_b, mat_b, mat_a)
                if bool_solutionSubModel
                    boundUp = min(boundUp, sub_obj + (transpose(vec_f) * vec_yBar)[1])
                else
                    (obj_ray, vec_uBar) = solve_ray(vec_uBar, vec_yBar, n_constraint, vec_b, mat_b, mat_a)
                end
                obj_mas = solve_master(vec_uBar, bool_solutionSubModel)
                for i = 1: n_y
                    vec_yBar[i] = getvalue(vec_y[i])
                end
                boundLow = max(obj_mas, boundLow)
                if bool_solutionSubModel
                    println("-------------------- $(timesIteration)-th Iteration with Sub Problem ",
                            "--------------------\n", "boundUp: $(boundUp), boundLow: $(boundLow), Sub: $(sub_obj).")
                else
                    println("-------------------- $(timesIteration)-th Iteration with Ray Problem ",
                            "--------------------\n", "boundUp: $(boundUp), boundLow: $(boundLow), Sub: $(obj_ray).")
                end
                timesIteration += 1
            end
        end
        println("----------------------------- Master Problem ----------------------------\n", model_mas)
        #-------------------------------------------------------------------
        println("----------------------------- Correct Ending ----------------------------\n")
    end
end
