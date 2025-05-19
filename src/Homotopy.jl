# detects weight for which we have sagbi

# ------------------- Input:
# poss_sagbi          Union{Vector{Vector{Expression}}
# ------------------- Output:
# w                   Vector{Int} 

function get_weight(poss_sagbi::Union{Vector{Vector{Expression}}, Vector{Expression}} )
    
    if typeof(poss_sagbi) == Vector{Expression}
        poss_sagbi = [poss_sagbi]
    end

    l = length(poss_sagbi)
    #vars = unique(vcat(variables.(poss_sagbi)...))
    vars = variables(vcat(poss_sagbi...))

    @unique_var t[1:l] 
    vars2 = tuple(vcat(t,vars)...)
    R, vars2 = Oscar.polynomial_ring(Oscar.QQ, vcat(["t$i" for i=1:l ],string.(vars)))
    
    poss_sagbiOscar = vcat([vars2[i].*([HC_to_oscar(f, gens(R)[l+1:end], vars) for f in poss_sagbi[i]]) for i=1:l ]...)

    w = weightVectorsRealizingSAGBI(poss_sagbiOscar, R)

    if typeof(w) == Nothing 
        return w
    end 

    return (w[l+1:end])
end


# solve a linear system in a sagbi basis

# ---------------------- Input:
# lin_sys                Union{Vector{Vector{Expression}}
# poss_sagbi             Union{Vector{Vector{Expression}}
# weight                 Vector{Int} (optional)
# degreeCheck = true     (optional)
# getBaseLocus = false   (optional)
# varyLinearPart = false (optional)
# ---------------------- Output:
# result                 
# sols                   

function sagbi_homotopy(lin_sys::Union{Vector{Vector{Expression}}, Vector{Expression}} , poss_sagbi::Union{Vector{Vector{Expression}}, Vector{Expression}} ; weight = nothing, degreeCheck = true, getBaseLocus = false, varyLinearPart = false)
    
    if typeof(lin_sys) == Vector{Expression}
        lin_sys = [lin_sys]
    end

    if typeof(poss_sagbi) == Vector{Expression}
        poss_sagbi = [poss_sagbi]
    end

    if (length(poss_sagbi) != length(lin_sys) )
        return "Error: length of linear equations does not match length of parameterization"
    end

    #detecting a weight
    if weight == nothing
        weight = get_weight(poss_sagbi)
        if weight == nothing
            return "Error: polynomials do not form a SAGBI basis."
        end
    end 


    if degreeCheck
        #degree of maps check
        d1 = degree_map(poss_sagbi) 
        d2 = degree_monomial_map(poss_sagbi, weight)
        if (d1 - d2 > 0)
            println("Error: degree of monomial parameterisation drops from ", d1, " to ", d2, ". SAGBI homotopy will not find all the solutions.")
        end
    end

    
    vars_lin_sys = variables.(lin_sys)
    #vars_sagbi = unique(vcat(variables.(poss_sagbi)...))
    vars_sagbi = variables(vcat(poss_sagbi...))

    
    @unique_var t
    deformed_sagbi = [ map(poly -> weight_deformation_for_poly(poly, vars_sagbi, weight, t = t), sagbi) for sagbi in poss_sagbi] 
    randomnumber = randn(ComplexF64)

    if varyLinearPart
        new_lin_sys = []
        for i in 1:length(lin_sys)
            var = variables(lin_sys[i])

            d = length(lin_sys[i])
            n = length(var)
            A0 = randn(ComplexF64,d,n)
            push!(new_lin_sys, (1-t).*(A0*var) + randomnumber*(t).*lin_sys[i])
        end
        
        new_sys = vcat([ subs(new_lin_sys[i], vars_lin_sys[i] => deformed_sagbi[i]) for i=1:length(new_lin_sys) ]...)
    else
        new_sys = vcat([ subs(lin_sys[i], vars_lin_sys[i] => deformed_sagbi[i]) for i=1:length(lin_sys) ]...)
    end

    F = System(new_sys, variables = vars_sagbi, parameters = [t])
    start_solutions = solutions(HomotopyContinuation.solve(F; target_parameters = [0]))
    
    result = HomotopyContinuation.solve(F, start_solutions; start_parameters = [0], target_parameters = [1])
    
    sols = solutions(result)
    
    if getBaseLocus
        base_sols = []
        for sagbi in poss_sagbi
            num_missing_sols = isolated_nsols_baselocus(sagbi)

            if num_missing_sols > 0
                missing_sols = get_base_locus(sagbi)
                for missing_sol in missing_sols
                    if norm( HomotopyContinuation.evaluate(F, missing_sol, [1]) , Inf) < 1e-10
                        push!(base_sols, missing_sol)
                    end
                end 
            end
        end
        base_sols = unique(base_sols)
        sols = vcat(sols, base_sols)
    end

    println("SAGBI homotopy successfully completed with ", length(sols), " solutions.")
    return (result, sols)
end