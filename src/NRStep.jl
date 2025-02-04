
struct Newton
    f::Function
    j::Function
    h::Function
    p
    x_0
    ϵ_1
    ϵ_2
    max_steps
end

function Newton(f::Function, j::Function, h::Function, p, x_0; ϵ_1 = 1E-6, ϵ_2 = 1E-6, max_steps = 500)
    return Newton(f, j, h, p, x_0, ϵ_1, ϵ_2, max_steps)
end

struct NewtonSolution
    converged
    x_final
    n_steps
    ϕ
    vec_norm
end

function calculate_jacobian(system::Newton, x)
    return system.j(x, system.p)
end

function calculate_function(system::Newton, x)
    return system.f(x, system.p)
end

function step(system::Newton, x)

    J_k = calculate_jacobian(system, x)
    f_k = calculate_function(system, x)

    p_k = J_k \ -f_k

    return p_k

end

function update_x(system::Newton, x)

    p_k = step(system, x)

    x_new = x + p_k
    x_new = reshape(x_new, length(x_new), 1)

    return x_new
end


function calculate_ϕ(system::Newton, x)

    f_k = system.f(x, system.p)

    return 0.5 * dot(f_k, f_k)

end

function calculate_vec_norm(x, x_k)

    soln_error = x - x_k
    sqrs_sum = sum([i^2 for i in soln_error])

    return sqrt(sqrs_sum)
    
end

function solve(system::Newton)

    n_steps = 1
    x = update_x(system, system.x_0)
    ϕ = calculate_ϕ(system, x)
    vec_norm = calculate_vec_norm(x, system.x_0)
    converged = false

    error_satisfied = ϕ <= system.ϵ_1 && vec_norm <= system.ϵ_2

    while !error_satisfied

        if n_steps == system.max_steps
            converged = false
            break
        end
        
        n_steps += 1
        
        x_k = x
        x = update_x(system, x)
        ϕ = calculate_ϕ(system, x)
        vec_norm = calculate_vec_norm(x, x_k)

        error_satisfied = ϕ <= system.ϵ_1 && vec_norm <= system.ϵ_2

    end

    return NewtonSolution(true, x, n_steps, ϕ, vec_norm)

end

function is_positive_def(system::Newton, x, p)

    H_x = system.h(x, p)
    λ = eigvals(H_x)

    if all(>(0), λ)
        return true
    else
        return false
    end

end

function find_global_min(f::Function, j::Function, h::Function, objective::Function, x_0)
    
    minima = []
    p = [0.0;0.0]

    for init_guess in x_0
        system = Newton(f, j, h, p, init_guess)
        x = solve(system)
        is_pos_def = is_positive_def(system, x.x_final, p)
        if is_pos_def
            push!(minima, x.x_final)
        end
    end

    f_x_opt = zeros(Float64, (length(minima), 1))
    for (idx,local_min) in enumerate(minima)
        f_x_opt[idx] = objective(local_min, p)
    end

    min_idx = findmin(f_x_opt)[2]

    return minima[min_idx]

end