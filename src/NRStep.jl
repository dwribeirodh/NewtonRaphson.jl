
struct Newton
    f::Function
    j::Function
    p
    x_0
    ϵ_1
    ϵ_2
    max_steps
end

function Newton(f::Function, j::Function, p, x_0; ϵ_1 = 1E-6, ϵ_2 = 1E-6, max_steps = 500)
    return Newton(f, j, p, x_0, ϵ_1, ϵ_2, max_steps)
end

struct NewtonSolution
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

    return x + p_k
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
    ϕ = system.ϵ_1 + 1
    vec_norm = system.ϵ_2 + 1



    error_tolerance = ϕ <= system.ϵ_1 && vec_norm <= system.ϵ_2

    while !error_tolerance || n_steps <= system.max_steps



    end



end
