using Revise
using NewtonRaphson

function f(x, p)

    x_1 = x[1]
    x_2 = x[2]

    f_x = [
        4.5+2*x_1-2*x_2+4*x_1^3-4*x_1*x_2;
        -4+4*x_2-2*x_1-2*x_1^2
    ]

    return f_x
    
end

function j(x, p)

    x_1 = x[1]
    x_2 = x[2]

    J_11 = 2 + 12*x_1^2 - 4*x_2
    J_12 = -2 - 4*x_1
    J_21 = J_12
    J_22 = 4

    return [J_11 J_12; J_21 J_22]

end

function h(x, p)

    x_1 = x[1]
    x_2 = x[2]

    H_11 = 12*x_1^2 - 4*x_2 + 2
    H_12 = -4*x_1 - 2
    H_21 = H_12
    H_22 = 4

    return [H_11 H_12; H_21 H_22]
end

function objective(x, p)

    x_1 = x[1]
    x_2 = x[2]

    return 4 + 4.5*x_1 - 4*x_2 + x_1^2 + 2*x_2^2 - 2*x_1*x_2 + x_1^4 - 2*x_1^2*x_2

end

x_0 = ([0.5; 0.5], [-10.0; 10.0],[2.0; 4.5])
p = [0.5 0.5]
x_opt = find_global_min(f, j, h, objective, x_0)

