function J = objective_fun(u, x0, x_eq, Q, R, p, Np, Ts)
    J = 0;
    x_ref = x_eq;
    for k = 1:Np
        [~, x] = ode45(@(t, x) tumor_growth_controlled(t, x, u, p), [0 Ts], x0); 
        x0 = x(end, :);
        J = J +  (x0 - x_ref') * Q * (x0 - x_ref')' +  R * u^2;  
    end
end
