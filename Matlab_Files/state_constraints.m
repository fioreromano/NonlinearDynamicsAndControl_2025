function [c, ceq] = state_constraints(u, x0_new, p, x_max, x_min)
    [~, x] = ode45(@(t, x) tumor_growth_controlled(t, x, u, p), [0 1], x0_new);
    
    c = [x_min - x(end, :)';  % Deve essere >= 0
         x(end, :)' - x_max]; % Deve essere <= 0
    
    ceq = []; % Nessun vincolo di uguaglianza
end