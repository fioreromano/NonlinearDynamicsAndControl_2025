function Vdot = V_dot_fun(xi, t, xi_eq, r1, r2, b1, b2, c1, c2, c3, c4, rho, alpha, s, d1)

P = [30 0 0; 0 70 0; 0 0 1];

Vdot = (xi-xi_eq)'* P *tumor_growth(t, xi, r1, r2, b1, b2, c1, c2, c3, c4, rho, alpha, s, d1);

end