function dxdt = tumor_growth_controlled(t, x, u, p)
    r1 = p(1);  r2 = p(2);  b1 = p(3);  b2 = p(4);
    c1 = p(5);  c2 = p(6);  c3 = p(7);  c4 = p(8);
    rho = p(9); alpha = p(10); s = p(11);
    d1 = p(12); gamma1 = p(13); gamma2 = p(14); 
    gamma3 = p(15); k = p(16);

    x1 = x(1);       % (normal cells)
    x2 = x(2);       % (tumor cells)
    x3 = x(3);       % (immune cells)
    x4 = x(4);       % (drug concentration)
    
    dx1dt = r2 * x1 * (1 - b2 * x1) - c4 * x2 * x1 - gamma1 * x4;
    dx2dt = r1 * x2 * (1 - b1 * x2) - c2 * x3 * x2 - c3 * x2 * x1 - gamma2 * x4;
    dx3dt = s + (rho * x3 * x2) / (alpha + x2) - c1 * x3 * x2 - d1 * x3 - gamma3 * x4;
    dx4dt = u - k * x4;
    
    dxdt = [dx1dt; dx2dt; dx3dt; dx4dt];
end

