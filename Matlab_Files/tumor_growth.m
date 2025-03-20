function dxdt = tumor_growth(t, x, r1, r2, b1, b2, c1, c2, c3, c4, rho, alpha, s, d1)

    x1 = x(1);       % (normal cells)
    x2 = x(2);       % (tumor cells)
    x3 = x(3);       % (immune cells)
    
    dx1dt = r2 * x1 * (1 - b2 * x1) - c4 * x2 * x1;
    dx2dt = r1 * x2 * (1 - b1 * x2) - c2 * x3 * x2 - c3 * x2 * x1;
    dx3dt = s + (rho * x3 * x2) / (alpha + x2) - c1 * x3 * x2 - d1 * x3;
    
    dxdt = [dx1dt; dx2dt; dx3dt];
end
