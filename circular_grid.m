function [x,y] = circular_grid(x0,y0,r,points)
    phi = linspace(0, 2*pi*(1-1/points), points);
    x = x0 + r*cos(phi);
    y = y0 + r*sin(phi);
end

