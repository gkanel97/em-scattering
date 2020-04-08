function rcs = rcs_distr(control_points, auxiliary_sources, far_field_points, radius)

    % To run it go to /path/to/rcs_dist/for_testing and run
    % ./run_rcs_distr.sh /Applications/MATLAB/MATLAB_Runtime/v94 <input arg>
    
    % Convert input from string to double
    cp = str2double(control_points);
    as = str2double(auxiliary_sources);
    ffp = str2double(far_field_points);

    % Compute the geometry of the problem
    r_cyl = str2double(radius);
    r_sym = 0.8*r_cyl;
    x0 = 0; y0 = 0;
    [x_cp, y_cp] = circular_grid(x0,y0,r_cyl,cp);
    [x_as, y_as] = circular_grid(x0,y0,r_sym,as);
    
    % Compute the coefficients of each auxiliary source
    % by solving the linear system E_s = -E_inc

    lambda = 1;
    k = 2*pi/lambda;
    nu = 0;
    K = 2;

    G = zeros(cp,as);
    for i = 1:cp
        for j = 1:as
            z = k*sqrt((y_cp(i)-y_as(j))^2 + (x_cp(i)-x_as(j))^2);
            G(i,j) = besselh(nu,K,z);
        end
    end

    E = -exp(1j*k*x_cp)';
    A = inv(G)*E;

    % Experimental scattering width

    % We need to define a surface in the far field upon which the far-field
    % control points will be placed
    rho = 100; % radius of the far-field circular surface
    [x_ff, y_ff] = circular_grid(x0,y0,rho,ffp);

    % Compute the electromagnetic field created by the auxiliary sources
    % on each far-field control point
    E_s = zeros(ffp,1);
    for i = 1:ffp
        for j = 1:as
            z = k*sqrt((y_ff(i)-y_as(j))^2 + (x_ff(i)-x_as(j))^2);
            E_s(i) = E_s(i) + A(j)*sqrt(2*1j/(pi*z))*exp(-1j*z);
        end
    end

    scatter_width = 2*pi*rho*(abs(E_s(1))).^2;
    rcs = scatter_width;
    disp(rcs);
end