function err_max = max_err_distr(control_points, auxiliary_sources, error_points, radius)
    
    % To run it go to /path/to/rcs_dist/for_testing and run
    % ./run_max_err_distr.sh /Applications/MATLAB/MATLAB_Runtime/v94 <input arg>

    % Initializations
    x0 = 0; y0 = 0;
    lambda = 1;
    cp = str2double(control_points);
    as = str2double(auxiliary_sources);
    errp = str2double(error_points);
    
    % Define the geometry of the scattering object
    r_cyl = 1; % normalized cylinder radius
    r_sym = str2double(radius); % normalized radius of the auxiliary surface

    % Compute the geometry of the problem
    x0 = 0; y0 = 0;
    [x_cp, y_cp] = circular_grid(x0,y0,r_cyl,cp);
    [x_as, y_as] = circular_grid(x0,y0,r_sym,as);

    % Compute the coefficients of each auxiliary source
    % by solving the linear system E_s = -E_inc

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
    A = inv(G) * E;

    % The error points are on the surface of the cylinder
    [x_ep, y_ep] = circular_grid(x0,y0,r_cyl,errp);

    E_inc = exp(1j*k*x_ep)'; % incident wave

    % Compute the electromagnetic field created by the auxiliary sources
    % on each error point
    E_s = zeros(errp,1);
    for i = 1:errp
        for j = 1:as
            z = k*sqrt((y_ep(i)-y_as(j))^2 + (x_ep(i)-x_as(j))^2);
            E_s(i) = E_s(i) + A(j)*besselh(nu,K,z);
        end
    end

    % Boundary condition: E_s + E_inc = 0 on the cylinder surface
    error = E_inc + E_s;
    err_max = max(abs(error));
    disp(err_max);
end
    
 