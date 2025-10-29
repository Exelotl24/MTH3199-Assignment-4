function J = approximate_jacobian(fun,X)
    % initialize Jacobian
    f0 = fun(X);
    J = zeros(length(f0), length(X));

    % initialize matrix
    e_n = zeros(length(X), 1);

    % define step size
    delta_X = 1e-6;

    % run through each point in vector X
    for n = 1:length(X)
        % set partial derivative to look at current point
        e_n(n) = 1;
        
        % calculate Jacobian numerator
        f_plus = fun(X+e_n*delta_X);
        f_minus = fun(X-e_n*delta_X);

        % update Jacobian matrix
        J(:,n) = (f_plus-f_minus)/(2*delta_X);

        e_n(n) = 0;

    end
end