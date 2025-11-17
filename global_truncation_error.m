function [h_list, e_global] = global_truncation_error(rate_func_in, solution, tspan, BT_struct)
    X0 = solution(tspan(1));
    
    delta_t = abs(tspan(end)-tspan(1));
    h_list = logspace(-5, -1, 25)*delta_t;
    e_global = zeros(size(h_list));


    for i = 1:length(h_list)
        [t_list,X_list,~, ~] = explicit_RK_fixed_step_integration(rate_func_in,tspan,X0,h_list(i),BT_struct);
        X = solution(t_list(end));
        
        e_global(i) = norm(X_list(end,:)'-X);
    end

end