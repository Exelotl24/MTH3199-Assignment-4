function [h_list, analytical_difference, e_local] = local_truncation_error(rate_func_in, solution, t_ref, BT_struct)
    X0 = solution(t_ref);

    h_list = logspace(-5, 1, 25);
    e_local = zeros(size(h_list));
    analytical_difference = zeros(size(h_list));


    for i = 1:length(h_list)
        X = solution(t_ref+h_list(i));
        analytical_difference(i) = norm(X-X0);
        e_local(i) = norm(explicit_RK_step(rate_func_in,t_ref,X0,h_list(i),BT_struct)-X);
    end

end