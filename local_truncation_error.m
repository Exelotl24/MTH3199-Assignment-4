function [h_list, analytical_difference, e_local_FE, e_local_MP, e_local_BE, e_local_IMP] = local_truncation_error(rate_func_in, solution, t_ref)
    X0 = solution(t_ref);

    h_list = logspace(-5, 1, 25);
    e_local_FE = zeros(size(h_list));
    e_local_MP = zeros(size(h_list));
    e_local_BE = zeros(size(h_list));
    e_local_IMP = zeros(size(h_list));
    analytical_difference = zeros(size(h_list));

    for i = 1:length(h_list)
        X = solution(t_ref+h_list(i));
        analytical_difference(i) = norm(X-X0);
        e_local_FE(i) = norm(forward_euler_step(rate_func_in, t_ref, X0, h_list(i))-X);
        e_local_MP(i) = norm(explicit_midpoint_step(rate_func_in, t_ref, X0, h_list(i))-X);
        e_local_BE(i) = norm(backward_euler_step(rate_func_in, t_ref, X0, h_list(i))-X);
        e_local_IMP(i) = norm(implicit_midpoint_step(rate_func_in, t_ref, X0, h_list(i))-X);
    end

end