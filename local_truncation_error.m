function [h_list, analytical_difference, e_local_1, e_local_2, e_local_3] = local_truncation_error(rate_func_in, solution, t_ref, BT_struct1, BT_struct2, BT_struct3)
    X0 = solution(t_ref);

    h_list = logspace(-5, 1, 25);
    e_local_1 = zeros(size(h_list));
    e_local_2 = zeros(size(h_list));
    e_local_3 = zeros(size(h_list));
    analytical_difference = zeros(size(h_list));


    for i = 1:length(h_list)
        X = solution(t_ref+h_list(i));
        analytical_difference(i) = norm(X-X0);
        e_local_1(i) = norm(explicit_RK_step(rate_func_in,t_ref,X0,h_list(i),BT_struct1)-X);
        e_local_2(i) = norm(explicit_RK_step(rate_func_in,t_ref,X0,h_list(i),BT_struct2)-X);
        e_local_3(i) = norm(explicit_RK_step(rate_func_in,t_ref,X0,h_list(i),BT_struct3)-X);
    end

end