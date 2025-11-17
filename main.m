clear
close all

main

function main()


    % Create Orbit Parameters
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    gravity_rate_func_wrapper = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);

    % Define Planet start
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    % Create BT struct for Forward Euler
    BT_struct_FE = struct();
    BT_struct_FE.C = [0];
    BT_struct_FE.B = [1];
    BT_struct_FE.A = [0];

    % Create BT struct for Explicit Midpoint
    BT_struct_EM = struct();
    BT_struct_EM.C = [0; 0.5];
    BT_struct_EM.B = [0, 1];
    BT_struct_EM.A = [0, 0; 0.5, 0];


    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
        5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
        1/5, 0, 0, 0,0,0,0;...
        3/40, 9/40, 0, 0, 0, 0,0;...
        44/45, -56/15, 32/9, 0, 0, 0,0;...
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
    Fehlberg = struct();
    Fehlberg.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
    Fehlberg.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
        25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    Fehlberg.A = [0,0,0,0,0,0;...
        1/4, 0,0,0,0,0;...
        3/32, 9/32, 0,0,0,0;...
        1932/2197, -7200/2197, 7296/2197, 0,0,0;...
        439/216, -8, 3680/513, -845/4104, 0,0;...
        -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
    HeunEuler = struct();
    HeunEuler.C = [0,1];
    HeunEuler.B = [1/2,1/2;1,0];
    HeunEuler.A = [0,0;1,0];
    FehlbergRK1 = struct();
    FehlbergRK1.C = [0,1/2,1];
    FehlbergRK1.B = [1/512, 255/256, 1/512;...
        1/256, 255/256, 0];
    FehlbergRK1.A = [0,0,0;1/2,0,0;1/256,255/256,0];
    Bogacki = struct();
    Bogacki.C = [0,1/2, 3/4, 1];
    Bogacki.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
    Bogacki.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3, 4/9, 0];

    ti = 0;             % starting time
    tf = 300;            % ending time

    tspan = [ti, tf];

    % Define X0 as starting point location & velocity of planet
    X0 = [x0;y0;dxdt0;dydt0];
    t_exact = linspace(tspan(1),tspan(2),300);
    V_list = compute_planetary_motion(t_exact,X0,orbit_params);
    
    
    
    h_ref= 0.1;
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(gravity_rate_func_wrapper,tspan,X0,h_ref,BT_struct_EM);

    figure(1);
    subplot(2,1,1)
    hold on
    plot(t_list,X_list(:,1),'r', 'DisplayName', 'Calculated X-pos');
    plot(t_list,X_list(:,2),'b', 'DisplayName', 'Calculated Y-pos');
    plot(t_exact,V_list(:,1),'k--', 'DisplayName', 'Exact X-pos');
    plot(t_exact,V_list(:,2),'k--', 'DisplayName', 'Exact Y-pos');

    xlabel('Time')
    ylabel('Position')
    title('Position Fixed RK Step Orbit')
    legend()

    subplot(2,1,2)
    hold on
    plot(t_list,X_list(:,3),'r', 'DisplayName', 'Calculated X-vel');
    plot(t_list,X_list(:,4),'b', 'DisplayName', 'Calculated Y-vel');
    plot(t_exact,V_list(:,3),'k--', 'DisplayName', 'Exact X-vel');
    plot(t_exact,V_list(:,4),'k--', 'DisplayName', 'Exact Y-vel');

    xlabel('Time')
    ylabel('Velocity')
    title('Velocity Fixed RK Step Orbit')
    legend()

    figure(2)
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(X_list(:,1),X_list(:,2),'r', 'DisplayName', 'Fixed RK');
    plot(V_list(:,1),V_list(:,2),'k--', 'DisplayName', 'Exact Orbit');
    title('Fixed RK Orbit vs Exact Solution')
    xlabel('Position (x)')
    ylabel('Position (y)')
    legend()


    % ------------------- LOCAL TRUNCATION ---------------------
    V0 = [x0;y0;dxdt0;dydt0];
    solution_func = @(t_in) compute_planetary_motion(t_in, V0, orbit_params);

    HeunEuler = struct();
    HeunEuler.C = [0,1];
    HeunEuler.B = [1/2,1/2;1,0];
    HeunEuler.A = [0,0;1,0];


    t_ref = 0.5;
    [h_list, analytical_difference, e_local_FE] = local_truncation_error(gravity_rate_func_wrapper, solution_func, t_ref, BT_struct_FE);
    [~, ~, e_local_MP] = local_truncation_error(gravity_rate_func_wrapper, solution_func, t_ref, BT_struct_EM);
    [~, ~, e_local_Heun] = local_truncation_error(gravity_rate_func_wrapper, solution_func, t_ref, HeunEuler);


    figure(3)
    loglog(h_list,analytical_difference, 'DisplayName', 'Analytical Difference');
    hold on
    loglog(h_list, e_local_FE, 'DisplayName', 'Forward Euler Error')
    loglog(h_list, e_local_MP, 'DisplayName', 'Explicit Midpoint Error')
    loglog(h_list, e_local_Heun, 'DisplayName', 'Heun Euler')
    title("Local Truncation Error vs. Step Size")
    xlabel('step size')
    ylabel('local truncation error')
    legend


    % % -------------- GLOBAL TRUNCATION ERROR -----------------
    % global_tspan = [0, 300];
    % 
    % disp('calculating global FE')
    % [h_list, e_global_FE] = global_truncation_error(gravity_rate_func_wrapper, solution_func, global_tspan, BT_struct_FE);
    % disp('calculating global MP')
    % [~, e_global_MP] = global_truncation_error(gravity_rate_func_wrapper, solution_func, global_tspan, BT_struct_EM);
    % disp('calculating global Heun')
    % [~, e_global_Heun] = global_truncation_error(gravity_rate_func_wrapper, solution_func, global_tspan, HeunEuler);
    % disp('calculating global DormandPrince')
    % [~, e_global_DP] = global_truncation_error(gravity_rate_func_wrapper, solution_func, global_tspan, DormandPrince);
    % disp('calculating global Fehlberg')
    % [~, e_global_Fehlberg] = global_truncation_error(gravity_rate_func_wrapper, solution_func, global_tspan, Fehlberg);
    % 
    % 
    % figure(4)
    % loglog(h_list, e_global_FE, 'DisplayName', 'Forward Euler Error')
    % hold on
    % loglog(h_list, e_global_MP, 'DisplayName', 'Explicit Midpoint Error')
    % loglog(h_list, e_global_Heun, 'DisplayName', 'Heun Euler')
    % loglog(h_list, e_global_DP, 'DisplayName', 'Dormand Prince')
    % loglog(h_list, e_global_Fehlberg, 'DisplayName', 'Fehlberg')
    % title("Global Truncation Error vs. Step Size")
    % xlabel('step size')
    % ylabel('global truncation error')
    % legend()


     % -------------- CONSERVATION OF PHYSICAL QUALITIES ----------------

     [EV, HV] = conservation(t_exact, V_list, orbit_params);        % figure(5)
     [EX, HX] = conservation(t_list, X_list, orbit_params);         % figure(6)

     % -------------- LOCAL TRUNCATION -----------------------

     h_list = logspace(-4, -1, 100);


    ti = 0;             % starting time
    tf = 300;           % ending time
    
    tspan = [ti, tf];
    
    % Define X0 as starting point location & velocity of planet
    X0 = [x0;y0;dxdt0;dydt0];


    
    % [XB1, XB2, ~] = RK_step_embedded(gravity_rate_func_wrapper,ti,X0,0.1,FehlbergRK1)

    h_list = logspace(-2,1,50);

%     XB1_list = zeros(4, length(h_list));
%     XB2_list = zeros(4, length(h_list));

    for i= 1:length(h_list)
%         [XB1, XB2, ~] = RK_step_embedded(gravity_rate_func_wrapper,ti,X0,h_list(i),DormandPrince);
%         XB1_list(:, i) = XB1;
%         XB2_list(:, i) = XB2;

        X0 = solution_func(t_ref);
        X = solution_func(t_ref+h_list(i));
        [XB1, XB2, ~] = RK_step_embedded(gravity_rate_func_wrapper,ti,X0,h_list(i),DormandPrince);

        analytical_difference(i) = norm(X-X0);
        e_local_XB1(i) = norm(XB1-X);
        e_local_XB2(i) = norm(XB2-X);
        XB_diff(i) = norm(XB2-XB1);
    end

    figure(7)
    loglog(h_list,e_local_XB1,'ro','MarkerFaceColor','r','markersize',3, 'DisplayName', 'XB1');
    hold on
    loglog(h_list,e_local_XB2,'bo','MarkerFaceColor','b','markersize',3, 'DisplayName', 'XB2');
    loglog(h_list,XB_diff,'go','MarkerFaceColor','g','markersize',3, 'DisplayName', 'XBdiff');
    xlabel('step-size')
    ylabel('local truncation error')
    legend()
    title('Local Truncation Error vs Step Size')


%     figure(8)
%     subplot(2, 2, 1)
%     loglog(h_list, XB1_list(1, :), 'r--', 'DisplayName', 'XB1 x-pos')
%     hold on
%     loglog(h_list, XB2_list(1, :), 'b--', 'DisplayName', 'XB2 x-pos')
%     legend()
%     title("x-position plot")
% 
%     subplot(2, 2, 2)
%     semilogx(h_list, XB1_list(2, :), 'DisplayName', 'XB1 y-pos')
%     hold on
%     semilogx(h_list, XB2_list(2, :), 'DisplayName', 'XB2 y-pos')
%     legend()
%     title("y-position plot")
% 
%     subplot(2, 2, 3)
%     semilogx(h_list, XB1_list(3, :), 'r--', 'DisplayName', 'XB1 x-velocity')
%     hold on
%     semilogx(h_list, XB2_list(3, :), 'b--', 'DisplayName', 'XB2 x-velocity')
%     legend()
%     title("x-velocity plot")
% 
%     subplot(2, 2, 4)
%     semilogx(h_list, XB1_list(4, :), 'DisplayName', 'XB1 y-velocity')
%     hold on
%     semilogx(h_list, XB2_list(4, :), 'DisplayName', 'XB2 y-velocity')
%     legend()
%     title("y-velocity plot")
% 
%     XB_diff = vecnorm(XB1_list-XB2_list);
% 
%     figure(9)
%     loglog(h_list, XB_diff)
% 

% ----------------------- ADAPTIVE RK STEP ------------------------------

p = 5;  % Dormand–Prince method
error_desired = 1e-5;
h_start = 0.1;

[~, X_var, ~, ~, ~] = explicit_RK_variable_step_integration(gravity_rate_func_wrapper, [0, 300], X0, h_ref, DormandPrince, p, error_desired);


% compare orbit to real solution
figure(10)
plot(0,0,'yo','markerfacecolor','y','markersize',8) % sun
hold on
plot(X_var(:,1), X_var(:,2), 'r', 'DisplayName', 'Adaptive RK')
plot(V_list(:,1), V_list(:,2), 'k--', 'DisplayName', 'Exact orbit')
axis equal; axis([-20, 20, -20, 20])
title('Adaptive RK Orbit vs Exact Solution')
xlabel('Position (x)')
ylabel('Position (y)')
legend

% --------------------- New define planet start --------------------------
    x0 = 3;
    y0 = 5;
    dxdt0 = 3.5;
    dydt0 = 0;

    ti = 0;             % starting time
    tf = 300;            % ending time

    tspan = [ti, tf];

    % Define X0 as starting point location & velocity of planet
    X0 = [x0;y0;dxdt0;dydt0];
    solution_func2 = @(t_in) compute_planetary_motion(t_in, X0, orbit_params);
    t_exact = linspace(tspan(1),tspan(2),300);
    V_list = compute_planetary_motion(t_exact,X0,orbit_params);
    
    h_ref= 0.01;

    [t_list,X_list,h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration(gravity_rate_func_wrapper, tspan, X0, h_ref, HeunEuler, p, error_desired);


    figure(11)
    axis([-20,40,-50,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5, 'DisplayName', 'Sun');
    plot(X_list(:,1),X_list(:,2),'r','DisplayName', 'Calculated'); 
    plot(V_list(:,1),V_list(:,2),'k--', 'DisplayName', 'Exact'); 
    legend()
    title('Calculated vs Exact Position')
    xlabel('Position (x)')
    ylabel('Position (y)')

   % ---------- ADAPTIVE VS FIXED STEP GLOBAL TRUNCATION ERROR -----------

    % Adaptive Step
    h_ref= 0.001;
    p = 5;          % DormandPrince
    error_desired_list = logspace(-12, -7, 50);
    global_trunc_error_adaptive = zeros(size(error_desired_list));
    avg_step_size_adaptive = zeros(size(error_desired_list));
    num_evals_adaptive = zeros(size(error_desired_list));
    step_failure_rate_list = zeros(size(error_desired_list));
    for i = 1:length(error_desired_list)
        [t_list, X_list, h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration(gravity_rate_func_wrapper, tspan, X0, h_ref, DormandPrince, p, error_desired_list(i));
        global_trunc_error_adaptive(i) = norm(X_list(end, :)'-solution_func2(t_list(end)));
        avg_step_size_adaptive(i) = h_avg;
        num_evals_adaptive(i) = num_evals;
        step_failure_rate_list(i) = step_failure_rate;
    end


    % Fixed Step
    step_size_list = logspace(-3, 0, 50);
    % step_size_list = avg_step_size_sweeperror;
    global_trunc_error_fixed = zeros(size(step_size_list));
    avg_step_size_fixed = zeros(size(step_size_list));
    num_evals_fixed = zeros(size(step_size_list));

    for i = 1:length(step_size_list)
        [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(gravity_rate_func_wrapper,tspan,X0,step_size_list(i),DormandPrince);
        global_trunc_error_fixed(i) = norm(X_list(end, :)'-solution_func2(t_list(end)));
        avg_step_size_fixed(i) = h_avg;
        num_evals_fixed(i) = num_evals;

    end

    figure(12)
    loglog(avg_step_size_adaptive, global_trunc_error_adaptive, 'DisplayName', 'Adaptive')
    hold on
    loglog(avg_step_size_fixed, global_trunc_error_fixed, 'DisplayName', 'Fixed Step')
    legend()
    xlabel("Avg Step Size")
    ylabel("Global Truncation Error")
    title("Step Size vs Global Truncation Error")

    figure(13)
    loglog(num_evals_adaptive, global_trunc_error_adaptive, 'DisplayName', 'Adaptive')
    hold on
    loglog(num_evals_fixed, global_trunc_error_fixed, 'DisplayName', 'Fixed Step')
    legend()
    xlabel("Number of Evaluations")
    ylabel("Global Truncation Error")
    title("# Evals vs Global Truncation Error")

    figure(14)
    semilogx(avg_step_size_adaptive, step_failure_rate_list, 'ko')
    xlabel("Average Step Size")
    ylabel("Step Failure Rate")
    title("Step Size vs. Failure Rate")


    % ------------ ANALYZING ADAPTIVE STEP FUNCTION ---------------------

    desired_error = 1e-6;
    [t_list, X_list, h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration(gravity_rate_func_wrapper, tspan, X0, h_ref, DormandPrince, p, desired_error);
    
    position_norms = vecnorm(X_list(:, 1:2)');
    % Position Plot
    figure(15)
    plot(X_list(:, 1), X_list(:, 2),'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2)
    hold on
    xlabel('Position (x)')
    ylabel('Position (Y)')
    title('Position vs Time (Adaptive)')
    plot(0,0,'yo','markerfacecolor','y','markersize',8) % sun



    velocity_norms = vecnorm(X_list(:, 3:4)');
    % Velocity Plot
    figure(16)
    plot(X_list(:, 3), X_list(:, 4),'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2)
    xlabel('Velocity (x)')
    ylabel('Velocity (y)')
    title('Velocity vs Time (Adaptive)')

    h_list = diff(t_list);
    figure(17)
    semilogx(h_list, position_norms(2:end),'bo')
    xlabel('Step Size')
    ylabel('Distance Planet to Sun')
    title('Distance vs Step Size')

end
