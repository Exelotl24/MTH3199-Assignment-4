% clear
% close all
% 
% main()

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

    ti = 0;             % starting time
    tf = 300;            % ending time

    tspan = [ti, tf];

    % Define X0 as starting point location & velocity of planet
    X0 = [x0;y0;dxdt0;dydt0];
    t_exact = linspace(tspan(1),tspan(2),300);
    V_list = compute_planetary_motion(t_exact,X0,orbit_params);
    
    
    
    h_ref= 0.1;
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(gravity_rate_func_wrapper,tspan,X0,h_ref,BT_struct_EM)

    figure(1);
    subplot(2,1,1)
    hold on
    plot(t_list,X_list(:,1),'r');
    plot(t_list,X_list(:,2),'b');
    plot(t_exact,V_list(:,1),'k--');
    plot(t_exact,V_list(:,2),'k--');

    xlabel('time')
    ylabel('position')

    subplot(2,1,2)
    hold on
    plot(t_list,X_list(:,3),'r');
    plot(t_list,X_list(:,4),'b');
    plot(t_exact,V_list(:,3),'k--');
    plot(t_exact,V_list(:,4),'k--');

    xlabel('time')
    ylabel('velocity')

    figure(2)
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(X_list(:,1),X_list(:,2),'r');
    plot(V_list(:,1),V_list(:,2),'k--');


    % ------------------- LOCAL TRUNCATION ---------------------
    V0 = [x0;y0;dxdt0;dydt0];
    solution_func = @(t_in) compute_planetary_motion(t_in, V0, orbit_params);

    HeunEuler = struct();
    HeunEuler.C = [0,1];
    HeunEuler.B = [1/2,1/2;1,0];
    HeunEuler.A = [0,0;1,0];


    t_ref = 0.5;
    [h_list, analytical_difference, e_local_FE, e_local_MP, e_local_Heun] = local_truncation_error(gravity_rate_func_wrapper, solution_func, t_ref, BT_struct_FE, BT_struct_EM, HeunEuler);

    figure()
    loglog(h_list,analytical_difference, 'DisplayName', 'Analytical Difference');
    hold on
    loglog(h_list, e_local_FE, 'DisplayName', 'Forward Euler Error')
    loglog(h_list, e_local_MP, 'DisplayName', 'Explicit Midpoint Error')
    loglog(h_list, e_local_Heun, 'DisplayName', 'Heun Euler')
    title("Local Truncation Error")
    legend


     % -------------- CONSERVATION OF PHYSICAL QUALITIES ----------------

     [EV, HV] = conservation(t_exact, V_list, orbit_params);
     [EX, HX] = conservation(t_list, X_list, orbit_params);


     % -------------- LOCAL TRUNCATION -----------------------

     h_list = logspace(-4, -1, 100);


    ti = 0;             % starting time
    tf = 300;           % ending time
    
    tspan = [ti, tf];
    
    % Define X0 as starting point location & velocity of planet
    X0 = [x0;y0;dxdt0;dydt0];

    
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

    figure()
    loglog(h_list,e_local_XB1,'ro','MarkerFaceColor','r','markersize',3);
    hold on
    loglog(h_list,e_local_XB2,'bo','MarkerFaceColor','b','markersize',3);
    loglog(h_list,XB_diff,'go','MarkerFaceColor','g','markersize',3);
%     figure()
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
%     figure()
%     loglog(h_list, XB_diff)
% 

% ----------------------- ADAPTIVE RK STEP ------------------------------

p = 5;  % Dormand–Prince method
error_desired = 1e-5;
h_start = 0.1;

[XB, num_evals_step, h_next, redo] = explicit_RK_variable_step(...
    gravity_rate_func_wrapper, 0, X0, h_start, DormandPrince, p, error_desired);

% disp('Single step test:');
% disp(['redo = ', num2str(redo)]);
% disp(['New step size h_next = ', num2str(h_next)]);
% disp(['Number of function evaluations = ', num2str(num_evals_step)]);
% disp('New state estimate XB = ');
disp(XB);


[t_var, X_var, h_avg_var, num_evals_var] = explicit_RK_variable_step_integration(...
    gravity_rate_func_wrapper, [0, 300], X0, h_ref, DormandPrince, p, error_desired);


% compare orbit to real solution
figure();
plot(0,0,'yo','markerfacecolor','y','markersize',8); % sun
hold on;
plot(X_var(:,1), X_var(:,2), 'r', 'DisplayName', 'Adaptive RK');
plot(V_list(:,1), V_list(:,2), 'k--', 'DisplayName', 'Exact orbit');
axis equal; axis([-20, 20, -20, 20]);
title('Adaptive RK Orbit vs Exact Solution');
xlabel('x position'); ylabel('y position');
legend show;


end

