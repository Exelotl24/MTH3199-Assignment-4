main()
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
    tf = 15;            % ending time

    tspan = [ti, tf];

    % Define X0 as starting point location & velocity of planet
    X0 = [x0;y0;dxdt0;dydt0];
    t_exact = linspace(tspan(1),tspan(2),100);
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

    t_ref = 0.5;
    [h_list, analytical_difference, e_local_FE, e_local_MP, e_local_BE, e_local_IMP] = local_truncation_error(gravity_rate_func_wrapper, solution_func, t_ref)


    figure()
    loglog(h_list,analytical_difference, 'DisplayName', 'Analytical Difference');
    hold on
    loglog(h_list, e_local_FE, 'DisplayName', 'Forward Euler Error')
    loglog(h_list, e_local_MP, 'DisplayName', 'Explicit Midpoint Error')
    % loglog(h_list, k_FE*h_list.^p_FE_01, 'DisplayName', 'FE linear')
    % loglog(h_list, k_MP_01*h_list.^p_MP_01, 'DisplayName', 'MP linear')
    title("Local Truncation Error")
    legend


     % -------------- CONSERVATION OF PHYSICAL QUALITIES ----------------

     [E_list, H_list] = conservation(t_list, V_list, orbit_params);

end

