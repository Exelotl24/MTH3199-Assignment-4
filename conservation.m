function [E_effectivity, H_effectivity] = conservation(t_list, V_list, orbit_params)
    E_list = .5*(orbit_params.m_planet)*(V_list(:, 3).^2 + V_list(:, 4).^2) - ... 
    (orbit_params.m_sun * orbit_params.m_planet * orbit_params.G)./(vecnorm(V_list(:,1:2), 2, 2));

    H_list = orbit_params.m_planet*(V_list(:, 1).*V_list(:, 4) - V_list(:, 2).*V_list(:, 3));
    
    E_error_list = E_list-E_list(1);
    H_error_list = H_list-H_list(1);

    figure()
    yyaxis left
    plot(t_list, E_error_list, '--', 'DisplayName', 'Energy Error')
    ylabel('Energy Error (J)')
    yyaxis right
    plot(t_list, H_error_list, 'DisplayName', 'Angular Momentum Error')
    ylabel('Angular Momentum Error (kg m2/s)')
    legend()
    xlabel('Time')
    title('Conservation of Energy')

    E_effectivity = norm(E_error_list);
    H_effectivity = norm(H_error_list);
end