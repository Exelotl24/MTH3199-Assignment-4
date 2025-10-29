function [E_list, H_list] = conservation(t_list, V_list, orbit_params)
    E_list = .5*(orbit_params.m_planet)*(V_list(3)^2 + V_list(4)^2) - ... 
    (orbit_params.m_sun * orbit_params.m_planet * orbit_params.G)/(norm(V_list(:,1:2)));
end