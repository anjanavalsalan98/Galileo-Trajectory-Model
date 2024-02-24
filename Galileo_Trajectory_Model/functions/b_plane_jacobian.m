function [bplane_partials,J] = b_plane_jacobian(mu,r,v)
    [bplane, ~,~,~] = b_plane(mu, r, v);
    BT = bplane(8);
    BR = bplane(9);

    delta_Vx = logspace(-15, 2, 1000);
    delta_Vy = logspace(-15, 2, 1000);
    delta_Vz = logspace(-15, 2, 1000);
    V_index = length(delta_Vx)/2;
   
    del_BT_Vx_arr = zeros(length(delta_Vx), 1);
    for i = 1:length(delta_Vx)
        delta_V = [v(1)+delta_Vx(i) v(2) v(3)];
        [bplane_plus, ~,~,~] = b_plane(mu, r, delta_V);
        BT_plus = bplane_plus(8);
        del_BT_Vx_arr(i) = (BT_plus - BT)/delta_Vx(i);
    end
    del_BT_Vx = del_BT_Vx_arr(V_index);

    del_BT_Vy_arr = zeros(length(delta_Vy), 1);
    for i = 1:length(delta_Vy)
        delta_V = [v(1) v(2)+delta_Vy(i) v(3)];
        [bplane_plus, ~,~,~] = b_plane(mu, r, delta_V);
        BT_plus = bplane_plus(8);
        del_BT_Vy_arr(i) = (BT_plus - BT)/delta_Vy(i);
    end
    del_BT_Vy = del_BT_Vy_arr(V_index);

    del_BT_Vz_arr = zeros(length(delta_Vz), 1);
    for i = 1:length(delta_Vz)
        delta_V = [v(1) v(2) v(3)+delta_Vz(i)];
        [bplane_plus,~,~, ~] = b_plane(mu, r, delta_V);
        BT_plus = bplane_plus(8);
        del_BT_Vz_arr(i) = (BT_plus - BT)/delta_Vz(i);
    end
    del_BT_Vz = del_BT_Vz_arr(V_index);
    
    del_BR_Vx_arr = zeros(length(delta_Vx), 1);
    for i = 1:length(delta_Vx)
        delta_V = [v(1)+delta_Vx(i) v(2) v(3)];
        [bplane_plus,~,~,~] = b_plane(mu, r, delta_V);
        BR_plus = bplane_plus(9);
        del_BR_Vx_arr(i) = (BR_plus - BR)/delta_Vx(i);
    end
    del_BR_Vx = del_BR_Vx_arr(V_index);

    del_BR_Vy_arr = zeros(length(delta_Vy), 1);
    for i = 1:length(delta_Vy)
        delta_V = [v(1) v(2)+delta_Vy(i) v(3)];
        [bplane_plus,~,~,~] = b_plane(mu, r, delta_V);
        BR_plus = bplane_plus(9);
        del_BR_Vy_arr(i) = (BR_plus - BR)/delta_Vy(i);
    end
    del_BR_Vy = del_BR_Vy_arr(V_index);

    del_BR_Vz_arr = zeros(length(delta_Vz), 1);
    for i = 1:length(delta_Vz)
        delta_V = [v(1) v(2) v(3)+delta_Vz(i)];
        [bplane_plus,~,~,~] = b_plane(mu, r, delta_V);
        BR_plus = bplane_plus(9);
        del_BR_Vz_arr(i) = (BR_plus - BR)/delta_Vz(i);
    end
    del_BR_Vz = del_BR_Vz_arr(V_index);

    bplane_partials(1,:) = delta_Vx; %Pertubations in velocity x-component
    bplane_partials(2,:) = delta_Vy; %Pertubations in velocity y-component
    bplane_partials(3,:) = delta_Vz; %Pertubations in velocity z-component
    bplane_partials(4,:) = del_BT_Vx_arr; %Change in BT wrt Vx pertubations
    bplane_partials(5,:) = del_BT_Vy_arr; %Change in BT wrt Vy pertubations
    bplane_partials(6,:) = del_BT_Vz_arr; %Change in BT wrt Vz pertubations
    bplane_partials(7,:) = del_BR_Vx_arr; %Change in BR wrt Vx pertubations
    bplane_partials(8,:) = del_BR_Vy_arr; %Change in BR wrt Vy pertubations
    bplane_partials(9,:) = del_BR_Vz_arr; %Change in BR wrt Vz pertubations
    
    J = [del_BT_Vx, del_BT_Vy, del_BT_Vz;
         del_BR_Vx, del_BR_Vy, del_BR_Vz]; %Jacobian matrix 


end


