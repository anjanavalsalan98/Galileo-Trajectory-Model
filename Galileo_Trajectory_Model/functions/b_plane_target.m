
% ***************************************************************************************
%    Title: rv2bplane.m
%    Author: David Eagle
%    Date: 20 Sep 2019
%    Code version: 2.0.0
%    Availability: https://www.mathworks.com/matlabcentral/fileexchange/39462-gravity-assist-trajectory-design-and-analysis

function [BT,BR, S_new, R_new, T_new] = b_plane_target(mu, r, v, S, inc, R_per)
    V_inf_mag = sqrt(norm(v)^2 - 2*mu/norm(r));
    fpa = 0; %Set flight path angle to 0 to target periapse radius on B-plane
    decl = pi/2 - acos(v(3)/norm(v));
    RA = atan3(v(2),v(1));

    b_target = cos(fpa) * sqrt((2*mu*R_per/V_inf_mag^2) + R_per^2);

    cos_theta = cos(inc)/cos(decl);
    sin_theta = -sqrt(1-cos_theta^2);

    BT = b_target * cos_theta;
    BR = b_target * sin_theta;

    sin_decl = sqrt(S(1)^2 + S(2)^2);
    cos_decl = sqrt(1-sin(decl)^2);
    S_new = [cos_decl * cos(RA), cos_decl * sin(RA), sin_decl];
    
    T_new(1) = S_new(2)/sqrt(S_new(1)^2 + S_new(2)^2);
    T_new(2) = -S_new(1)/sqrt(S_new(1)^2 + S_new(2)^2);
    T_new(3) = 0;

    R_new = cross(S_new,T_new);
    R_new = R_new/norm(R_new);
end

