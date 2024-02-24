
% ***************************************************************************************
%    Title: rv2bplane.m
%    Author: David Eagle
%    Date: 20 Sep 2019
%    Code version: 2.0.0
%    Availability: https://www.mathworks.com/matlabcentral/fileexchange/39462-gravity-assist-trajectory-design-and-analysis
% ***************************************************************************************/

function [bplane, S,R, T] = b_plane(mu,r, v)
    flight_path_angle = asin(dot(r,v)/norm(dot(r,v)));

    h = cross(r,v); % Normalize by dividing by the cross
    h_mag = norm(h);

    R_mag = norm(r);
    V_mag = norm(v);

    c3 = V_mag^2 - 2*mu/R_mag;
    V_hyp_mag = sqrt(c3); %Hyperbolic speed

    p = (h_mag^2)/mu; %semi-latus rectum
    a =  R_mag / (2 - R_mag * V_mag^2 / mu); %semi-major axis
    ecc = sqrt(1 - p / a); %orbital eccentricity

    radius_rate = dot(r,v)/R_mag;

    TA_cos = (p - R_mag) / (ecc * R_mag);
    TA_sin = radius_rate * h_mag / (ecc * mu);


    b = sqrt(p* abs(a));
    c = sqrt(a^2 + b^2);

    z = R_mag / h_mag * v - radius_rate / h_mag * r;

    P = TA_cos * r / R_mag - TA_sin * z;
    Q = TA_sin * r / R_mag + TA_cos * z;

    S = -a / c * P + b / c * Q;
    S = S/norm(S);
    B =  b^2 / c * P + a * b / c * Q;

    decl = asin(S(3)); %declination of asymptote
    RA = atan3(S(2), S(1)); %right ascension of asymptote
    r_per = (ecc - 1) * mu / c3;
    
    T(1) = S(2)/sqrt(S(1)^2 + S(2)^2);
    T(2) = -S(1)/sqrt(S(1)^2 + S(2)^2);
    T(3) = 0;

    R = cross(S,T);
    R = R/norm(R);

    BT = dot(B,T);
    BR = dot(B,R);

    theta = atan3(BR, BT);

    B_mag = r_per * sqrt(1 + 2 * mu / (r_per * c3));
   
    bplane(1) = V_hyp_mag;  %Hyperbolic speed
    bplane(2) = flight_path_angle; %Flight path angle
    bplane(3) = decl; %declination of asymptote
    bplane(4) = RA; %right ascension of asymptote
    bplane(5) = R_mag; %Magnitude of r
    bplane(6) = r_per; %periapsis radius of hyperbola
    bplane(7) = theta; %b-plane angle
    bplane(8) = BT; %B dot T
    bplane(9) = BR; %B dot R
    bplane(10) = B_mag; %b magnitude
    bplane(11:13) = B; %b-plane vector

end