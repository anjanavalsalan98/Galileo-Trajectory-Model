%                    ********************************
%                    *  ____       _ _ _            *
%                    * / ___| __ _| (_) | ___  ___  *
%                    *| |  _ / _` | | | |/ _ \/ _ \ *
%                    *| |_| | (_| | | | |  __/ (_) |*
%                    * \____|\__,_|_|_|_|\___|\___/ *
%                    ********************************

%~~~~~~~~~~~~~~~~~~~~~~~~~ Environment Set Up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear;
clc;
close all;
disp('TRAJECTORY MODEL OF GALILEO INTERPLANETARY CRUISE');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MICE Kernels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
directory = pwd;
path = strcat(directory,'/mice/src/mice/');
path_lib = strcat(directory,'/mice/lib/');
path_fn = strcat(directory,'/functions/');
addpath(path);
addpath(path_lib);
addpath(path_fn);

cspice_furnsh( { strcat(directory,'/kernels/spk/s030916a.bsp')});
cspice_furnsh( { strcat(directory,'/kernels/spk/s960730a.bsp')});
cspice_furnsh( { strcat(directory,'/kernels/spk/s970311a.bsp')});
cspice_furnsh( { strcat(directory,'/kernels/spk/s980326a.bsp')});
cspice_furnsh( { strcat(directory,'/kernels/spk/s000131a.bsp')});
cspice_furnsh( { strcat(directory,'/kernels/lsk/naif0008.tls')});
cspice_furnsh( { strcat(directory,'/kernels/lsk/mk98264a.tls')});
cspice_furnsh( { strcat(directory,'/kernels/pck/pck00007.tpc')});
cspice_furnsh( { strcat(directory,'/kernels/pck/pk96030a.tpc')});
cspice_furnsh( { strcat(directory,'/kernels/pck/mips_010314.tpc')});
cspice_furnsh( { strcat(directory,'/kernels/sclk/mk00062a.tsc')});
cspice_furnsh( { strcat(directory,'/kernels/spk/de432s.bsp')});
cspice_furnsh( { strcat(directory,'/kernels/spk/20000951.bsp')});
cspice_furnsh( { strcat(directory,'/kernels/pck/gm_de431.tpc')});

%~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initial Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sun_ID = '10';
earth_ID = '399';
gaspra_ID = '20000951';
venus_ID = '299';
ida_ID = '2431010';
jupiter_ID = '599';
gal_ID = '-77';
europa_ID = '502';
io_ID = '501';

mass_gaspra = 3e15;
mass_ida = 100e15;
UC = 6.67430e-11;

mu_sun = cspice_bodvrd( 'SUN', 'GM', 1 );
mu_earth = cspice_bodvrd( 'EARTH', 'GM', 1 );
mu_venus = cspice_bodvrd( 'VENUS', 'GM', 1 );
mu_jupiter = cspice_bodvrd( 'JUPITER', 'GM', 1 );
mu_gaspra = mass_gaspra * UC;
mu_ida = mass_ida * UC;

rad_venus = mean(cspice_bodvrd('VENUS','RADII',3));
rad_earth = mean(cspice_bodvrd('EARTH','RADII',3));
rad_gaspra = mean(cspice_bodvrd('GASPRA','RADII',3));
rad_ida = mean(cspice_bodvrd('IDA','RADII',3));
rad_jupiter = mean(cspice_bodvrd('JUPITER','RADII',3));

cruise_table_sz = [1 14];
cruise_table_types = ["int16", "string", repmat("double", 1, 12)];
cruise_table_headings = ["Phase", "Segment", ...
                         "vhel_planet_x", "vhel_planet_y", "vhel_planet_z", "Vhel_planet_mag", ...
                         "vhel_gal_x", "vhel_gal_y", "vhel_gal_z", "Vhel_gal_mag", ...
                         "vinf_gal_x", "vinf_gal_y", "vinf_gal_z", "Vinf_gal_mag"];
cruise_table = table('Size', cruise_table_sz, 'VariableTypes', cruise_table_types, 'VariableNames', cruise_table_headings);

flyby_table_sz = [1 10];
flyby_table_types = ["int16", "string", repmat("double", 1, 4), "string", repmat("double", 1, 3)];
flyby_table_headings = ["Phase", "Segment", ...
                        "tcm_vx", "tcm_vy", "tcm_vz", "tcm_Vmag", ...
                        "UTC_ca", "alt_ca", "inc_ca", ...
                        "flyby_deltaV"];
flyby_table = table('Size', flyby_table_sz, 'VariableTypes', flyby_table_types, 'VariableNames', flyby_table_headings);

%**************************************************************************
%                                  PHASE 1: 
%                               EARTH TO VENUS
%**************************************************************************
disp('-> Phase 1: Earth to Venus ✓');

%--------------------------------------------------------------------------
% 1.1 Earth to Venus Cruise
%--------------------------------------------------------------------------

%Getting states of Earth at Departure and Venus at Arrival
utc_earth_dep = '1989-10-19T01:40:00';
utc_venus_arr = '1990-2-10T05:59:00';
[earth_dep_state, et_earth_dep] = spice_states(utc_earth_dep, utc_earth_dep, 1, sun_ID, earth_ID);
[venus_arr_state, et_venus_arr] = spice_states(utc_venus_arr, utc_venus_arr, 1, sun_ID, venus_ID);

%Calculating time of flight between Earth departure and Venus arrival
JD_earth_dep = julian_date(utc_earth_dep);
JD_venus_arr = julian_date(utc_venus_arr);
tof_e2v = (JD_venus_arr - JD_earth_dep)*(24*60*60);

%Solving Lambert's Problem to find velocity vectors at Earth departure and Venus arrival
num_rev_e2v = 0;
[v_hel_Edep, v_hel_Varr] = glambert(mu_sun, earth_dep_state, venus_arr_state, tof_e2v, num_rev_e2v);

%Orbit Propagation using Two Body Problem
[e2v_states, e2v_etvec] = twobody(mu_sun, et_earth_dep, et_venus_arr, [earth_dep_state(1:3) v_hel_Edep]);
e2v_r = e2v_states(:,1:3);
e2v_v = e2v_states(:,4:6);

%--------------------------------------------------------------------------
% 1.2 Earth Departure
%--------------------------------------------------------------------------

%Calulating the launch heliocentric delta-v (at Earth departure)
delta_v_Edep = v_hel_Edep - earth_dep_state(4:6);
delta_V_Edep_mag = norm(v_hel_Edep) - norm(earth_dep_state(4:6));

%Calculating Injection Energy
%Square of asymptotic speed wrt Earth on departure hyperbola
C3_Edep = delta_V_Edep_mag^2;

%Adding data to cruise table
row_c = 1;
cruise_table(row_c,:) = {1, "Edep", ...
                       earth_dep_state(4), earth_dep_state(5), earth_dep_state(6), norm(earth_dep_state(4:6)), ...
                       v_hel_Edep(1), v_hel_Edep(2), v_hel_Edep(3), norm(v_hel_Edep), ...
                       delta_v_Edep(1), delta_v_Edep(2), delta_v_Edep(3), norm(delta_v_Edep)};

%--------------------------------------------------------------------------
% 1.3 Venus Arrival
%--------------------------------------------------------------------------

%Calculating incoming v-infinity vector
v_inf_Varr = v_hel_Varr - venus_arr_state(4:6);
V_inf_Varr_mag = norm(v_inf_Varr);

%Finding the state vector of spacecraft at Venus SOI
for i = 1:length(e2v_etvec)
    utc_venus = cspice_et2utc(e2v_etvec(i),'ISOC', 3);
    [state_venus, et_venus] = spice_states(utc_venus, utc_venus, 1, sun_ID, venus_ID);
    venus_SOI = norm(state_venus(1:3)) * (mu_venus/mu_sun)^(2/5);

    position_gal = e2v_r(i,:) - state_venus(1:3);
    velocity_gal = e2v_v(i,:) - state_venus(4:6);

    dist_venus_gal = norm(position_gal);

    if dist_venus_gal <= venus_SOI
        break;
    end
end

if i == length(e2v_etvec)
    error('Galileo does not intercept Venus SOI');
else
    gal_VSOI = [position_gal velocity_gal];
    et_VSOI_gal = e2v_etvec(i);
end

%Add data to cruise table
row_c = row_c + 1;
cruise_table(row_c,:) = {1, "Varr", ...
                       venus_arr_state(4), venus_arr_state(5), venus_arr_state(6), norm(venus_arr_state(4:6)), ...
                       v_hel_Varr(1), v_hel_Varr(2), v_hel_Varr(3), norm(v_hel_Varr), ...
                       v_inf_Varr(1), v_inf_Varr(2), v_inf_Varr(3), V_inf_Varr_mag};

%--------------------------------------------------------------------------
% 1.4 B-Plane Targeting of Venus Flyby 
%--------------------------------------------------------------------------

%Getting current b-plane parameters
[bplane_Vfb,S, R, T] = b_plane(mu_venus,gal_VSOI(1:3), gal_VSOI(4:6));

%Getting target b-plane
Vfb_inc = deg2rad(138.1);
Vfb_alt = 16123;
Vfb_R_per = Vfb_alt + rad_venus;
[BT_star,BR_star, ~, ~, ~] = b_plane_target(mu_venus, gal_VSOI(1:3), gal_VSOI(4:6), S, Vfb_inc, Vfb_R_per);

%Calculating b-plane Jacobian and its pseudoinverse
[bplane_Vfb_partials, J] = b_plane_jacobian(mu_venus, gal_VSOI(1:3), gal_VSOI(4:6));
J_inv = transpose(J) * (J * transpose(J))^(-1);

%Calculating delta B
BT = bplane_Vfb(8);
BR = bplane_Vfb(9);
delta_B_venus = [BT_star - BT;
                 BR_star - BR];

%Calculating delta V for TCM-1
tcm1_Vfb = transpose(J_inv * delta_B_venus);

%Getting updated state of Galileo at Venus SOI after TCM-1
gal_VSOI_tcm1 = [gal_VSOI(1:3) gal_VSOI(4:6)+tcm1_Vfb];

%--------------------------------------------------------------------------
% 1.5 Closest Approach Parameters of Venus FlyBy 
%--------------------------------------------------------------------------

%Specifying start and stop propagation times in Venus SOI
utc_Vfb_stop = '1990-2-10T09:59:00';
et_Vfb_stop = cspice_str2et(utc_Vfb_stop);
et_Vfb_start = et_VSOI_gal;

%Propagating Galileo trajectory in Venus SOI without TCM-1
[gal_Vfb_states, gal_Vfb_etvec] = twobody(mu_venus, et_Vfb_start, et_Vfb_stop, gal_VSOI);
gal_Vfb_r = gal_Vfb_states(:,1:3);
gal_Vfb_v = gal_Vfb_states(:,4:6);

[Vca_inc, Vca_alt, Vca_et] = closest_approach(gal_Vfb_etvec, gal_Vfb_r, gal_Vfb_v, rad_venus);
Vca_utc = cspice_et2utc(Vca_et, 'ISOC', 3);
Vca_inc_deg = rad2deg(Vca_inc);

%Propagating Galileo trajectory in Venus SOI with TCM-1
[gal_Vfb_tcm1_states, gal_Vfb_tcm1_etvec] = twobody(mu_venus, et_Vfb_start, et_Vfb_stop, gal_VSOI_tcm1);
gal_Vfb_tcm1_r = gal_Vfb_tcm1_states(:,1:3);
gal_Vfb_tcm1_v = gal_Vfb_tcm1_states(:,4:6);

[Vca_tcm1_inc,Vca_tcm1_alt, Vca_tcm1_et] = closest_approach(gal_Vfb_tcm1_etvec, gal_Vfb_tcm1_r, gal_Vfb_tcm1_v ,rad_venus);
Vca_tcm1_utc = cspice_et2utc(Vca_tcm1_et, 'ISOC', 3);
Vca_tcm1_inc_deg = rad2deg(Vca_tcm1_inc);

%Add data to flyby table
row_fb = 1;
flyby_table(row_fb,:) = {1, "Vfb", ...
                       tcm1_Vfb(1), tcm1_Vfb(2), tcm1_Vfb(3), norm(tcm1_Vfb), ...
                       Vca_tcm1_utc, Vca_tcm1_alt, Vca_tcm1_inc_deg, 0};

%--------------------------------------------------------------------------
% 1.6 Plotting Graphs for Phase 1: Earth to Venus
%--------------------------------------------------------------------------

%Plotting Earth to Venus Cruise (SPICE vs Model)
[gal_spice_e2v_states, gal_spice_e2v_etvec] = spice_states(utc_earth_dep, utc_venus_arr, 2500, sun_ID , gal_ID);
gal_spice_e2v_r = gal_spice_e2v_states(:,1:3);
gal_spice_e2v_v = gal_spice_e2v_states(:,4:6);

figure(1);
marker = 20;
plot3(e2v_r(:,1), e2v_r(:,2), e2v_r(:,3),'-b');
hold on;
plot3(gal_spice_e2v_r(:,1), gal_spice_e2v_r(:,2), gal_spice_e2v_r(:,3),'--r');
plot3(0, 0, 0, '.','Color', '[1, 0.8, 0.1]','MarkerSize',marker*3);
plot3(earth_dep_state(1), earth_dep_state(2), earth_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(venus_arr_state(1), venus_arr_state(2), venus_arr_state(3), '.','Color', '[0.7, 0, 1]','MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model','SPICE', 'Sun', 'Earth', 'Venus');
title('Earth to Venus Cruise (SPICE vs Model)')

%Plotting Partial Derivatives wrt Velocity Perturbations
delta_Vx = bplane_Vfb_partials(1,:); 
delta_Vy = bplane_Vfb_partials(2,:); 
delta_Vz = bplane_Vfb_partials(3,:);
del_BT_Vx_arr = bplane_Vfb_partials(4,:);
del_BT_Vy_arr = bplane_Vfb_partials(5,:); 
del_BT_Vz_arr = bplane_Vfb_partials(6,:);
del_BR_Vx_arr = bplane_Vfb_partials(7,:); 
del_BR_Vy_arr = bplane_Vfb_partials(8,:);
del_BR_Vz_arr = bplane_Vfb_partials(9,:);

figure(2);
t = tiledlayout(3,2);
title(t,'Venus B-Plane Partial Derivatives');
nexttile;
plot(log10(delta_Vx), del_BT_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vx), del_BR_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BT_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BR_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BT_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_z}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BR_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_z}$','Interpreter','latex');

%Plotting b-plane with and without TCM
W = null(S);
[P,Q] = meshgrid(-2e5:1e5:2e5); 
X = 0+W(1,1)*P+W(1,2)*Q; 
Y = 0+W(2,1)*P+W(2,2)*Q; 
Z = 0+W(3,1)*P+W(3,2)*Q;

R_plot = [0 0 0; R.*1e5];
T_plot = [0 0 0; T.*1e5];
S_plot = [0 0 0; S.*1e5];

midpoint_S = S_plot(round(length(S_plot) / 2), :);
midpoint_T = T_plot(round(length(T_plot) / 2), :);
midpoint_R = R_plot(round(length(R_plot) / 2), :);

figure(3);
plot3(gal_Vfb_r(:,1), gal_Vfb_r(:,2), gal_Vfb_r(:,3),'-b');
hold on;
plot3(gal_Vfb_tcm1_r(:,1), gal_Vfb_tcm1_r(:,2), gal_Vfb_tcm1_r(:,3),'-r');
plot3(0,0,0,'.','Color', '[0.7, 0, 1]','MarkerSize', marker);
surf(X,Y,Z, 'EdgeColor', 'none', 'FaceColor', '[0.3, 1, 0.8]','FaceAlpha', 0.3);
plot3(S_plot(:,1),S_plot(:,2),S_plot(:,3),'-k');
plot3(T_plot(:,1),T_plot(:,2),T_plot(:,3),'-k');
plot3(R_plot(:,1),R_plot(:,2),R_plot(:,3),'-k');
text(midpoint_S(1), midpoint_S(2), midpoint_S(3), ' S', 'FontSize', 8, 'Color', 'black');
text(midpoint_T(1), midpoint_T(2), midpoint_T(3), ' T', 'FontSize', 8, 'Color', 'black');
text(midpoint_R(1), midpoint_R(2), midpoint_R(3), ' R', 'FontSize', 8, 'Color', 'black');
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Without TCM','With TCM','Venus','B-plane');
title('B-plane of Venus With and Without TCM');

%Plotting Galileo closest approach to Venus
utc_Vca_start = '1990-2-10T01:59:00';
et_Vca_start = cspice_str2et(utc_Vca_start);
utc_Vca_stop = cspice_et2utc(et_Vfb_stop, 'ISOC',1);

gal_Vca_tcm1_r = slice_flyby(gal_Vfb_tcm1_r, gal_Vfb_tcm1_etvec, et_Vca_start);

[gal_spice_Vca_states, gal_spice_Vca_etvec] = spice_states(utc_Vca_start, utc_Vca_stop, 2500, venus_ID , gal_ID);
gal_spice_Vca_r = gal_spice_Vca_states(:,1:3);
gal_spice_Vca_v = gal_spice_Vca_states(:,4:6);

[gal_spice_Vca_inc,gal_spice_Vca_alt,gal_spice_Vca_et] = closest_approach(gal_spice_Vca_etvec,gal_spice_Vca_r,gal_spice_Vca_v,rad_venus);
gal_spice_Vca_utc = cspice_et2utc(gal_spice_Vca_et, 'ISOC', 3);
gal_spice_Vca_inc_deg = rad2deg(gal_spice_Vca_inc);

figure(4)
plot3(gal_Vca_tcm1_r(:,1), gal_Vca_tcm1_r(:,2),gal_Vca_tcm1_r(:,3) ,'-b');
hold on;
plot3(gal_spice_Vca_r(:,1), gal_spice_Vca_r(:,2), gal_spice_Vca_r(:,3), '--r');
plot3(0,0,0,'.', 'Color', '[0.7, 0, 1]', 'MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model', 'SPICE', 'Venus');
title('Closest Approach to Venus - SPICE vs Model');

%**************************************************************************
%                                  PHASE 2: 
%                               VENUS TO EARTH-1
%**************************************************************************
disp('-> Phase 2: Venus to Earth-1 ✓');

%--------------------------------------------------------------------------
% 2.1 Venus to Earth-1 Cruise
%--------------------------------------------------------------------------

%Getting states of Venus at departure and Earth-1 at arrival
utc_venus_dep = '1990-2-10T05:59:00';
utc_earth1_arr = '1990-12-08T20:35:00';
[venus_dep_state, et_venus_dep] = spice_states(utc_venus_dep, utc_venus_dep, 1, sun_ID, venus_ID);
[earth1_arr_state, et_earth1_arr] = spice_states(utc_earth1_arr, utc_earth1_arr, 1, sun_ID, earth_ID);

%Calculating time of flight between Venus departure and Earth-1 arrival
JD_venus_dep = julian_date(utc_venus_dep);
JD_earth1_arr = julian_date(utc_earth1_arr);
tof_v2e1 = (JD_earth1_arr - JD_venus_dep)*(24*60*60);

%Solving Lambert's Problem to find velocity vectors at Venus departure and Earth-1 arrival
num_rev_v2e1 = 0;
[v_hel_Vdep, v_hel_E1arr] = glambert(mu_sun, venus_dep_state, earth1_arr_state, tof_v2e1, num_rev_v2e1);

%Orbit Propagation using Two Body Problem ignoring pertubations
[v2e1_states, v2e1_etvec] = twobody(mu_sun, et_venus_dep, et_earth1_arr, [venus_dep_state(1:3) v_hel_Vdep]);
v2e1_r = v2e1_states(:,1:3);
v2e1_v = v2e1_states(:,4:6);

%--------------------------------------------------------------------------
% 2.2 Venus Departure
%--------------------------------------------------------------------------

%Calulating heliocentric velocity change due to flyby
delta_v_Vdep = v_hel_Vdep -v_hel_Varr;
delta_v_Vdep_mag = norm(v_hel_Vdep) - norm(v_hel_Varr);

%Calculating outgoing v-infinity vector
v_inf_Vdep = v_hel_Vdep - venus_dep_state(4:6);
V_inf_Vdep_mag = norm(v_inf_Vdep);

%Add data to cruise table
row_c = row_c+1;
cruise_table(row_c,:) = {2, "Vdep", ...
                       venus_dep_state(4), venus_dep_state(5), venus_dep_state(6), norm(venus_dep_state(4:6)), ...
                       v_hel_Vdep(1), v_hel_Vdep(2), v_hel_Vdep(3), norm(v_hel_Vdep), ...
                       v_inf_Vdep(1), v_inf_Vdep(2), v_inf_Vdep(3), V_inf_Vdep_mag};

%Add data to flyby table
flyby_table.flyby_deltaV(row_fb) = delta_v_Vdep_mag;

%--------------------------------------------------------------------------
% 2.3 Earth-1 Arrival
%--------------------------------------------------------------------------

%Calculating incoming v-infinity vector
v_inf_E1arr = v_hel_E1arr - earth1_arr_state(4:6);
V_inf_E1arr_mag = norm(v_inf_E1arr);

%Finding the state vector of spacecraft at Earth-1 SOI
for i = 1:length(v2e1_etvec)
    utc_earth1 = cspice_et2utc(v2e1_etvec(i),'ISOC', 3);
    [state_earth1, et_earth1] = spice_states(utc_earth1, utc_earth1, 1, sun_ID, earth_ID);
    earth1_SOI = norm(state_earth1(1:3)) * (mu_earth/mu_sun)^(2/5);

    position_gal = v2e1_r(i,:) - state_earth1(1:3);
    velocity_gal = v2e1_v(i,:) - state_earth1(4:6);

    dist_earth1_gal = norm(position_gal);

    if dist_earth1_gal <= earth1_SOI
        break;
    end
end

if i == length(v2e1_etvec)
    error('Galileo does not intercept Earth-1 SOI');
else
    gal_E1SOI = [position_gal velocity_gal];
    et_E1SOI_gal = v2e1_etvec(i);
end

%Add data to cruise table
row_c = row_c + 1;
cruise_table(row_c,:) = {2, "E1arr", ...
                       earth1_arr_state(4), earth1_arr_state(5), earth1_arr_state(6), norm(earth1_arr_state(4:6)), ...
                       v_hel_E1arr(1), v_hel_E1arr(2), v_hel_E1arr(3), norm(v_hel_E1arr), ...
                       v_inf_E1arr(1), v_inf_E1arr(2), v_inf_E1arr(3), V_inf_E1arr_mag};

%--------------------------------------------------------------------------
% 2.4 B-Plane Targeting of Earth-1 Fly-by 
%--------------------------------------------------------------------------

%Getting current b-plane parameters
[bplane_E1fb,S, R, T] = b_plane(mu_earth,gal_E1SOI(1:3), gal_E1SOI(4:6));

%Getting target b-plane
E1fb_inc = deg2rad(142.9); %ecliptic frame
E1fb_alt = 960;
E1fb_R_per = E1fb_alt + rad_earth;
[BT_star,BR_star, ~, ~, ~] = b_plane_target(mu_earth, gal_E1SOI(1:3), gal_E1SOI(4:6), S, E1fb_inc, E1fb_R_per);

%Calculating b-plane Jacobian and its pseudoinverse
[bplane_E1fb_partials, J] = b_plane_jacobian(mu_earth, gal_E1SOI(1:3), gal_E1SOI(4:6));
J_inv = transpose(J) * (J * transpose(J))^(-1);

%Calculating delta B
BT = bplane_E1fb(8);
BR = bplane_E1fb(9);
delta_B_earth1 = [BT_star - BT;
                 BR_star - BR];

%Calculating delta V for TCM-1
tcm2_E1fb = transpose(J_inv * delta_B_earth1);

%Getting updated state of Galileo at Venus SOI after TCM-1
gal_E1SOI_tcm2 = [gal_E1SOI(1:3) gal_E1SOI(4:6)+tcm2_E1fb];

%--------------------------------------------------------------------------
% 2.5 Closest Approach Parameters of Earth-1 Fly-By 
%--------------------------------------------------------------------------

%Specifying start and stop propagation times in Earth-1 SOI
utc_E1fb_stop = '1990-12-09T01:35:00';
et_E1fb_stop = cspice_str2et(utc_E1fb_stop);
et_E1fb_start = et_E1SOI_gal;

%Propagating Galileo trajectory in Earth-1 SOI without TCM-2
[gal_E1fb_states, gal_E1fb_etvec] = twobody(mu_earth, et_E1fb_start, et_E1fb_stop, gal_E1SOI);
gal_E1fb_r = gal_E1fb_states(:,1:3);
gal_E1fb_v = gal_E1fb_states(:,4:6);

[E1ca_inc, E1ca_alt, E1ca_et] = closest_approach(gal_E1fb_etvec, gal_E1fb_r, gal_E1fb_v, rad_earth);
E1ca_utc = cspice_et2utc(E1ca_et, 'ISOC', 3);
E1ca_inc_deg = rad2deg(E1ca_inc);

%Propagating Galileo trajectory in Venus SOI with TCM-1
[gal_E1fb_tcm2_states, gal_E1fb_tcm2_etvec] = twobody(mu_earth, et_E1fb_start, et_E1fb_stop, gal_E1SOI_tcm2);
gal_E1fb_tcm2_r = gal_E1fb_tcm2_states(:,1:3);
gal_E1fb_tcm2_v = gal_E1fb_tcm2_states(:,4:6);

[E1ca_tcm2_inc,E1ca_tcm2_alt, E1ca_tcm2_et] = closest_approach(gal_E1fb_tcm2_etvec, gal_E1fb_tcm2_r, gal_E1fb_tcm2_v ,rad_earth);
E1ca_tcm2_utc = cspice_et2utc(E1ca_tcm2_et, 'ISOC', 3);
E1ca_tcm2_inc_deg = rad2deg(E1ca_tcm2_inc);

%Add data to flyby table
row_fb = row_fb+1;
flyby_table(row_fb,:) = {2, "E1fb", ...
                       tcm2_E1fb(1), tcm2_E1fb(2), tcm2_E1fb(3), norm(tcm2_E1fb), ...
                       E1ca_tcm2_utc, E1ca_tcm2_alt, E1ca_tcm2_inc_deg, 0};

%--------------------------------------------------------------------------
% 2.6 Plotting Graphs for Phase 2: Venus to Earth-1
%--------------------------------------------------------------------------

%Plotting Venus to Earth-1 Cruise (SPICE vs Model)
[gal_spice_v2e1_states, gal_spice_v2e1_etvec] = spice_states(utc_venus_dep, utc_earth1_arr, 2500, sun_ID , gal_ID);
gal_spice_v2e1_r = gal_spice_v2e1_states(:,1:3);
gal_spice_v2e1_v = gal_spice_v2e1_states(:,4:6);

figure(5);
plot3(v2e1_r(:,1), v2e1_r(:,2), v2e1_r(:,3),'-b');
hold on;
plot3(gal_spice_v2e1_r(:,1), gal_spice_v2e1_r(:,2), gal_spice_v2e1_r(:,3),'--r');
plot3(0, 0, 0, '.','Color', '[1, 0.8, 0.1]','MarkerSize',marker*3);
plot3(earth1_arr_state(1), earth1_arr_state(2), earth1_arr_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(venus_dep_state(1), venus_dep_state(2), venus_dep_state(3), '.','Color', '[0.7, 0, 1]','MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model','SPICE', 'Sun', 'Earth', 'Venus');
title('Venus to Earth-1 Cruise (SPICE vs Model)');

%Plotting Partial Derivatives wrt Velocity Perturbations
delta_Vx = bplane_E1fb_partials(1,:); 
delta_Vy = bplane_E1fb_partials(2,:); 
delta_Vz = bplane_E1fb_partials(3,:);
del_BT_Vx_arr = bplane_E1fb_partials(4,:);
del_BT_Vy_arr = bplane_E1fb_partials(5,:); 
del_BT_Vz_arr = bplane_E1fb_partials(6,:);
del_BR_Vx_arr = bplane_E1fb_partials(7,:); 
del_BR_Vy_arr = bplane_E1fb_partials(8,:);
del_BR_Vz_arr = bplane_E1fb_partials(9,:);

figure(6);
t = tiledlayout(3,2);
title(t,'Earth-1 B-Plane Partial Derivatives');
nexttile;
plot(log10(delta_Vx), del_BT_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vx), del_BR_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BT_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BR_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BT_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_z}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BR_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_z}$','Interpreter','latex');

%Plotting b-plane with and without TCM
W = null(S);
[P,Q] = meshgrid(-2e5:1e5:2e5); 
X = 0+W(1,1)*P+W(1,2)*Q; 
Y = 0+W(2,1)*P+W(2,2)*Q; 
Z = 0+W(3,1)*P+W(3,2)*Q;

R_plot = [0 0 0; R.*1e5];
T_plot = [0 0 0; T.*1e5];
S_plot = [0 0 0; S.*1e5];

midpoint_S = S_plot(round(length(S_plot) / 2), :);
midpoint_T = T_plot(round(length(T_plot) / 2), :);
midpoint_R = R_plot(round(length(R_plot) / 2), :);

figure(7);
plot3(gal_E1fb_r(:,1), gal_E1fb_r(:,2), gal_E1fb_r(:,3),'-b');
hold on;
plot3(gal_E1fb_tcm2_r(:,1), gal_E1fb_tcm2_r(:,2), gal_E1fb_tcm2_r(:,3),'-r');
plot3(0,0,0,'.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
surf(X,Y,Z, 'EdgeColor', 'none', 'FaceColor', '[0.3, 1, 0.8]','FaceAlpha', 0.5);
plot3(S_plot(:,1),S_plot(:,2),S_plot(:,3),'-k');
plot3(T_plot(:,1),T_plot(:,2),T_plot(:,3),'-k');
plot3(R_plot(:,1),R_plot(:,2),R_plot(:,3),'-k');
text(midpoint_S(1), midpoint_S(2), midpoint_S(3), ' S', 'FontSize', 8, 'Color', 'black');
text(midpoint_T(1), midpoint_T(2), midpoint_T(3), ' T', 'FontSize', 8, 'Color', 'black');
text(midpoint_R(1), midpoint_R(2), midpoint_R(3), ' R', 'FontSize', 8, 'Color', 'black');
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Without TCM','With TCM','Earth','B-plane');
title('B-plane of Earth-1 With and Without TCM');

%Plotting Galileo closest approach to Earth-1
utc_E1ca_start = '1990-12-08T16:35:00';
et_E1ca_start = cspice_str2et(utc_E1ca_start);
utc_E1ca_stop = cspice_et2utc(et_E1fb_stop, 'ISOC',1);

gal_E1ca_tcm2_r = slice_flyby(gal_E1fb_tcm2_r, gal_E1fb_tcm2_etvec, et_E1ca_start);

[gal_spice_E1ca_states, gal_spice_E1ca_etvec] = spice_states(utc_E1ca_start, utc_E1ca_stop, 2500, earth_ID , gal_ID);
gal_spice_E1ca_r = gal_spice_E1ca_states(:,1:3);
gal_spice_E1ca_v = gal_spice_E1ca_states(:,4:6);

[gal_spice_E1ca_inc, gal_spice_E1ca_alt, gal_spice_E1ca_et] = closest_approach(gal_spice_E1ca_etvec, gal_spice_E1ca_r, gal_spice_E1ca_v, rad_earth);
gal_spice_E1ca_utc = cspice_et2utc(gal_spice_E1ca_et, 'ISOC', 3);
gal_spice_E1ca_inc_deg = rad2deg(gal_spice_E1ca_inc);

figure(8)
plot3(gal_E1ca_tcm2_r(:,1), gal_E1ca_tcm2_r(:,2),gal_E1ca_tcm2_r(:,3) ,'-b');
hold on;
plot3(gal_spice_E1ca_r(:,1), gal_spice_E1ca_r(:,2), gal_spice_E1ca_r(:,3), '--r');
plot3(0,0,0,'.', 'Color', '[0.1, 0.8, 1]', 'MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model', 'SPICE', 'Earth');
title('Closest Approach to Earth-1 - SPICE vs Model');

%**************************************************************************
%                                  PHASE 3: 
%                               EARTH-1 TO EARTH-2
%**************************************************************************
disp('-> Phase 3: Earth-1 to Earth-2 ✓');

%--------------------------------------------------------------------------
% 3.1 Earth-1 to Earth-2 Cruise
%--------------------------------------------------------------------------

%Getting states of Earth-1 at departure and Earth-2 at arrival
utc_earth1_dep = '1990-12-08T20:35:00';
utc_gaspra_arr = '1991-10-29T22:37:44';
utc_earth2_arr = '1992-12-08T15:35:00';
[earth1_dep_state, et_earth1_dep] = spice_states(utc_earth1_dep, utc_earth1_dep, 1, sun_ID, earth_ID);
[gaspra_arr_state, et_gaspra_arr] = spice_states(utc_gaspra_arr, utc_gaspra_arr, 1, sun_ID, gaspra_ID);
[gaspra_dep_state, et_gaspra_dep] = spice_states(utc_gaspra_arr, utc_gaspra_arr, 1, sun_ID, gaspra_ID);
[earth2_arr_state, et_earth2_arr] = spice_states(utc_earth2_arr, utc_earth2_arr, 1, sun_ID, earth_ID);

%Calculating time of flight between Earth-1 departure and Gaspra arrival
JD_earth1_dep = julian_date(utc_earth1_dep);
JD_gaspra_arr = julian_date(utc_gaspra_arr);
tof_e12g = (JD_gaspra_arr - JD_earth1_dep)*(24*60*60);

%Calculating time of flight between Gaspra departure and Earth-2 arrival
JD_gaspra_dep = julian_date(utc_gaspra_arr);
JD_earth2_arr = julian_date(utc_earth2_arr);
tof_g2e2 = (JD_earth2_arr - JD_gaspra_dep)*(24*60*60);

%Solving Lambert's Problem to find velocity vectors at Earth-1 departure and Gaspra arrival
num_rev_e12g = 0;
[v_hel_E1dep, v_hel_Garr] = glambert(mu_sun, earth1_dep_state, gaspra_arr_state, tof_e12g, num_rev_e12g);

%Solving Lambert's Problem to find velocity vectors at Gaspra departure and Earth-2 arrival
num_rev_g2e2 = 0;
[v_hel_Gdep, v_hel_E2arr] = glambert(mu_sun, gaspra_arr_state, earth2_arr_state, tof_g2e2, num_rev_g2e2);

%Orbit Propagation using Two Body Problem ignoring pertubations
[e12g_states, e12g_etvec] = twobody(mu_sun, et_earth1_dep, et_gaspra_arr, [earth1_dep_state(1:3) v_hel_E1dep]);
e12g_r = e12g_states(:,1:3);
e12g_v = e12g_states(:,4:6);

%Orbit Propagation using Two Body Problem ignoring pertubations
[g2e2_states, g2e2_etvec] = twobody(mu_sun, et_gaspra_dep, et_earth2_arr, [gaspra_dep_state(1:3) v_hel_Gdep]);
g2e2_r = g2e2_states(:,1:3);
g2e2_v = g2e2_states(:,4:6);

%--------------------------------------------------------------------------
% 3.2 Earth-1 Departure
%--------------------------------------------------------------------------

%Calulating heliocentric velocity change due to flyby
delta_v_E1dep = v_hel_E1dep -v_hel_E1arr;
delta_v_E1dep_mag = norm(v_hel_E1dep) - norm(v_hel_E1arr);

%Calculating outgoing v-infinity vector
v_inf_E1dep = v_hel_E1dep - earth1_dep_state(4:6);
V_inf_E1dep_mag = norm(v_inf_E1dep);

%Add data to cruise table
row_c = row_c+1;
cruise_table(row_c,:) = {3, "E1dep", ...
                       earth1_dep_state(4), earth1_dep_state(5), earth1_dep_state(6), norm(earth1_dep_state(4:6)), ...
                       v_hel_E1dep(1), v_hel_E1dep(2), v_hel_E1dep(3), norm(v_hel_E1dep), ...
                       v_inf_E1dep(1), v_inf_E1dep(2), v_inf_E1dep(3), V_inf_E1dep_mag};

%Add data to flyby table
flyby_table.flyby_deltaV(row_fb) = delta_v_E1dep_mag;

%--------------------------------------------------------------------------
% 3.3 Gaspra Arrival
%--------------------------------------------------------------------------

%Calculating incoming v-infinity vector
v_inf_Garr = v_hel_Garr - gaspra_arr_state(4:6);
V_inf_Garr_mag = norm(v_inf_Garr);

%Add data to cruise table
row_c = row_c + 1;
cruise_table(row_c,:) = {3, "Garr", ...
                       gaspra_arr_state(4), gaspra_arr_state(5), gaspra_arr_state(6), norm(gaspra_arr_state(4:6)), ...
                       v_hel_Garr(1), v_hel_Garr(2), v_hel_Garr(3), norm(v_hel_Garr), ...
                       v_inf_Garr(1), v_inf_Garr(2), v_inf_Garr(3), V_inf_Garr_mag};

%--------------------------------------------------------------------------
% 3.4 Gaspra Departure
%--------------------------------------------------------------------------

%Calulating heliocentric velocity change due to flyby
delta_v_Gdep = v_hel_Gdep -v_hel_Garr;
delta_v_Gdep_mag = norm(v_hel_Gdep) - norm(v_hel_Garr);

%Calculating outgoing v-infinity vector
v_inf_Gdep = v_hel_Gdep - gaspra_dep_state(4:6);
V_inf_Gdep_mag = norm(v_inf_Gdep);

%Add data to cruise table
row_c = row_c+1;
cruise_table(row_c,:) = {3, "Gdep", ...
                       gaspra_dep_state(4), gaspra_dep_state(5), gaspra_dep_state(6), norm(gaspra_dep_state(4:6)), ...
                       v_hel_Gdep(1), v_hel_Gdep(2), v_hel_Gdep(3), norm(v_hel_Gdep), ...
                       v_inf_Gdep(1), v_inf_Gdep(2), v_inf_Gdep(3), V_inf_Gdep_mag};

%--------------------------------------------------------------------------
% 3.5 Earth-2 Arrival
%--------------------------------------------------------------------------

%Calculating incoming v-infinity vector
v_inf_E2arr = v_hel_E2arr - earth2_arr_state(4:6);
V_inf_E2arr_mag = norm(v_inf_E2arr);

%Finding the state vector of spacecraft at Earth-2 SOI
for i = 1:length(g2e2_etvec)
    utc_earth2 = cspice_et2utc(g2e2_etvec(i),'ISOC', 3);
    [state_earth2, et_earth2] = spice_states(utc_earth2, utc_earth2, 1, sun_ID, earth_ID);
    earth2_SOI = norm(state_earth2(1:3)) * (mu_earth/mu_sun)^(2/5);

    position_gal = g2e2_r(i,:) - state_earth2(1:3);
    velocity_gal = g2e2_v(i,:) - state_earth2(4:6);

    dist_earth2_gal = norm(position_gal);

    if dist_earth2_gal <= earth2_SOI
        break;
    end
end

if i == length(g2e2_etvec)
    error('Galileo does not intercept Earth-2 SOI');
else
    gal_E2SOI = [position_gal velocity_gal];
    et_E2SOI_gal = g2e2_etvec(i);
end

%Add data to cruise table
row_c = row_c + 1;
cruise_table(row_c,:) = {3, "E2arr", ...
                       earth2_arr_state(4), earth2_arr_state(5), earth2_arr_state(6), norm(earth2_arr_state(4:6)), ...
                       v_hel_E2arr(1), v_hel_E2arr(2), v_hel_E2arr(3), norm(v_hel_E2arr), ...
                       v_inf_E2arr(1), v_inf_E2arr(2), v_inf_E2arr(3), V_inf_E2arr_mag};

%--------------------------------------------------------------------------
% 3.6 B-Plane Targeting of Earth-2 Fly-by 
%--------------------------------------------------------------------------

%Getting current b-plane parameters
[bplane_E2fb,S, R, T] = b_plane(mu_earth,gal_E2SOI(1:3), gal_E2SOI(4:6));

%Getting target b-plane
E2fb_inc = deg2rad(138.6); %ecliptic frame
E2fb_alt = 300;
E2fb_R_per = E2fb_alt + rad_earth;
[BT_star,BR_star, ~, ~, ~] = b_plane_target(mu_earth, gal_E2SOI(1:3), gal_E2SOI(4:6), S, E2fb_inc, E2fb_R_per);

%Calculating b-plane Jacobian and its pseudoinverse
[bplane_E2fb_partials, J] = b_plane_jacobian(mu_earth, gal_E2SOI(1:3), gal_E2SOI(4:6));
J_inv = transpose(J) * (J * transpose(J))^(-1);

%Calculating delta B
BT = bplane_E2fb(8);
BR = bplane_E2fb(9);
delta_B_earth2 = [BT_star - BT;
                 BR_star - BR];

%Calculating delta V for TCM-3
tcm3_E2fb = transpose(J_inv * delta_B_earth2);

%Getting updated state of Galileo at Gaspra SOI after TCM-3
gal_E2SOI_tcm3 = [gal_E2SOI(1:3) gal_E2SOI(4:6)+tcm3_E2fb];

%--------------------------------------------------------------------------
% 3.7 Closest Approach Parameters of Earth-2 Fly-By 
%--------------------------------------------------------------------------

%Specifying start and stop propagation times in Earth-2 SOI
utc_E2fb_stop = '1992-12-08T19:35:00'; 
et_E2fb_stop = cspice_str2et(utc_E2fb_stop);
et_E2fb_start = et_E2SOI_gal;

%Propagating Galileo trajectory in Gaspra SOI without TCM-3
[gal_E2fb_states, gal_E2fb_etvec] = twobody(mu_earth, et_E2fb_start, et_E2fb_stop, gal_E2SOI);
gal_E2fb_r = gal_E2fb_states(:,1:3);
gal_E2fb_v = gal_E2fb_states(:,4:6);

[E2ca_inc, E2ca_alt, E2ca_et] = closest_approach(gal_E2fb_etvec, gal_E2fb_r, gal_E2fb_v, rad_earth);
E2ca_utc = cspice_et2utc(E2ca_et, 'ISOC', 3);
E2ca_inc_deg = rad2deg(E2ca_inc);

%Propagating Galileo trajectory in Gaspra SOI with TCM-3
[gal_E2fb_tcm3_states, gal_E2fb_tcm3_etvec] = twobody(mu_earth, et_E2fb_start, et_E2fb_stop, gal_E2SOI_tcm3);
gal_E2fb_tcm3_r = gal_E2fb_tcm3_states(:,1:3);
gal_E2fb_tcm3_v = gal_E2fb_tcm3_states(:,4:6);

[E2ca_tcm3_inc,E2ca_tcm3_alt, E2ca_tcm3_et] = closest_approach(gal_E2fb_tcm3_etvec, gal_E2fb_tcm3_r, gal_E2fb_tcm3_v ,rad_earth);
E2ca_tcm3_utc = cspice_et2utc(E2ca_tcm3_et, 'ISOC', 3);
E2ca_tcm3_inc_deg = rad2deg(E2ca_tcm3_inc);

%Add data to flyby table
row_fb = row_fb+1;
flyby_table(row_fb,:) = {3, "E2fb", ...
                       tcm3_E2fb(1), tcm3_E2fb(2), tcm3_E2fb(3), norm(tcm3_E2fb), ...
                       E2ca_tcm3_utc, E2ca_tcm3_alt, E2ca_tcm3_inc_deg, 0};

%--------------------------------------------------------------------------
% 3.8 Plotting Graphs for Phase 3: Earth-1 to Earth-2
%--------------------------------------------------------------------------

%Plotting Earth-1 to Earth-2 Cruise (SPICE vs Model)
[gal_spice_e12e2_states, gal_spice_e12e2_etvec] = spice_states(utc_earth1_dep, utc_earth2_arr, 2500, sun_ID , gal_ID);
gal_spice_e12e2_r = gal_spice_e12e2_states(:,1:3);
gal_spice_e12e2_v = gal_spice_e12e2_states(:,4:6);

figure(9);
plot3(e12g_r(:,1), e12g_r(:,2), e12g_r(:,3),'-b');
hold on;
plot3(g2e2_r(:,1), g2e2_r(:,2), g2e2_r(:,3),'-b');
plot3(gal_spice_e12e2_r(:,1), gal_spice_e12e2_r(:,2), gal_spice_e12e2_r(:,3),'--r');
plot3(0, 0, 0, '.','Color', '[1, 0.8, 0.1]','MarkerSize',marker*3);
plot3(earth1_dep_state(1), earth1_dep_state(2), earth1_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker*2);
plot3(earth2_arr_state(1), earth2_arr_state(2), earth2_arr_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker*2);
plot3(gaspra_dep_state(1), gaspra_dep_state(2), gaspra_dep_state(3), '.','Color', '[1, 0.4, 0.8]','MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model','','SPICE', 'Sun', 'Earth','', 'Gaspra');
title('Earth-1 to Earth-2 Cruise (SPICE vs Model)');

%Plotting Partial Derivatives wrt Velocity Perturbations
delta_Vx = bplane_E2fb_partials(1,:); 
delta_Vy = bplane_E2fb_partials(2,:); 
delta_Vz = bplane_E2fb_partials(3,:);
del_BT_Vx_arr = bplane_E2fb_partials(4,:);
del_BT_Vy_arr = bplane_E2fb_partials(5,:); 
del_BT_Vz_arr = bplane_E2fb_partials(6,:);
del_BR_Vx_arr = bplane_E2fb_partials(7,:); 
del_BR_Vy_arr = bplane_E2fb_partials(8,:);
del_BR_Vz_arr = bplane_E2fb_partials(9,:);

figure(10);
t = tiledlayout(3,2);
title(t,'Earth-2 B-Plane Partial Derivatives');
nexttile;
plot(log10(delta_Vx), del_BT_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vx), del_BR_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BT_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BR_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BT_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_z}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BR_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_z}$','Interpreter','latex');

%Plotting b-plane with and without TCM
W = null(S);
[P,Q] = meshgrid(-2e5:1e5:2e5); 
X = 0+W(1,1)*P+W(1,2)*Q; 
Y = 0+W(2,1)*P+W(2,2)*Q; 
Z = 0+W(3,1)*P+W(3,2)*Q;

R_plot = [0 0 0; R.*1e5];
T_plot = [0 0 0; T.*1e5];
S_plot = [0 0 0; S.*1e5];

midpoint_S = S_plot(round(length(S_plot) / 2), :);
midpoint_T = T_plot(round(length(T_plot) / 2), :);
midpoint_R = R_plot(round(length(R_plot) / 2), :);

figure(11);
plot3(gal_E2fb_r(:,1), gal_E2fb_r(:,2), gal_E2fb_r(:,3),'-b');
hold on;
plot3(gal_E2fb_tcm3_r(:,1), gal_E2fb_tcm3_r(:,2), gal_E2fb_tcm3_r(:,3),'-r');
plot3(0,0,0,'.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
surf(X,Y,Z, 'EdgeColor', 'none', 'FaceColor', '[0.3, 1, 0.8]','FaceAlpha', 0.5);
plot3(S_plot(:,1),S_plot(:,2),S_plot(:,3),'-k');
plot3(T_plot(:,1),T_plot(:,2),T_plot(:,3),'-k');
plot3(R_plot(:,1),R_plot(:,2),R_plot(:,3),'-k');
text(midpoint_S(1), midpoint_S(2), midpoint_S(3), ' S', 'FontSize', 8, 'Color', 'black');
text(midpoint_T(1), midpoint_T(2), midpoint_T(3), ' T', 'FontSize', 8, 'Color', 'black');
text(midpoint_R(1), midpoint_R(2), midpoint_R(3), ' R', 'FontSize', 8, 'Color', 'black');
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Without TCM','With TCM','Earth','B-plane');
title('B-plane of Earth-2 With and Without TCM');

%Plotting Galileo closest approach to Gaspra
utc_E2ca_start = '1992-12-08T11:35:00'; 
et_E2ca_start = cspice_str2et(utc_E2ca_start);
utc_E2ca_stop = cspice_et2utc(et_E2fb_stop, 'ISOC',1);

gal_E2ca_tcm3_r = slice_flyby(gal_E2fb_tcm3_r, gal_E2fb_tcm3_etvec, et_E2ca_start);

[gal_spice_E2ca_states, gal_spice_E2ca_etvec] = spice_states(utc_E2ca_start, utc_E2ca_stop, 2500, earth_ID , gal_ID);
gal_spice_E2ca_r = gal_spice_E2ca_states(:,1:3);
gal_spice_E2ca_v = gal_spice_E2ca_states(:,4:6);

[gal_spice_E2ca_inc, gal_spice_E2ca_alt, gal_spice_E2ca_et] = closest_approach(gal_spice_E2ca_etvec, gal_spice_E2ca_r, gal_spice_E2ca_v, rad_earth);
gal_spice_E2ca_utc = cspice_et2utc(gal_spice_E2ca_et, 'ISOC', 3);
gal_spice_E2ca_inc_deg = rad2deg(gal_spice_E2ca_inc);

figure(12)
plot3(gal_E2ca_tcm3_r(:,1), gal_E2ca_tcm3_r(:,2),gal_E2ca_tcm3_r(:,3) ,'-b');
hold on;
plot3(gal_spice_E2ca_r(:,1), gal_spice_E2ca_r(:,2), gal_spice_E2ca_r(:,3), '--r');
plot3(0,0,0,'.', 'Color', '[0.1, 0.8, 1]', 'MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model', 'SPICE', 'Earth');
title('Closest Approach to Earth-2 - SPICE vs Model');


%**************************************************************************
%                                  PHASE 4: 
%                               EARTH-2 TO JUPITER
%**************************************************************************
disp('-> Phase 4: Earth-2 to Jupiter ✓');

%--------------------------------------------------------------------------
% 4.1 Earth-2 to Jupiter Cruise
%--------------------------------------------------------------------------

%Getting states of Earth-2 at departure, Ida and Jupiter at arrival
utc_earth2_dep = '1992-12-08T15:35:00';
utc_ida_arr = '1993-08-28T16:53:04';
utc_jupiter_arr = '1995-12-07T21:54:00';
[earth2_dep_state, et_earth2_dep] = spice_states(utc_earth2_dep, utc_earth2_dep, 1, sun_ID, earth_ID);
[ida_arr_state, et_ida_arr] = spice_states(utc_ida_arr, utc_ida_arr, 1, sun_ID, ida_ID);
[ida_dep_state, et_ida_dep] = spice_states(utc_ida_arr, utc_ida_arr, 1, sun_ID, ida_ID);
[jupiter_arr_state, et_jupiter_arr] = spice_states(utc_jupiter_arr, utc_jupiter_arr, 1, sun_ID, jupiter_ID);

%Calculating time of flight between Earth-2 departure and Ida arrival
JD_earth2_dep = julian_date(utc_earth2_dep);
JD_ida_arr = julian_date(utc_ida_arr);
JD_ida_dep = julian_date(utc_ida_arr);
JD_jupiter_arr = julian_date(utc_jupiter_arr);
tof_e22i = (JD_ida_arr - JD_earth2_dep)*(24*60*60);
tof_i2j = (JD_jupiter_arr - JD_ida_dep)*(24*60*60);

%Solving Lambert's Problem to find velocity vectors at Earth-2 departure and Ida arrival
num_rev_e22i = 0;
[v_hel_E2dep, v_hel_Iarr] = glambert(mu_sun, earth2_dep_state, ida_arr_state, tof_e22i, num_rev_e22i);

%Solving Lambert's Problem to find velocity vectors at Ida departure and Jupiter arrival
num_rev_i2j = 0;
[v_hel_Idep, v_hel_Jarr] = glambert(mu_sun, ida_dep_state, jupiter_arr_state, tof_i2j, num_rev_i2j);

%Orbit Propagation using Two Body Problem ignoring pertubations
[e22i_states, e22i_etvec] = twobody(mu_sun, et_earth2_dep, et_ida_arr, [earth2_dep_state(1:3) v_hel_E2dep]);
e22i_r = e22i_states(:,1:3);
e22i_v = e22i_states(:,4:6);

%Orbit Propagation using Two Body Problem ignoring pertubations
[i2j_states, i2j_etvec] = twobody(mu_sun, et_ida_dep, et_jupiter_arr, [ida_dep_state(1:3) v_hel_Idep]);
i2j_r = i2j_states(:,1:3);
i2j_v = i2j_states(:,4:6);

%--------------------------------------------------------------------------
% 4.2 Earth-2 Departure
%--------------------------------------------------------------------------

%Calulating heliocentric velocity change due to flyby
delta_v_E2dep = v_hel_E2dep -v_hel_E2arr;
delta_v_E2dep_mag = norm(v_hel_E2dep) - norm(v_hel_E2arr);

%Calculating outgoing v-infinity vector
v_inf_E2dep = v_hel_E2dep - earth2_dep_state(4:6);
V_inf_E2dep_mag = norm(v_inf_E2dep);

%Add data to cruise table
row_c = row_c+1;
cruise_table(row_c,:) = {4, "E2dep", ...
                       earth2_dep_state(4), earth2_dep_state(5), earth2_dep_state(6), norm(earth2_dep_state(4:6)), ...
                       v_hel_E2dep(1), v_hel_E2dep(2), v_hel_E2dep(3), norm(v_hel_E2dep), ...
                       v_inf_E2dep(1), v_inf_E2dep(2), v_inf_E2dep(3), V_inf_E2dep_mag};

%Add data to flyby table
flyby_table.flyby_deltaV(row_fb) = delta_v_E2dep_mag;

%--------------------------------------------------------------------------
% 4.3 Ida Arrival
%--------------------------------------------------------------------------

%Calculating incoming v-infinity vector
v_inf_Iarr = v_hel_Iarr - ida_arr_state(4:6);
V_inf_Iarr_mag = norm(v_inf_Iarr);

%Add data to cruise table
row_c = row_c + 1;
cruise_table(row_c,:) = {4, "Iarr", ...
                       ida_arr_state(4), ida_arr_state(5), ida_arr_state(6), norm(ida_arr_state(4:6)), ...
                       v_hel_Iarr(1), v_hel_Iarr(2), v_hel_Iarr(3), norm(v_hel_Iarr), ...
                       v_inf_Iarr(1), v_inf_Iarr(2), v_inf_Iarr(3), V_inf_Iarr_mag};


%--------------------------------------------------------------------------
% 4.4 Ida Departure
%--------------------------------------------------------------------------

%Calulating heliocentric velocity change due to flyby
delta_v_Idep = v_hel_Idep -v_hel_Iarr;
delta_v_Idep_mag = norm(v_hel_Idep) - norm(v_hel_Iarr);

%Calculating outgoing v-infinity vector
v_inf_Idep = v_hel_Idep - ida_dep_state(4:6);
V_inf_Idep_mag = norm(v_inf_Idep);

%Add data to cruise table
row_c = row_c+1;
cruise_table(row_c,:) = {4, "Idep", ...
                       ida_dep_state(4), ida_dep_state(5), ida_dep_state(6), norm(ida_dep_state(4:6)), ...
                       v_hel_Idep(1), v_hel_Idep(2), v_hel_Idep(3), norm(v_hel_Idep), ...
                       v_inf_Idep(1), v_inf_Idep(2), v_inf_Idep(3), V_inf_Idep_mag};

%--------------------------------------------------------------------------
% 4.5 Jupiter Arrival
%--------------------------------------------------------------------------

%Calculating incoming v-infinity vector
v_inf_Jarr = v_hel_Jarr - jupiter_arr_state(4:6);
V_inf_Jarr_mag = norm(v_inf_Jarr);

%Finding the state vector of spacecraft at Jupiter SOI
for i = 2210:length(i2j_etvec)
    utc_jupiter = cspice_et2utc(i2j_etvec(i),'ISOC', 3);
    [state_jupiter, et_jupiter] = spice_states(utc_jupiter, utc_jupiter, 1, sun_ID, jupiter_ID);
    jupiter_SOI = norm(state_jupiter(1:3)) * (mu_jupiter/mu_sun)^(2/5);

    position_gal = i2j_r(i,:) - state_jupiter(1:3);
    velocity_gal = i2j_v(i,:) - state_jupiter(4:6);

    dist_jupiter_gal = norm(position_gal);

    if dist_jupiter_gal <= jupiter_SOI
        break;
    end
end

if i == length(i2j_etvec)
    error('Galileo does not intercept Jupiter SOI');
else
    i2jsoi_r = i2j_r(1:i+100,:);
    i2jsoi_v = i2j_v(1:i+100,:);
    i2jsoi_etvec = i2j_etvec(1:i+100);

    gal_JSOI = [position_gal velocity_gal];
    et_JSOI_gal = i2j_etvec(i);
end

%Add data to cruise table
row_c = row_c + 1;
cruise_table(row_c,:) = {4, "Jarr", ...
                       jupiter_arr_state(4), jupiter_arr_state(5), jupiter_arr_state(6), norm(jupiter_arr_state(4:6)), ...
                       v_hel_Jarr(1), v_hel_Jarr(2), v_hel_Jarr(3), norm(v_hel_Jarr), ...
                       v_inf_Jarr(1), v_inf_Jarr(2), v_inf_Jarr(3), V_inf_Jarr_mag};


%--------------------------------------------------------------------------
% 4.6 B-Plane Targeting of Jupiter Fly-by 
%--------------------------------------------------------------------------

%Getting current b-plane parameters
[bplane_Jfb,S, R, T] = b_plane(mu_jupiter ,gal_JSOI(1:3), gal_JSOI(4:6));

%Getting target b-plane
Jfb_inc = deg2rad(26.5); %ecliptic frame
Jfb_alt = 2.1610e5;
Jfb_R_per = Jfb_alt + rad_jupiter;
[BT_star,BR_star, ~, ~, ~] = b_plane_target(mu_jupiter, gal_JSOI(1:3), gal_JSOI(4:6), S, Jfb_inc, Jfb_R_per);

%Calculating b-plane Jacobian and its pseudoinverse
[bplane_Jfb_partials, J] = b_plane_jacobian(mu_jupiter, gal_JSOI(1:3), gal_JSOI(4:6));
J_inv = transpose(J) * (J * transpose(J))^(-1);

%Calculating delta B
BT = bplane_Jfb(8);
BR = bplane_Jfb(9);
delta_B_jupiter = [BT_star - BT;
                 BR_star - BR];

%Calculating delta V for TCM-4
tcm4_Jfb = transpose(J_inv * delta_B_jupiter);

%Getting updated state of Galileo at Jupiter after TCM-4
gal_JSOI_tcm4 = [gal_JSOI(1:3) gal_JSOI(4:6)+tcm4_Jfb];

%--------------------------------------------------------------------------
% 4.7 Closest Approach Parameters of Jupiter Fly-By 
%--------------------------------------------------------------------------

%Specifying start and stop propagation times in Jupiter SOI
utc_Jfb_stop = '1995-11-25T21:40:00'; 
et_Jfb_stop = cspice_str2et(utc_Jfb_stop);
et_Jfb_start = et_JSOI_gal;

%Propagating Galileo trajectory in Jupiter SOI without TCM-4
[gal_Jfb_states, gal_Jfb_etvec] = twobody(mu_jupiter, et_Jfb_start, et_Jfb_stop, gal_JSOI);
gal_Jfb_r = gal_Jfb_states(:,1:3);
gal_Jfb_v = gal_Jfb_states(:,4:6);

[Jca_inc, Jca_alt, Jca_et] = closest_approach(gal_Jfb_etvec, gal_Jfb_r, gal_Jfb_v, rad_jupiter);
Jca_utc = cspice_et2utc(Jca_et, 'ISOC', 3);
Jca_inc_deg = rad2deg(Jca_inc);

%Propagating Galileo trajectory in Jupiter SOI with TCM-4
[gal_Jfb_tcm4_states, gal_Jfb_tcm4_etvec] = twobody(mu_jupiter, et_Jfb_start, et_Jfb_stop, gal_JSOI_tcm4);
gal_Jfb_tcm4_r = gal_Jfb_tcm4_states(:,1:3);
gal_Jfb_tcm4_v = gal_Jfb_tcm4_states(:,4:6);

[Jca_tcm4_inc,Jca_tcm4_alt, Jca_tcm4_et] = closest_approach(gal_Jfb_tcm4_etvec, gal_Jfb_tcm4_r, gal_Jfb_tcm4_v ,rad_jupiter);
Jca_tcm4_utc = cspice_et2utc(Jca_tcm4_et, 'ISOC', 3);
Jca_tcm4_inc_deg = rad2deg(Jca_tcm4_inc);

%Add data to flyby table
row_fb = row_fb+1;
flyby_table(row_fb,:) = {4, "Jca", ...
                       tcm4_Jfb(1), tcm4_Jfb(2), tcm4_Jfb(3), norm(tcm4_Jfb), ...
                       Jca_tcm4_utc, Jca_tcm4_alt, Jca_tcm4_inc_deg, 0};

%--------------------------------------------------------------------------
% 4.8 Plotting Graphs for Phase 4: Earth-2 to Jupiter
%--------------------------------------------------------------------------

%Plotting Earth-2 to Jupiter Cruise (SPICE vs Model)
[gal_spice_e22j_states, gal_spice_e22j_etvec] = spice_states(utc_earth2_dep, utc_jupiter_arr, 2500, sun_ID , gal_ID);
gal_spice_e22j_r = gal_spice_e22j_states(:,1:3);
gal_spice_e22j_v = gal_spice_e22j_states(:,4:6);

figure(13);
plot3(e22i_r(:,1), e22i_r(:,2), e22i_r(:,3),'-b');
hold on;
plot3(i2j_r(:,1), i2j_r(:,2), i2j_r(:,3),'-b');
plot3(gal_spice_e22j_r(:,1), gal_spice_e22j_r(:,2), gal_spice_e22j_r(:,3),'--r');
plot3(0, 0, 0, '.','Color', '[1, 0.8, 0.1]','MarkerSize',marker*3);
plot3(earth2_dep_state(1), earth2_dep_state(2), earth2_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker*2);
plot3(jupiter_arr_state(1), jupiter_arr_state(2), jupiter_arr_state(3), '.','Color', '[0.9, 0.6, 0]','MarkerSize', marker*2);
plot3(ida_dep_state(1), ida_dep_state(2), ida_dep_state(3), '.','Color', '[0.4, 0.9, 0.2]','MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model','','SPICE', 'Sun', 'Earth','Jupiter', 'Ida');
title('Earth-2 to Jupiter Cruise (SPICE vs Model)');

%Plotting Partial Derivatives wrt Velocity Perturbations
delta_Vx = bplane_Jfb_partials(1,:); 
delta_Vy = bplane_Jfb_partials(2,:); 
delta_Vz = bplane_Jfb_partials(3,:);
del_BT_Vx_arr = bplane_Jfb_partials(4,:);
del_BT_Vy_arr = bplane_Jfb_partials(5,:); 
del_BT_Vz_arr = bplane_Jfb_partials(6,:);
del_BR_Vx_arr = bplane_Jfb_partials(7,:); 
del_BR_Vy_arr = bplane_Jfb_partials(8,:);
del_BR_Vz_arr = bplane_Jfb_partials(9,:);

figure(14);
t = tiledlayout(3,2);
title(t,'Jupiter B-Plane Partial Derivatives');
nexttile;
plot(log10(delta_Vx), del_BT_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vx), del_BR_Vx_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_x}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BT_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vy), del_BR_Vy_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_y}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BT_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_T / \delta \Delta V_z}$','Interpreter','latex');

nexttile;
plot(log10(delta_Vz), del_BR_Vz_arr);
grid minor;
xlabel('Perturbation');
ylabel('$\bf{\delta B_R / \delta \Delta V_z}$','Interpreter','latex');

%Plotting b-plane with and without TCM
utc_europa_arr = '1995-12-7T14:30:00';
utc_io_arr = '1995-12-7T17:54:00';
[europa_arr_state, et_europa_arr] = spice_states(utc_europa_arr, utc_europa_arr, 1, jupiter_ID, europa_ID);
[io_arr_state, et_io_arr] = spice_states(utc_io_arr, utc_io_arr, 1, jupiter_ID, io_ID);

W = null(S);
[P,Q] = meshgrid(-2e6:1e6:2e6); 
X = 0+W(1,1)*P+W(1,2)*Q; 
Y = 0+W(2,1)*P+W(2,2)*Q; 
Z = 0+W(3,1)*P+W(3,2)*Q;

R_plot = [0 0 0; R.*1e6];
T_plot = [0 0 0; T.*1e6];
S_plot = [0 0 0; S.*1e6];

midpoint_S = S_plot(round(length(S_plot) / 2), :);
midpoint_T = T_plot(round(length(T_plot) / 2), :);
midpoint_R = R_plot(round(length(R_plot) / 2), :);

figure(15);
plot3(gal_Jfb_r(2300:end,1), gal_Jfb_r(2300:end,2), gal_Jfb_r(2300:end,3),'-b');
hold on;
plot3(gal_Jfb_tcm4_r(2300:end,1), gal_Jfb_tcm4_r(2300:end,2), gal_Jfb_tcm4_r(2300:end,3),'-r');
plot3(0,0,0,'.','Color', '[0.9, 0.6, 0]','MarkerSize', marker);
plot3(europa_arr_state(1), europa_arr_state(2), europa_arr_state(3), '.m','MarkerSize', marker);
plot3(io_arr_state(1), io_arr_state(2), io_arr_state(3), '.c','MarkerSize', marker);
surf(X,Y,Z, 'EdgeColor', 'none', 'FaceColor', '[0.3, 1, 0.8]','FaceAlpha', 0.5);
plot3(S_plot(:,1),S_plot(:,2),S_plot(:,3),'-k');
plot3(T_plot(:,1),T_plot(:,2),T_plot(:,3),'-k');
plot3(R_plot(:,1),R_plot(:,2),R_plot(:,3),'-k');
text(midpoint_S(1), midpoint_S(2), midpoint_S(3), ' S', 'FontSize', 8, 'Color', 'black');
text(midpoint_T(1), midpoint_T(2), midpoint_T(3), ' T', 'FontSize', 8, 'Color', 'black');
text(midpoint_R(1), midpoint_R(2), midpoint_R(3), ' R', 'FontSize', 8, 'Color', 'black');
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Without TCM','With TCM','Jupiter','Europa','Io','B-plane');
title('B-plane of Jupiter With and Without TCM');

%Plotting Galileo in Jupiter SOI
figure(16)
plot3(gal_Jfb_r(:,1), gal_Jfb_r(:,2),gal_Jfb_r(:,3) ,'-b');
hold on;
plot3(gal_Jfb_tcm4_r(:,1), gal_Jfb_tcm4_r(:,2),gal_Jfb_tcm4_r(:,3) ,'-r');
plot3(0,0,0,'.', 'Color', '[0.9, 0.6, 0]', 'MarkerSize', marker);
plot3(europa_arr_state(1), europa_arr_state(2), europa_arr_state(3), '.m','MarkerSize', marker);
plot3(io_arr_state(1), io_arr_state(2), io_arr_state(3), '.c','MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Without TCM','With TCM','Jupiter', 'Europa', 'Io');
title('TBD');

%Plotting Galileo closest approach to Jupiter
utc_Jca_gal = '1995-11-25T05:40:00'; 
et_Jca_gal = cspice_str2et(utc_Jca_gal);
utc_Jca_spice_start = '1995-12-7T10:30:00';
utc_Jca_spice_stop = '1995-12-08T03:54:00';

gal_Jca_tcm4_r = slice_flyby(gal_Jfb_tcm4_r, gal_Jfb_tcm4_etvec, et_Jca_gal);

[gal_spice_Jca_states, gal_spice_Jca_etvec] = spice_states(utc_Jca_spice_start, utc_Jca_spice_stop, 2500, jupiter_ID , gal_ID);
gal_spice_Jca_r = gal_spice_Jca_states(:,1:3);
gal_spice_Jca_v = gal_spice_Jca_states(:,4:6);

[gal_spice_Jca_inc, gal_spice_Jca_alt, gal_spice_Jca_et] = closest_approach(gal_spice_Jca_etvec, gal_spice_Jca_r, gal_spice_Jca_v, rad_jupiter);
gal_spice_Jca_utc = cspice_et2utc(gal_spice_Jca_et, 'ISOC', 3);
gal_spice_Jca_inc_deg = rad2deg(gal_spice_Jca_inc);

figure(17)
plot3(gal_Jca_tcm4_r(:,1), gal_Jca_tcm4_r(:,2),gal_Jca_tcm4_r(:,3) ,'-b');
hold on;
plot3(gal_spice_Jca_r(:,1), gal_spice_Jca_r(:,2), gal_spice_Jca_r(:,3), '--r');
plot3(0,0,0,'.', 'Color', '[0.9, 0.6, 0]', 'MarkerSize', marker);
plot3(europa_arr_state(1), europa_arr_state(2), europa_arr_state(3), '.m','MarkerSize', marker);
plot3(io_arr_state(1), io_arr_state(2), io_arr_state(3), '.c','MarkerSize', marker);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Model', 'SPICE', 'Jupiter', 'Europa', 'Io');
title('Closest Approach to Jupiter - SPICE vs Model');


%**************************************************************************
%                                  SAVING
%                         GALILEO TRAJECTORY DATA 
%**************************************************************************
disp('-> Saving Cruise and Flyby Data ✓');

%Save cruise table 
filepath = strcat(directory,'/cruise.xlsx');
writetable(cruise_table, filepath);

%Save flyby table 
filepath = strcat(directory,'/flyby.xlsx');
writetable(flyby_table, filepath);

%**************************************************************************
%                                  PLOTTING
%                         GALILEO INTERPLANETARY CRUISE 
%**************************************************************************
disp('-> Plotting Model Interplanetary Cruise ✓');

figure(18);
plot3(e2v_r(:,1), e2v_r(:,2), e2v_r(:,3), '-k');

hold on;
plot3(v2e1_r(:,1), v2e1_r(:,2), v2e1_r(:,3), '-k');
plot3(e12g_r(:,1), e12g_r(:,2), e12g_r(:,3), '-k');
plot3(g2e2_r(:,1), g2e2_r(:,2), g2e2_r(:,3), '-k');
plot3(e22i_r(:,1), e22i_r(:,2), e22i_r(:,3), '-k');
plot3(i2j_r(:,1), i2j_r(:,2), i2j_r(:,3), '-k');

plot3(0, 0, 0, '.','Color', '[1, 0.8, 0.1]','MarkerSize',marker*3);
plot3(venus_dep_state(1), venus_dep_state(2), venus_dep_state(3), '.','Color', '[0.7, 0, 1]','MarkerSize', marker);
plot3(gaspra_dep_state(1), gaspra_dep_state(2), gaspra_dep_state(3),'.','Color', '[1, 0.4, 0.8]','MarkerSize',marker*0.75);
plot3(earth_dep_state(1), earth_dep_state(2), earth_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(earth1_dep_state(1), earth1_dep_state(2), earth1_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(earth2_dep_state(1), earth2_dep_state(2), earth2_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(ida_dep_state(1), ida_dep_state(2), ida_dep_state(3),'.','Color', '[0.4, 0.9, 0.2]','MarkerSize',marker*0.75);
plot3(jupiter_arr_state(1), jupiter_arr_state(2), jupiter_arr_state(3), '.','Color', '[0.9, 0.6, 0]','MarkerSize', marker*2);
hold off;
view(3);
grid on;
axis equal;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Trajectory of Galileo Mission Interplanetary Cruise');
legend('Model','','','','','','Sun', 'Venus','Gaspra', 'Earth','','','Ida','Jupiter');

%**************************************************************************
%                     ACTUAL TRAJECTORY VS MODEL TRAJECTORY
%                         GALILEO INTERPLANETARY CRUISE 
%**************************************************************************
disp('-> Plotting Spice vs Model ✓');

%Concatenating the phases for the model trajectory
gal_model_x = [e2v_r(:,1); v2e1_r(:, 1); e12g_r(:, 1); g2e2_r(:, 1); e22i_r(:,1); i2j_r(:, 1)];
gal_model_y = [e2v_r(:,2); v2e1_r(:, 2); e12g_r(:, 2); g2e2_r(:, 2); e22i_r(:,2); i2j_r(:, 2)];
gal_model_z = [e2v_r(:,3); v2e1_r(:, 3); e12g_r(:, 3); g2e2_r(:, 3); e22i_r(:,3); i2j_r(:, 3)];
gal_model_v = vecnorm([e2v_v; v2e1_v; e12g_v; g2e2_v; e22i_v; i2jsoi_v],2,2);

%Concatenating the times for the model trajectory
gal_et = [e2v_etvec; v2e1_etvec; e12g_etvec; g2e2_etvec; e22i_etvec; i2jsoi_etvec];
gal_dt = NaT(length(gal_et),1);
for i= 1:length(gal_dt)
    gal_utc = cspice_et2utc(gal_et(i), 'ISOC', 3);
    gal_dt(i,1) = datetime(gal_utc);
end

%Generating the Spice states for Galileo from Earth departure to Jupiter arrival
max_pts = length(gal_model_x);
[gal_spice_states, et_gal_spice] = spice_states(utc_earth_dep, utc_jupiter_arr, max_pts, sun_ID, gal_ID);
gal_spice_r = gal_spice_states(:,1:3);
gal_spice_Vmag = vecnorm(gal_spice_states(:,4:6),2,2);

gal_spice_dt = NaT(length(et_gal_spice),1);
for i= 1:length(gal_spice_dt)
    gal_spice_utc = cspice_et2utc(et_gal_spice(i), 'ISOC', 3);
    gal_spice_dt(i,1) = datetime(gal_spice_utc);

    if gal_spice_dt(i,1) >= gal_dt(end,1)
        break;
    end
end

utc_JSOI_gal = cspice_et2utc(et_JSOI_gal, 'ISOC', 3);

%Generating Spice States for Celestial bodies
%Earth
utc_earth_year =  convertStringsToChars(string(datetime(utc_earth_dep) + days(365.25), "yyyy-MM-dd'T'HH:mm:ss"));
[earth_spice_states, et_earth_spice] = spice_states(utc_earth_dep, utc_earth_year, max_pts, sun_ID, earth_ID);
earth_spice_r = earth_spice_states(:,1:3);

%Venus
utc_venus_year =  convertStringsToChars(string(datetime(utc_venus_dep) + days(224.7), "yyyy-MM-dd'T'HH:mm:ss"));
[venus_spice_states, et_venus_spice] = spice_states(utc_venus_dep, utc_venus_year, max_pts, sun_ID, venus_ID);
venus_spice_r = venus_spice_states(:,1:3);

%Jupiter
utc_jupiter_year =  convertStringsToChars(string(datetime(utc_ida_arr) + days(4332.59), "yyyy-MM-dd'T'HH:mm:ss"));
[jupiter_spice_states, et_jupiter_spice] = spice_states(utc_ida_arr, utc_jupiter_year, max_pts, sun_ID, '5');
jupiter_spice_r = jupiter_spice_states(:,1:3);

%Plotting the Model Trajectory and Spice trajectory for Galileo
figure(19);
plot(gal_dt,gal_model_v,'-b');
hold on;
plot(gal_spice_dt,gal_spice_Vmag,'--r');
plot([datetime(utc_venus_arr) datetime(utc_venus_arr)], [0 45], '-.', 'Color', '[0.7, 0, 1]');
plot([datetime(utc_earth1_arr) datetime(utc_earth1_arr)], [0 45], '-.', 'Color', '[0.1, 0.8, 1]');
plot([datetime(utc_earth2_arr) datetime(utc_earth2_arr)], [0 45], '-.', 'Color', '[0.1, 0.8, 1]');
plot([datetime(utc_JSOI_gal) datetime(utc_JSOI_gal)], [0 45], '-.', 'Color', '[0.9, 0.6, 0]');
hold off;
grid on;
xlabel('UTC');
ylabel('V_{mag} (km/s)');
title('Model vs Spice Velocities Throughout Galileo Mission Interplanetary Cruise');
legend('Galileo - Model','Galileo - Spice', 'Venus Flyby', 'Earth-1 Flyby', 'Earth-2 Flyby', 'Jupiter SOI Entry');

%Plotting the Model Trajectory and Spice trajectory for Galileo
figure(20);
plot3(gal_model_x, gal_model_y, gal_model_z, '-b');
hold on
plot3(gal_spice_r(:,1), gal_spice_r(:,2), gal_spice_r(:,3), '--r');

plot3(0, 0, 0, '.','Color', '[1, 0.8, 0.1]','MarkerSize',marker*3);
plot3(venus_dep_state(1), venus_dep_state(2), venus_dep_state(3), '.','Color', '[0.7, 0, 1]','MarkerSize', marker);
plot3(gaspra_dep_state(1), gaspra_dep_state(2), gaspra_dep_state(3),'.','Color', '[1, 0.4, 0.8]','MarkerSize',marker*0.75);
plot3(earth_dep_state(1), earth_dep_state(2), earth_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(earth1_dep_state(1), earth1_dep_state(2), earth1_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(earth2_dep_state(1), earth2_dep_state(2), earth2_dep_state(3), '.','Color', '[0.1, 0.8, 1]','MarkerSize', marker);
plot3(ida_dep_state(1), ida_dep_state(2), ida_dep_state(3),'.','Color', '[0.4, 0.9, 0.2]','MarkerSize',marker*0.75);
plot3(jupiter_arr_state(1), jupiter_arr_state(2), jupiter_arr_state(3), '.','Color', '[0.9, 0.6, 0]','MarkerSize', marker*2);

plot3(earth_spice_r(:,1), earth_spice_r(:,2), earth_spice_r(:,3), '--','Color', '[0, 0.9, 0.8]');
plot3(venus_spice_r(:,1), venus_spice_r(:,2), venus_spice_r(:,3), '--m');

hold off;
view(3);
grid on;
axis equal;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Model vs Spice Trajectory of Galileo Mission Interplanetary Cruise');
legend('Galileo - Model','Galileo - Spice','Sun', 'Venus','Gaspra', 'Earth','','','Ida','Jupiter');

%**************************************************************************
%                                   CLEAN
%                                     UP
%**************************************************************************

%Unload kernels
cspice_kclear;