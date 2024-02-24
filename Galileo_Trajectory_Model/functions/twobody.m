function [states, et_vec] = twobody(mu, start_et, stop_et, state_i)
    tspan = linspace(start_et, stop_et, 2500);

    options = odeset('RelTol', 2.5e-14, 'AbsTol',1e-20);
    [et_vec,states] = ode113(@(t,y) eqn_motion_twobody(t, y, mu), tspan, state_i, options);
end

function [dydt] = eqn_motion_twobody(~, y, mu)
    x1 = y(1);
    x2 = y(2);
    x3 = y(3);

    r = norm([x1 x2 x3]);
    mu_r = mu/r^3;

    z1 = y(4);
    z2 = y(5);
    z3 = y(6);

    z4 = - mu_r * x1;
    z5 = - mu_r * x2;
    z6 = - mu_r * x3;

    dydt = [z1; z2; z3; z4; z5; z6];
end