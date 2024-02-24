function [v1, v2] = lambert(mu, r1, r2, tof, prograde)
    R1_mag = norm(r1);
    R2_mag = norm(r2);

    delta_TA = delta_TA_fn(r1, r2, prograde);
    
    A = A_fn(R1_mag, R2_mag, delta_TA);

    z = z_soln_fn(0, tof, A, R1_mag, R2_mag, mu);

    y = y_fn(z, R1_mag, R2_mag, A);

    [f, g, ~, dg] = fandg_fn(z, R1_mag, R2_mag,y,A, mu);

    v1 = 1/g * (r2 - f * r1);
    v2 = 1/g * (dg * r2 - r1);
end

function [s] = S(z)
    if z > 0
        s = (sqrt(z)-sin(sqrt(z))) / (sqrt(z))^3;
    elseif z < 0 
        s = (sinh(sqrt(-z))-sqrt(-z)) / (sqrt(-z))^3;
    else
        s = 1/6;
    end
end

function [c] = C(z)
    if z > 0
        c = (1-cos(sqrt(z))) / z;
    elseif z < 0
        c = (cosh(sqrt(-z))-1) / (-z);
    else
        c = 1/2;
    end
end

function [delta_TA] = delta_TA_fn(r1, r2, prograde)
    R1_mag = norm(r1);
    R2_mag = norm(r2);
    cang = acos(dot(r1,r2)/(R1_mag * R2_mag));
    r1xr2 = cross(r1,r2);

    if prograde == 1
        if r1xr2(3) >= 0
            out = cang;
        else
            out = 2 * pi - cang;
        end
    else
        if r1xr2(3) < 0
            out = cang;
        else
            out = 2 * pi - cang;
        end
    end
    delta_TA = out;
end

function [A] = A_fn(R1_mag, R2_mag, delta_TA)
    A = sin(delta_TA)*sqrt((R1_mag*R2_mag)/(1 - cos(delta_TA)));
end

function [y] = y_fn(z, R1_mag, R2_mag, A)
    y = R1_mag + R2_mag + A*(((z*S(z))-1)/(sqrt(C(z))));
end

function [f, g, df, dg] = fandg_fn(z, R1_mag, R2_mag,y,A, mu)
    f = 1 - y/R1_mag;
    g = A * sqrt(y/mu);
    df = (sqrt(mu))/(R1_mag * R2_mag) * sqrt(y/C(z)) * (z * S(z) - 1);
    dg = 1 - y/R2_mag;
end

function [F] = F_fn(z, tof, A, y, mu)
    F = ((y/C(z))^1.5) * S(z) + (A * sqrt(y)) - (sqrt(mu) * tof);
end

function [dF] = dF_fn(z, A, y, y0)
    if z ~= 0
        dFi = ((y/C(z))^1.5) * (1/(2*z)*(C(z)-(3/2)*S(z)/C(z))+(3/4)*S(z)*S(z)/C(z));
        dFi = dFi + (A/8)*(3*(S(z)/C(z)) * sqrt(y) + A*(sqrt(C(z)/y)));
    else
        dFi = (sqrt(2)/40)*(y0^1.5) + (A/8)*(sqrt(y0) + A*sqrt(1/(2*y0)));
    end
    dF = dFi;
end

function [z] = z_soln_fn(z0, tof, A, R1_mag, R2_mag, mu)
    z1 = z0;
    for i = 1:100
        y =  y_fn(z1, R1_mag, R2_mag, A);
        y0 = y_fn(0, R1_mag, R2_mag, A);
        ratio = F_fn(z1, tof, A, y, mu)/dF_fn(z1, A, y, y0);
        zi = z1 - ratio;
        z1 = zi;
    end
    z = zi;
end
