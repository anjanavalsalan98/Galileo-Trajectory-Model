function [osc] = pv2osc(r,v,mu)
    eps = 10.e-8;
    R_mag = norm(r);
    V_mag = norm(v);
    h = cross(r,v);
    h_mag = norm(h);
    ecc = (1/mu)*((V_mag^2 - mu/R_mag)*r-(dot(r,v)*v));
    ecc_mag = norm(ecc);
    epsilon = ((V_mag^2)/2)-(mu/R_mag);
    inc = acos(h(3)/h_mag);
    n = cross([0; 0; 1],h);
    n_mag = norm(n);
    if abs(epsilon) > eps
        a = -mu/(2*epsilon);
    else
        a = Inf;
    end
    if (ecc_mag > eps) && (inc > eps)
        Omega = acos(n(1)/n_mag);
        if n(2)<0
            Omega = 2*pi-Omega;
        end
        w = acos(dot(n,ecc)/(n_mag*ecc_mag));
        if ecc(3)<0
            w = 2*pi-w;
        end
        f = acos(dot(ecc,r)/(ecc_mag*R_mag));
        if (dot(v,r)/R_mag) < 0
            f = 2*pi-f;
        end
    elseif (ecc_mag < eps) && (inc<eps)
        Omega = 0;
        w = 0;
        f = acos(r(1)/R_mag);
        if r(3)<= 0
            f = 2*pi-f;
        end
    
    elseif (ecc_mag > eps) && (inc < eps)
        Omega = 0;
        w = acos(ecc(1)/ecc_mag);
        if ecc(2)<0
            w = 2*pi-w;
        end
        f = acos(dot(ecc,r)/(ecc_mag*R_mag));
        if dot(ecc,r)<0
            f = 2*pi-f;
        end
    elseif (ecc_mag<eps) && (inc>eps)
        Omega = acos(n(1)/n_mag);
        if n(2)<0
            Omega = 2*pi-Omega;
        end
        w = 0;
        f = acos(dot(n,r)/(R_mag*n_mag));
        if r(3)<=0
            f = 2*pi-f;
        end
    end
   
    osc = [a, ecc_mag, inc, Omega, w, f, h_mag, epsilon];
end
