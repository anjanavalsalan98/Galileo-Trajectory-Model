function ang = atan3(y, x)
    ang = mod(atan2(y, x), 2*pi);
end
