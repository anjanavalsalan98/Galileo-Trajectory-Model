function [JD] = julian_date(utc)
    date = datetime(utc);
    [Y, M, D] = ymd(date);
    [h, m ,s] = hms(date);

    F = h/24 + m/(24*60) + s/(24*60*60);
    JD0 = 367*Y-floor(7/4 * (Y + floor((M+9)/12)))+floor(275*M/9)+D+1721013.5;
    JD = JD0 + F;
end

