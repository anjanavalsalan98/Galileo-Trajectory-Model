function [inc,alt,et] = closest_approach(et_vec,r_states,v_states,planet_radius)
    min_dist = norm(r_states(1, :));

    for i = 1:length(r_states)
        dist = norm(r_states(i, :));
    
        if dist < min_dist 
            min_dist = dist;
            h = cross(r_states(i, :),v_states(i, :));
            inc = acos(h(3)/norm(h));
            et = et_vec(i);
        end
    end
    
    alt = min_dist - planet_radius;
end

