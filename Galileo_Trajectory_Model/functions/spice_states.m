function [states, et_vec] = spice_states(start_utc,stop_utc,num_states, obs, targ)
    start_et = cspice_str2et(start_utc);
    stop_et = cspice_str2et(stop_utc);

    if num_states > 1
        time_interval = (stop_et-start_et)/(num_states-1.);
        et_vec = (0:(num_states-1))*time_interval+start_et;
    else
        et_vec = start_et;
    end
    
    [states_t, ~] = cspice_spkezr( targ, et_vec, 'J2000', 'LT', obs );

    states = transpose(states_t);

end