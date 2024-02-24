function [states_ca] = slice_flyby(states_fb, etvec_fb, et_pre_ca)
    et_min = abs(etvec_fb(1));
    index = 1;

    for i = 1:length(etvec_fb)
        et_diff = abs(etvec_fb(i) - et_pre_ca);
        if et_diff < et_min
            et_min = et_diff;
            index = i;
        end
    end

    states_ca = states_fb(index:end, :);
end

