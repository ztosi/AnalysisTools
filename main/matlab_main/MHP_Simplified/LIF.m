function [v_m, lastSpk, spks] = LIF(v_m, th, lastSpk, i_e, i_i, ref, time)
global N N_e tau_e tau_i v_l v_r r_m i_bg dt noiseVar;   
       
       % Linear LIF equations....
       dv_m = dt .* ( (v_l-v_m) + r_m.* (i_e + i_i + i_bg + randn(N,1) .* noiseVar));
       dv_m(1:N_e) = dv_m(1:N_e)./tau_e; % exc
       dv_m(N_e+1:end) = dv_m(N_e+1:end) ./tau_i; % inh
       v_m = v_m + dv_m; % integrate
       v_m(ref) = v_r; % keep ref at the reset
       
        % find spikes
       spks = double(v_m >= th); 
       v_m(spks==1) = v_r; %reset voltage for neurons that spiked
       lastSpk(spks==1) = time; % set new last spike times
       
end