clear;
global N N_i N_e q_e q_i tau_e tau_i rp_i rp_e v_l r_m i_bg v_r dt alpha beta lowFR noiseVar;
f=figure
set(f, 'Position', [100, 100, 1400, 500]);
%set(f, 'color', [0 0 0]);
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 70 70];
vid = VideoWriter('~/SmoothMANAVid.avi', 'Motion JPEG AVI');
open(vid);
f.InvertHardcopy = 'off';

%% Main Sim constants
N = 500;
N_i = 100;
N_e = 400;
dt = 0.5; %ms - integration time step (euler)
time_f = 5E5; %ms -- simulation duration
onsetTime = 20000;
%% LIF Neuron parameters and vars
i_e = zeros(N,1);
i_i = zeros(N,1);
th = -50 .* ones(N,1);%mV
lastSpk = zeros(N,1);
i_bg = 18; % nA
v_l = -70; %mV
v_r = -55; %mV
v_m = zeros(N,1)+v_r;
rp_e = 3;%ms
rp_i = 2;%ms
tau_e = 30; %ms
tau_i = 20; %ms
q_e = 3; %ms
q_i = 6; %ms
r_m = 1.0; %M. ohm
noiseVar = 0.01;

%% Spike recording
asdf = cell(N+2,1);
spks_tmp = zeros(N, 10000, 'uint8');
spkInd = 1;
spks = zeros(N,1);

%% MHP and HP stuff
ra_tau = 0.001 * dt;
prefFRs = zeros(N,1);
lambda = 5E-5;
eta = 0.0001;
eta_f = 0.0001;
eta_decay = 1E-5;
alpha = 2.5;
lowFR = 1.5; %Hz
beta = 30;

%% Weight Matrices & input
wtMn = 1;
wtStd = 1;
adj_mat = (rand(N) < 0.1);
wtMat =  0.1 .* adj_mat .* abs((randn(N).*wtStd+wtMn));
%dlyMat = adj_mat.* 
wtE = sparse(wtMat(:, 1:N_e));
wtI = sparse(wtMat(:, N_e+1:end));
adj_mat = double(adj_mat)';
adj_mat(adj_mat(:)==0) = NaN;
inD = sum(adj_mat, 'omitnan')';
noInp = 100;
inpMat = 2.*sparse((rand(N, noInp) < 0.2) .* (abs(randn(N, noInp).*wtStd+wtMn)));
inpMnFR = 10;
ef= zeros(N,1);
estFR = ones(N,1).*0.001;
pfRecord = zeros(N, 10000);
%% Run 
for time = 0:dt:time_f
    
    % post synaptic potential sum decay (using conv.) 
    i_e = i_e - dt*(i_e./q_e);
    i_i = i_i - dt*(i_i./q_i);
    
    inp = rand(100, 1) < (inpMnFR/1000 * dt); % Generate Poisson inputs
    
    % Accumulate inputs e/i
    i_e = i_e + inpMat * inp;
    i_e = i_e + wtE * spks(1:N_e);
    i_i = i_i + wtI * spks(N_e+1:end);
    ref = time - lastSpk;
    ref(1:N_e) = ref(1:N_e) <= rp_e;
    ref(N_e+1:end) = ref(N_e+1:end) <= rp_i; 
    ref = ref==1;
    i_e(ref) = 0;
    i_i(ref) = 0;
    
    % LIF Neuron volts
    [v_m, lastSpk, spks] = LIF(v_m, th, lastSpk, i_e, -i_i, ref, time);
%     if (time > onsetTime)
        ef = ef + spks;
        tauA = 10000./sqrt(estFR.*1000);
        ef = ef - dt .* ef./tauA;
        estFR = estFR + dt .* (ef./tauA - estFR);
%     else
%         rec = rec + spks;
%         estFR = rec./time;
%         ef = estFR;
%     end

    
    % Cache spikes
    spks_tmp(:, spkInd) = spks;
    spkInd = spkInd+1;
    % Dump to asdf every 10000 iterations (minimizes array resizes)
    if(spkInd > 10000)
        spkInd = 1;
        for ii=1:N
            asdf{ii} = [asdf{ii}, find(spks_tmp(ii,:))*dt + time-(10000*dt)];
        end
    end
    
    if (time > onsetTime)
        if mod(time,200) == 0
            subplot(122);
            cla;
            hold on;
            title('Target Firing Rate Distribution');
            plot(-3:0.01:5, pdf(orig, -3:0.01:5), 'k', 'LineWidth', 2);
            kd = fitdist(log(prefFRs), 'Normal');
            plot(-3:0.01:5, pdf(kd, -3:0.01:5), '--k', 'LineWidth', 2);
            kd = fitdist(log(prefFRs), 'Kernel');
            plot(-3:0.01:5, pdf(kd, -3:0.01:5), 'LineWidth', 2);
            kd = fitdist(log(1000*estFR), 'Kernel');
            plot(-3:0.01:5, pdf(kd, -3:0.01:5), 'Color', [.05 .55 .17], 'LineWidth', 2);
            %overlayedHistogram(orig, log(prefFRs), log(estFR*1000), 'Normalization', 'pdf');
            legend( 'Original EFRs', 'Lognormal Fit', 'Target FRs', 'Estimated FRs');
            xlim([-3 5]);
            ylim([0 .6]);
            xlabel('Ln(firing rates)');
            hold off;
        end
        prefFRs = metaHP(adj_mat, inD, prefFRs, 1000 .* estFR, eta, time);
        if mod(time,200) == 0
            writeVideo(vid,getframe(gcf));
            cla;
        end
 
        
        if (uint32(mod(time, 10))==0)
            
            pfRecord(:, uint32((time-onsetTime)/10)+1) = prefFRs;
            
        end
        prefFRs(prefFRs < 0.01) = 0.01;
        prefFRs(prefFRs > 100) = 100;
        th = homPlast(th, 1000.* estFR, prefFRs, lambda);
        eta = eta + dt * (eta_f-eta) * eta_decay;
        
        
    else
        prefFRs = estFR * 1000;
        prefFRs(prefFRs < 0.01) = 0.01;
        prefFRs(prefFRs>100) = 100;
        if time == onsetTime
            orig = fitdist(log(prefFRs), 'kernel');
        end
    end
    
    if (mod(uint32(time/dt), uint32(1000/dt))==0)
        disp(time);
    end
    
    
end

close(vid);




