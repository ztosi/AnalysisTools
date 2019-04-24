%lognormalConsistency
%want to see if Zoe Tosi's rule, by itself, leads to a self-consistent solution
%it may not, suggesting that some other rule is needed. And/or it may
%suggest that only certain types of connecetivity will give stable
%solutions

%set size of simulation
N = 800;   %number of neurons
FR = rand(N,1); %vector giving random [0 - 1] initial firing rates of these neurons

%set connectivity 
%make it random at first, to be possibly changed later
C = zeros(N, N);
%pick sparsity from 0 (no connections) to 1 (all to all connections)
sparsity = 0.1;
%populate it with connections at this sparsity
Mask = rand(N, N); 
I = find(Mask < sparsity);
C(I) = 1;
%figure; spy2(C); grid off

%Now set up self-consistency of Zoe's rule
alpha = 0.01;  %learning rate parameter to prevent oscillations during learning, provides stability if small
beta = 0.2;   %slope of the logistic, or squashing, function
gamma = 100; %sets maximum FR allowed
iterations = 110; %number of time steps
FRall = zeros(N, iterations);
FRall (:, 1) = FR;
FreezeTime = floor(iterations/2);

for i=1:iterations-1  %go through each time step
    
    %just to see what happens, reset the connectivity matrix halfway
    %through
   % if i == floor(iterations/2)
        % C = ones(N, N);
%         %populate it with connections at this sparsity
% %         Mask = rand(N, N);
% %         I = find(Mask < sparsity);
% %         C(I) = 1;
     %end
    for j= 1:N %go through each neuron
        %get average input to neuron j
        avg = (1/N) * (C(:, j)' * FRall(:, i));
        
        %update its firing rate according to the rule, but do not allow it
        %to go below zero (all FR must be positive)
        holder = FRall(j, i) + alpha*(FRall(j, i) - avg)*(FRall(j, i)); % + rand; %without adding rand, you do not get lognormal distribution
       % FRall(j, i+1) = holder;
        if holder > 0
            FRall(j, i+1) = holder; %gamma * (1 + exp(-beta * (holder - 0.5)))^(-1);   %was = holder;
        else
            FRall(j, i+1) = rand;
        end
        
        %try FR homeostasis to see if this stabilizes/ turn this off if you want to continue dispersing firing rates 
%         if i >= floor(iterations/2)
%             FRall(j, i+1) = FRall(j, i+1) + 1.75*alpha*(FRall(j, FreezeTime) - FRall(j, i)); % + rand;
%         end
 
    end
end
            
%playing with a logistic equation to squash firing rates
% for x=1:100
%     y(x) = gamma * (1 + exp(-beta * (x - 50)))^(-1);  
% end
% figure; plot(y);

 
%plot firing rates over time
timepoint = iterations; %100; 
figure; plot(FRall(1, 1:timepoint)); hold on
for k=1:N
   plot(FRall(k, 1:timepoint)); 
   hold on 
end

%plot resulting firing rate distribution
timepoint = iterations; 
[Num,EDGES] = histcounts(FRall(:, 1:timepoint),'BinMethod', 'sqrt');
figure; histogram(FRall(:, 1:timepoint));

timepoint = iterations;
[Num,EDGES] = histcounts(log10(FRall(:, timepoint)));
figure; histogram(log10(FRall(:, timepoint)));











