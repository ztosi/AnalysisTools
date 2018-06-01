function [ norm_ran ] = songMotifs( bin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% find frequecy of 3 neuron motifs as in Figure 4A 
% of Song et al 2005

%clear all
%load wgts_1_16ms.mat; load PDF_1_16_30ms.mat;
%bin = PDF(:,:,45);
[A B] = size(bin);

[f F]=motif3struct_bin(bin);
song_order = [3, 1, 2, 4, 6, 9,  5, 7, 8, 10,11,12,13];
f = f(song_order);

clear F
% compare with expected values of motifs
edges = length(find(bin~=0));
%p_nc = length(find(W==0))/(A*(A-1));
num_nc = (A * (B-1)) - edges; % number of 0s in the binary matrix
j = 1;
[C1 C2] = find(bin~=0);

for ww = 1:length(C1)
     if bin(C1(ww),C2(ww))*bin(C2(ww),C1(ww))==0
        j = j +1;
     end
end
num_uni = j-1;
num_bi =  (trace(bin^2)/2) * 2;
% calculate probabilities of uni directional
% and bidirectional connections
p_nc = num_nc/(num_nc+num_uni+num_bi);
p_uni = (0.5*num_uni)/(num_nc+num_uni+num_bi);
p_bi = num_bi/(num_nc+num_uni+num_bi);
% p_uni = (j-1)/(A*(A-1));
% p_bi = (trace(W^2)/2)/((A*(A-1))/2);

% expected counts of 3 neuron motifs
N_exp = exp_count(A,p_nc,p_uni,p_bi,song_order);

norm_ran = f./N_exp';
%save('C:\Users\Nigam\Desktop\15.mat','f','N_exp','norm_ran')
hbars = bar(norm_ran);
set(hbars(1),'BaseValue',1);
%set(gca,'Yscale','log')

end

