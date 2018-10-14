% load all shuffled data

clear all
for ii = 1:100
    load(['TEdelays_rand_all',num2str(ii)]);
    
    TEdelays_rand(:,:,:,ii)=TEdelays_rand_all;
    clear TEdelays_rand_all
end
TEdelays_rand_all = TEdelays_rand; clear TEdelays_rand;
save('TEdelays_rand_all.mat','TEdelays_rand_all','-v7.3')



%clear all; te = [];
%for ii = 1:4
%    load(['TEdelays00_rand_all',num2str(ii)]);
%    te =cat(4,te,TEdelays00_rand_all);
%    clear TEdelays_rand
%end
%TEdelays00_rand_all = te; clear te
%save('TEdelays00_rand_all.mat','TEdelays00_rand_all')

