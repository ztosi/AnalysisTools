% load all shuffled data

clear all
load TEdelays_rand_all

TE_rand_avg = mean(TEdelays_rand_all,4);
TE_rand_std = std(TEdelays_rand_all,[],4);
sel = TEdelays_rand_all(:,:,1:16,:);
err_sel = std(sel,[],4);

clearvars -except TE_rand_avg TE_rand_std err_sel

load TEdelays.mat
temp = TEdelays(:,:,1:16);
[mte max_ind] = max(temp,[],3);
[A B] = size(mte);
TE_shuff = zeros(A,B);
temp1 = TE_rand_avg(:,:,1:16);
for aa = 1:A
    for bb = 1:B
                mte_shuf(aa,bb) = temp1(aa,bb,max_ind(aa,bb));
                fin_err(aa,bb) = err_sel(aa,bb,max_ind(aa,bb));
    end 
end

% wgt = zeros(A,B);
wgt = mte - mte_shuf; wgt(wgt<0)=0;
save('wgts_1_16ms.mat','wgt','mte','mte_shuf','fin_err','max_ind')

