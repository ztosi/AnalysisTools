% delays analysis. Look at the delays of 
% triplet connections which could be due to 
% common drive

clear all; load PDF_1_16_30ms.mat
load wgts_1_16ms.mat; [q1 q2 q3]= size(PDF);iti=500;
pdf_cor = zeros(q1,q2,iti);
clearvars -except wgt mte PDF max_ind q1 cnd iti
thr = 35;
tic
parfor nn = 1:iti
    mte = PDF(:,:,thr).*wgt;[A B] = size(mte);
    bin = zeros(A,B);
    bin(find(mte~=0))=1;delay = max_ind.*bin;
    % randomly pick nodes to start
    ind = datasample(1:A,A,'replace',false);
    out_r(nn) = sum(mte(ind(1),:));
    for ii = 1:length(ind)
        tar = find(mte(ind(ii),:)~=0); 
        if length(tar)>1
            actdel = delay(tar,tar);
            del_sr_tar = delay(ind(ii),tar);L1 = length(del_sr_tar);
            cmn_drv_del = zeros(L1,L1);
            for mm = 1:L1
                cmn_drv_del(mm,:)= del_sr_tar-del_sr_tar(mm);
            end
            cmn_drv_del(cmn_drv_del<0)=999;
            cmn_drv_del(cmn_drv_del==0)=999;
            % check to see if there are delays satisfying
        % common drive in the actual delay matrix
            diff_d = cmn_drv_del - actdel;
            %[I J] = find(diff_d ==0);
            [I J] = find(diff_d <1 | diff_d >-1);
            if length(I)>0
                sr_tar = [tar(I);tar(J)];
                [M N] = size(sr_tar);
                for qq = 1:N
                    % eliminate common drive
                    mte(sr_tar(1,qq),sr_tar(2,qq))=0;
                    % eliminate chain
                    %mte(ind(ii),sr_tar(2,qq))=0;
                    %eliminate strongest of the 2 spurious
%                     tt1 = mte(ii,sr_tar(2,qq));
%                     tt2 = mte(sr_tar(1,qq),sr_tar(2,qq));
%                     if tt1<tt2
%                         mte(ii,sr_tar(2,qq))=0;
%                     else
%                         mte(sr_tar(1,qq),sr_tar(2,qq))=0;
%                     end
                end
            end
        end

    end
    temp = zeros(A,B);temp(mte~=0)=1;
    pdf_cor(:,:,nn) = temp; %clear temp
    cnd_aft(nn) = (length(find(mte~=0))/(A*(A-1)))*100;   
    %clear mte
end
toc

pdf = mean(pdf_cor,3);
%subplot(1,2,1); imagesc(pdf); axis square;
pdf(find(pdf~=1))=0;

%subplot(1,2,2); imagesc(PDF(:,:,45)-pdf); axis square
per_chn = ((cnd(thr)-mean(cnd_aft))/cnd(thr))*100
cnd_act = cnd(thr);err = std(cnd_aft);
cd_cor = mean(cnd_aft);
[cnd_act cd_cor err]

save('1.mat','cnd_act','cd_cor','err','pdf')
% [A B]= hist(cnd_aft,100)
% bar(B,A./sum(A)); hold on
% plot([cnd(45) cnd(45)],[0 0.5])
% plot([0.01:0.01:0.8],cnd); hold on
% plot([0.01:0.01:0.8],cnd_aft,'r');
%save('pdf_cor.mat','pdf_cor')