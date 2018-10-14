clear all
    load('TEdelays00_rand_all.mat')
    TEdelays00_rand_all(:,1,:)=0; 
    [A B C] = size(TEdelays00_rand_all);
     for ww = 1:C
        for ii = 2:A
            temp = circshift(TEdelays00_rand_all(ii,1:ii,ww),[0 -1]);
            TEdelays00_rand_all(ii,1:ii,ww) = temp;
%             clear temp
        end
     end
  save(['TEdelays00_rand_all.mat'],'TEdelays00_rand_all')

