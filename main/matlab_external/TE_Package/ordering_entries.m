clear all
for qq = 1:100
   
    load(['TEdelays_rand_all', num2str(qq),'.mat'])
    TEdelays_rand_all(:,1,:)=0; 
    [A B C] = size(TEdelays_rand_all);
     for ww = 1:C
        for ii = 2:A
            temp = circshift(TEdelays_rand_all(ii,1:ii,ww),[0 -1]);
            TEdelays_rand_all(ii,1:ii,ww) = temp;
            clear temp
        end
     end
     save(['TEdelays_rand_all',num2str(qq),'.mat'],'TEdelays_rand_all')
   
end
