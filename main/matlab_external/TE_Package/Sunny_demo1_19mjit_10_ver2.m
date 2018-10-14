clear all; close all;
addpath ./TEpackage_old;
addpath ./TEpackage

for data_num10 = [10];
    data_num10
    save ./data_num10;
    clear all
    load ./data_num10;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if data_num10 ==  10;
        data_name = '2013-01-09-0';
        load ./data000/2013-01-09-0/asdf.mat
        %load ./data000/2013-01-10-3/xy.mat;
        %clear location;    location = [x' y'];
    end
     matlabpool
    tic
    % ================================================
    % Calculation of distances between neurons
    % ================================================
    %clear distance
    %for i = 1:size(location,1)
        %for j = 1:size(location,1)
            %distance(i,j) = ...
                %sqrt((location(i,1)-location(j,1)).*(location(i,1)-location(j,1)) + ...
                %(location(i,2)-location(j,2)).*(location(i,2)-location(j,2)));
        %end
    %end
    
    % ================================================-
    % Change the bin-size.
    % ================================================-
    asdf_new = ChangeBinning(asdf_raw, 1); % 2);
    
    % ================================================-
    % Calculation of TE for each delay
    % ================================================-
    delay0 = [1:30];       delay1 = delay0+1;
     %delay0 = 1;
    %clear TEdelays TEstds;
    %[TEdelays, TEstds] = ASDFTEdelays(asdf_new, delay0);
%      matlabpool open
    [TEdelays] = ASDFTE_parallel(asdf_new, delay0);

%     [TEdelays, TEstds] = ASDFTE(asdf_new, delay0);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str = ['mkdir ./data/',num2str(data_name),'_graph/;'];   eval(str);
 %   %str = ['save ./data/',num2str(data_name),'_graph/distance distance;'];   eval(str);
    str = ['save ./data/',num2str(data_name),'_graph/TEdelays TEdelays;'];   eval(str);

    clear TEdelays_rand_all TEdelays0_rand
    % ================================================-
    %  ISI shuffling (100 times) for testing the signiifice
    %  of TE values for each each neuron pair.
    % ================================================-
     shuffle_num_all = 10;
% %    TEdelays = zeros(asdf_new{end}(1,1),asdf_new{end}(1,1));
% %     TE_rand_all  = zeros(size(TEdelays,1), size(TEdelays,2), shuffle_num_all);
% %     TEdelays_rand_all   = TEdelays;
%     %TEstds_rand   = TEstds;
%     
% %     for kk = 1:size(asdf_new,1)
% %         for jj = 1:size(asdf_new{kk},2)-1
% %             asdf_int{kk}(jj) = asdf_new{kk}(jj+1)-asdf_new{kk}(jj);
% %         end
% %     end
%    c
     clear TEdelays_rand_all TEstds_rand_all ;
     n_neu = asdf_new{end}(1,1); time_ax = length(delay0);
     for shuffle_num = 1:shuffle_num_all
         TEdelays_rand_all = zeros(n_neu,n_neu,time_ax);
         clear asdf_rand 
         asdf_rand = asdf_new;
         for kk = 1:size(asdf_new,1)-2
             for jj = 1:size(asdf_new{kk},2)-5
                  jitter_bin0 = randperm(19) - 10;
 %                 jitter_bin0 = randperm(11) - 6;
                 jitter_bin = jitter_bin0(1);
                 if asdf_rand{kk}(jj) >= 2
                     asdf_rand{kk}(jj) = asdf_new{kk}(jj) + jitter_bin;
                     asdf_rand{kk}(jj+1) = asdf_new{kk}(jj+1) - jitter_bin;
                 end
             end
             asdf_rand{kk} = sort(asdf_rand{kk});
             asdf_mod{1} = asdf_rand{kk}; ind = setdiff(1:size(asdf_new,1)-2,kk);
             temp = ASDFSubsample(asdf_new,ind); 
             temp{end}(1,1) = temp{end}(1,1) +1;
              asdf_mod = [asdf_mod; temp]; %asdf_mod{end}(1,1)
% %             size(asdf_mod)
%        
% %         clear TEdelays_rand TEstds_rand;
%         %[TEdelays_rand, TEstds_rand] = ASDFTEdelays(asdf_rand, delay0);
         [TEdelays_rand] = ASDFTE_parallel_mod(asdf_mod, delay0);
         clear asdf_mod
           TEdelays_rand_all(kk,:,:) = TEdelays_rand;
         clear TEdelays_rand
         end
% 
     str = ['save ./data/',num2str(data_name),'_graph/TEdelays_rand_all',num2str(shuffle_num),' TEdelays_rand_all;']; eval(str);
%     %str = ['save ./data/',num2str(data_name),'_graph/TEstds_rand_all',num2str(shuffle_num),' TEstds_rand;'];     eval(str);
%     
     end
    clear TEdelays_rand_all 
%     % ================================================-
% %     figure(1);
% %     subplot(2,1,1); hist(nonzeros(TEdelays(:,:,:)),100);
% %     subplot(2,1,2); hist(nonzeros(TEdelays_rand(:,:,:)),100);
% %     
% %     clear TEpeak0 TEdelays0
% %     
% % 
% %     close all;
%     
end
toc
 matlabpool close


