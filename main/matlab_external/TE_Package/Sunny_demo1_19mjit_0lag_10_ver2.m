clear all; close all;

addpath ./TEpackage_old;
addpath ./TEpackage
tic
for data_num10 =  [10];
    
    save ./data_num10;
    clear all
    load ./data_num10;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if data_num10 ==  10;
        data_name = '2013-01-09-0';
        load ./data000/2013-01-09-0/asdf.mat
        %load ./data000/2013-01-10-1/xy.mat;
        %clear location;    location = [x' y'];
    end
    
    % ================================================-
    % Change the bin-size.
    % ================================================-
    asdf_new = ChangeBinning(asdf_raw, 1); % 2);
    matlabpool open
    % ================================================-
    % Calculation of TE for each delay
    % ================================================-
    delay0 = [0];       delay1 = delay0+1;
    clear TEdelays00 TEstds00;
    %[TEdelays00, TEstds00] = ASDFTEdelays(asdf_new, delay0);
     [TEdelays00] = ASDFTE_parallel(asdf_new, delay0);
     str = ['save ./data/',num2str(data_name),'_graph/TEdelays00 TEdelays00;'];               eval(str);


    clear TEdelays00_rand_all TEdelays000_rand
    % ================================================-
    %  ISI shuffling (100 times) for testing the signiifice
    %  of TE values for each each neuron pair.
    % ================================================-
    shuffle_num_all = 20;
    
%    TE_rand_all  = zeros(size(TEdelays00,1), size(TEdelays00,2), shuffle_num_all);
%     TEdelays00_rand_all   = TEdelays00;
    
    %TEstds00_rand   = TEstds00;
    
%    for kk = 1:size(asdf_new,1)
%        for jj = 1:size(asdf_new{kk},2)-1
%            asdf_int{kk}(jj) = asdf_new{kk}(jj+1)-asdf_new{kk}(jj);
%        end
%    end
    
    clear TEdelays00_rand_all TEstds00_rand_all ;
    
       n_neu = size(asdf_new,1)-2; time_ax = length(delay0);
       TEdelays00_rand_all = zeros(n_neu,n_neu,shuffle_num_all);
      for shuffle_num = 1:shuffle_num_all
        %disp('Shuffling#: ',num2str(shuffle_num));
        clear asdf_rand;
        
         asdf_rand = asdf_new;
        for kk = 1:size(asdf_new,1)-2 
              
            for jj = 1:size(asdf_new{kk},2)-5
                jitter_bin0 = randperm(19) - 10;
                %jitter_bin0 = randperm(11) - 6;
                
                jitter_bin = jitter_bin0(1);
                if asdf_rand{kk}(jj) >= 15
                    asdf_rand{kk}(jj) = asdf_new{kk}(jj) + jitter_bin;
                    asdf_rand{kk}(jj+1) = asdf_new{kk}(jj+1) - jitter_bin;
                end
            end
            asdf_mod{1} = asdf_rand{kk}; ind = setdiff(1:size(asdf_new,1)-2,kk);
            temp = ASDFSubsample(asdf_new,ind); 
            temp{end}(1,1) = temp{end}(1,1) +1;
             asdf_mod = [asdf_mod; temp]; %asdf_mod{end}(1,1)
%             size(asdf_mod)
        
        
        %clear TEdelays00_rand TEstds00_rand;
        %[TEdelays00_rand, TEstds00_rand] = ASDFTEdelays(asdf_rand, delay0);
         [TEdelays00_rand] = ASDFTE_parallel_mod_0delay(asdf_mod, delay0);
          clear asdf_mod
         TEdelays00_rand_all(kk,:,shuffle_num) = TEdelays00_rand;
%          TEdelays00_rand_all(kk,:,:) = TEdelays00_rand(:,:,:);
%           TEdelays00_rand_all = vertcat(TEdelays00_rand_all,TEdelays00_rand);
        
        end
        %TEstds00_rand_all(:,:,:,shuffle_num) = TEstds00_rand(:,:,:);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str = ['mkdir ./data/',num2str(data_name),'_graph/;'];   eval(str);
 %   str = ['save ./data/',num2str(data_name),'_graph/TEdelays00 TEdelays00;'];               eval(str);
    str = ['save ./data/',num2str(data_name),'_graph/TEdelays00_rand_all TEdelays00_rand_all;'];  eval(str);
%     str = ['save ./data/',num2str(data_name),'_graph/TEdelays00_rand TEdelays00_rand;'];  eval(str);
    %str = ['save ./data/th8_19msjit/',num2str(data_name),'_graph/TEstds00_rand_all TEstds00_rand_all;'];      eval(str);
    
    close all;
    
% end
end
toc

matlabpool close

