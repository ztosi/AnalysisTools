%% Find Rich/non-Rich neurons in all avalanches of particular lengths

% Detects rich/non-rich neurons in sequences of all avalanches to
% see if there is a repeating profile in the sequence of nodes.
%
% Hadi Hafizi, Sept. 2015

%clear all

%dataset = '02-0';
BinSize = 1;

%cd(dataset)
% load asdf.mat;  
%load PDF_1_16_30ms.mat; load wgts_1_16ms.mat
%W = PDF(:,:,45).*wgt;oute = sum(W,2);
[Rich Rich_idx] = sort(kOut,'descend');
num_rich = ceil(0.2*length(Rich));
Rich_idx = Rich_idx(1:num_rich); nonRich_idx = setdiff(1:length(Rich),Rich_idx);
% inte=sum(W,1); % do the same thing for indegree
% [inRichness, inRichness_idx] = sort(inte,'descend');
% innum_rich = ceil(0.2*length(inRichness));
% inRich_idx = Richness_idx(1:innum_rich); innonRich_idx = setdiff(1:length(inRichness),inRich_idx);
% nr_inte = inte(innonRich_idx);
% [innr_vals, innr_idx] = sort(nr_inte,'descend');

%load(['aval_list',dataset,'_bs',num2str(BinSize),'ms.mat'])

% choose particular avalnche lengths and plot the distribution of rich
% neurons along the avalanche sequence

len_chosen = 500;%5:5:70;% sze_chosen = [70];
Rich_sum_all = zeros(max(len_chosen),length(len_chosen));
NRich_sum_all = zeros(max(len_chosen),length(len_chosen));
for i=1:length(str_aval)
    leni(i)=length(unique(str_aval{i}(:,2)));
end

for ChosenLenNum=1:length(len_chosen)
    SampChosenLenNum = 1;
    for AllAvalNum=1:length(str_aval)
        if leni(AllAvalNum)>=len_chosen(ChosenLenNum)
            AvalChosenLen{SampChosenLenNum} = str_aval{AllAvalNum};
            %         AvalRich_ChosenLen = [str_aval{ii}];
            SampChosenLenNum = SampChosenLenNum + 1;
        end
    end
    
    if SampChosenLenNum==1
        warning('No avalanche of length %s',num2str(len_chosen(ChosenLenNum)))
    else
    
    Rich_Propor = zeros(len_chosen(ChosenLenNum),1);
    NRich_sum = zeros(len_chosen(ChosenLenNum),1);
    for SampChosenLenNum=1:length(AvalChosenLen)
        AvalTemp = AvalChosenLen{SampChosenLenNum};
        AvalTempNR = [AvalChosenLen{SampChosenLenNum}];
% % % % % % %         [crap iAvalTemp ib] = intersect(AvalTemp(:,1),
% Rich_idx); %%%%%%% "intersect" does not count repetition and this screws
% up what I'm trying to do. I should put this line inside the time loop.
%         AvalTemp(iAvalTemp,1) = 1;
%         AvalTemp(setdiff(1:size(AvalTemp,1),iAvalTemp),1) = 0;
%         AvalTempNR(iAvalTemp,1) = 0;
%         AvalTempNR(setdiff(1:size(AvalTemp,1),iAvalTemp),1) = 1;
        for tAval=1:len_chosen(ChosenLenNum)
            Rich_Propor(tAval) = length(intersect(AvalTemp(find(AvalTemp(:,2)==tAval),1),...
                Rich_idx))/length(find(AvalTemp(:,2)==tAval));
%             Rich_sum(tAval) = sum(AvalTemp(AvalTemp(:,2)==tAval),1);%/length(AvalTemp(AvalTemp(:,2)==ww));
%             NRich_sum(tAval) = sum(AvalTempNR(AvalTempNR(:,2)==tAval),1);%/length(AvalTempNR(AvalTempNR(:,2)==ww));
        end
%         Rich_sum_all(:,ChosenLenNum) = Rich_sum_all(:,ChosenLenNum) + padarray(zscore(Rich_sum),[max(len_chosen)-length(Rich_sum) 0],'post');
%         NRich_sum_all(:,ChosenLenNum) = NRich_sum_all(:,ChosenLenNum) + padarray(zscore(NRich_sum),[max(len_chosen)-length(NRich_sum) 0],'post');
        Rich_sum_all(:,ChosenLenNum) = Rich_Propor/length(AvalChosenLen) + Rich_sum_all(:,ChosenLenNum);
        clear AvalTemp AvalTempNR iAvalTemp
    end
%     Rich_sum_all(:,ChosenLenNum) = Rich_sum_all(:,ChosenLenNum)/length(AvalChosenLen);
%     Rich_mean_all(:,ChosenLenNum) = mean(Rich_sum_all,2);
%     NRich_sum_all(:,ChosenLenNum) = NRich_sum_all(:,ChosenLenNum)/length(AvalChosenLen);
    clear AvalChosenLen
    end
end

figure;
plot(Rich_sum_all);

cd ..

%%
% 
% 
% figure;
% surf(1:max(len_chosen),len_chosen,Rich_sum_all'); colormap('jet');colorbar; shading('interp')
% title('Rich Node Profile through Avalanches (Rich)','FontSize',16);
% xlabel(['Time bin [',num2str(BinSize),'ms]'],'FontSize',16);
% ylabel('Avalanche Length','FontSize',16);
% zlabel([{['z-score of Rich Node Counts']} {['Normalized by the Number of Samples']}],'FontSize',16);
% set(gca,'FontSize',16);
% view(2)
% cd Figures
% savefig(['AvalRichProfile',num2str(BinSize),'msAllSizes.fig']);
% cd ..
% 
% figure;
% surf(1:max(len_chosen),len_chosen,NRich_sum_all'); colormap('jet');colorbar; shading('interp')
% title('Rich Node Profile through Avalanches (nonRich)','FontSize',16);
% xlabel(['Time bin [',num2str(BinSize),'ms]'],'FontSize',16);
% ylabel('Avalanche Length','FontSize',16);
% zlabel([{['z-score of Rich Node Counts']} {['Normalized by the Number of Samples']}],'FontSize',16);
% set(gca,'FontSize',16);
% view(2)
% cd Figures
% savefig(['AvalnonRichProfile',num2str(BinSize),'msAllSizes.fig']);
% cd ..
% 
% 
% cd ..
