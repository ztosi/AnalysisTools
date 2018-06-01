
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program perform (stat. test) & (binalizing) on TE for 2013 cortex data.
% jittering size: 19ms
% Author: Masanori Shimono         2011 originally written for 2010 data 
%         edited  Masanori Shimono 2013
%         cleaned Masanori Shimono 2014
% contact address: nori417@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off;
clear all; close all;
rmpath(genpath('./function/'));
addpath ./TEpackage; addpath ./TEpackage_old;

data_num03 = 1;%[1:15]; 
data_num03
    save ./data_num032 data_num03;
    clear all
    load ./data_num032;

%% %%%=====================================================================
                if data_num03 ==  1; %load ./analysis/asdf2013spon/2013-01-02-0/data000/asdf.mat
                        load ./data000/2013-01-09-0/asdf.mat
                        %load ./data000/2013-01-05-0/xy.mat
                end
                

                %% %%%=================================================================
                if data_num03 == 1;       
                     data_name = '2013-01-09-0';%data_name = '13010200';
                
                end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str = ['load ./data/',num2str(data_name),'_graph/TEdelays ;'];               eval(str);
    str = ['load ./data/',num2str(data_name),'_graph/TEdelays_rand_all ;'];               eval(str);
    str = ['load ./data/',num2str(data_name),'_graph/TEdelays00 ;'];               eval(str);
    str = ['load ./data/',num2str(data_name),'_graph/TEdelays00_rand_all ;'];               eval(str);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %str = ['load ./data/',num2str(data_name),'_graph_CC/cchmap_all;'];               eval(str);

    % ================================================-
    %pp   =  size(cchmap_all,3);
    %pp2 =  (pp+1)/2;

    %cchmap_all(:,:,pp2) = 0 * cchmap_all(:,:,pp2);

    %for i = 1:size(cchmap_all, 1)
   %     peaknum_ac(i) =     max(squeeze(cchmap_all(i,i,:)));
    %end

    %cell_label = ones(1,size(location,1));
    % I don't consider cell categ here, but all connected neurons may be excitatory.
    
    %     cell_label = cell_label12;
    %     cell_label = cell_labelCC;
    % 0: non-classified
    % 1: pyramidal cells
    % 2: interneurons

    %for i =  1: size(cchmap_all,1)
        % refrectory period
                %if peaknum_ac(i) <= 2 | cchmap_all(i,i,pp2-1) >= 0.1*max(cchmap_all(i,i,:)) | cchmap_all(i,i,pp2+1) >= 0.1*max(cchmap_all(i,i,:)) 
            %cell_label(1,i) = 0;
            %cchmap_all(i, i, :) = 0*cchmap_all(i, i, :);
        %end
    %end

    % ================================================-
   % cell_label_all = cell_label;
    %cell_label_new = cell_label(1, cell_label~=0);
    %clear cell_label;
   % cell_label = cell_label_new;
   % clear cell_label_new ;

    % ================================================-
    %clear location_all
    %location_all = location;
    %location_new = location(cell_label~=0,:);
    %clear location;
    %location = location_new;
    %clear location_new;

    % ================================================
    % Calculation of distances between neurons
    % ================================================
    %clear distance
    %for i = 1:size(location,1)
     %   for j = 1:size(location,1)
      %      distance(i,j) = ...
       %     sqrt((location(i,1)-location(j,1)).*(location(i,1)-location(j,1)) + ...
        %            (location(i,2)-location(j,2)).*(location(i,2)-location(j,2)));
        %end
   % end
    
    %     % ================================================
    %     % Change the bin-size.
    %     % ================================================
    %     asdf_new = ChangeBinning(asdf_raw, 1); % 2);

    % ================================================
    % Search the delay showing the maxmum
    % ================================================
    delay0 = [1:30];    %   delay1 = delay0+1;
    delay0_peaksearch = [1:16];    %   delay1 = delay0+1;
    delay0_baseline = [17:30];    %   delay1 = delay0+1;
    window_max = max(delay0);

    TEdelays = TEdelays( :, :, delay0);
    TEdelays_rand_all = TEdelays_rand_all( :, :, delay0, :);

    window_size = 3;

    %distance_max = 1200;
    %grid_num = distance_max/50;
    
    clear TEpeaks* TEdelays0_w* TEdelays0_rand_w*
%    matlabpool open
    for i = 1:size(TEdelays,1)
        i
        for j = 1:size(TEdelays,2)

            aap = delay0(TEdelays(i,j,delay0) == max(TEdelays(i,j,delay0_peaksearch)));
            TEdelays0(i,j) = aap(1);
            
            TEdelays0_window(i,j,:)          = delay0_peaksearch; % aap(1);
            TEdelays0_window_rem(i,j,:) = delay0_baseline; %  delay00;

            TEpeak0(i,j)           = TEdelays(i,j,TEdelays0(i,j));
            TEpeak0_sur(i,j)   = squeeze(sum(TEdelays(i,j,TEdelays0_window(i,j,:))));
            TEpeak0_cont(i,j) = squeeze(sum(TEdelays(i,j,TEdelays0_window_rem(i,j,:))));

            for k = 1 % :3 % size(TEdelays_rand_all,4)
                delay00 = delay0;
                aap = squeeze(delay0(TEdelays_rand_all(i,j,delay0,k) == max(TEdelays_rand_all(i,j,delay0_peaksearch,k))));
                TEdelays0_rand(i,j,k) = aap(1);

                TEdelays0_rand_window(i,j,:,k)          = delay0_peaksearch;
                TEdelays0_rand_window_rem(i,j,:,k) = delay0_baseline;

                TEpeak0_rand_all(i,j,k)           = squeeze(TEdelays_rand_all(i,j,TEdelays0_rand(i,j,k), k));
                TEpeak0_rand_all_sur(i,j,k)   = squeeze(sum(TEdelays_rand_all(i,j,TEdelays0_rand_window(i,j,:,k),k)));
                TEpeak0_rand_all_cont(i,j,k) = squeeze(sum(TEdelays_rand_all(i,j,TEdelays0_rand_window_rem(i,j,:,k),k)));

            end
        end
    end

    
%% 

    for i = 1:size(TEpeak0_rand_all, 3)
        CI_rand(:,:,i) = TEpeak0_rand_all_sur(:,:,i)./(TEpeak0_rand_all_sur(:,:,i)+TEpeak0_rand_all_cont(:,:,i));
    end

    CI = TEpeak0_sur./(TEpeak0_cont+TEpeak0_sur) ;
    %distance2 = reshape(distance, 1,size(distance,1)*size(distance,2));
%    save('CI.mat','CI','CI_rand','TEpeak0','TEpeak0_rand_all')


%    clear pdf2  pdf_sim2  pdf2_non dis_grid dis_grid TE_density TEsim_density TE_density_non
     thr = [0.01:0.01:0.80];
    for mm = 1:length(thr)
        [pdf, errorrate, imageregion] = TECIConnectivity(TEpeak0, CI, TEpeak0_rand_all, CI_rand, [25 25],  [-5 1; 0.6 1], thr(mm));
         %pdf(TEdelays00>TEpeak0) = 0;
         %str = ['mkdir  ./data/',num2str(data_name),'_graph/',num2str(thr(mm)),'_sunny;'];           eval(str);
         %str = ['save ./data/',num2str(data_name),'_graph/',num2str(thr(mm)),'_sunny/pdf pdf;'];     eval(str);
         [q1 q2] = size(pdf); cnd(mm) = (length(find(pdf~=0))/(q1*(q1-1)))*100;
         PDF(:,:,mm) = pdf; clear pdf
         ER(:,:,mm) = errorrate; clear errorrate; 
         clear pdf q1 q2
      end
     save('PDF_1_16_30ms.mat','PDF','cnd','errorrate')
%    matlabpool close
%% Plot CI-TE map

% figure(177);
% 
%   filter2 = (errorrate>0.005);
%   
%     CI_rand(isnan(CI_rand)) = 0;
%     CI_rand_thr = median(nonzeros(CI_rand));
%     logTEpeak0_rand_all = log10(nonzeros(TEpeak0_rand_all));
%     logte_thr = median(logTEpeak0_rand_all);
%     
%     x_range = [-7: 4/25: -3-4/25];
%     y_range = [0:1/25:1-1/25];
%   
%   filter2(x_range <= logte_thr, :) = 1/2;
%   filter2(:, y_range <= CI_rand_thr) = 1/2;
%               
% subplot(3,4,3);
% imagesc([-7: 4/25: -3-4/25], [1:-1/25:1/25], log10(errorrate')); 
% 
% 
% subplot(3,4,4);
% imagesc([-7: 4/25: -3-4/25], [1:-1/25:1/25], log10(errorrate')); colorbar;
% caxis([-2.3 -0.7]);
% % title([' (p = .005) '],'fontsize', 12,'fontname','Arial');
% 
% figure(177);
% subplot(3,4,9);
% imagesc([-7: 4/25: -3-4/15], [1:-1/25:1/25], filter2'); 
% title(['Filter (p = .005)'],'fontsize', 10,'fontname','Arial');
% % colorbar;
% % caxis([-2.3 -0.7]);
% % title([' (p = .005) '],'fontsize', 12,'fontname','Arial');
% 
% subplot(3,4,5)
% imagesc(pdf); colormap gray;
% title(['Connectivity map'],'fontsize', 12,'fontname','Arial');
% 
% %% DIstance dependency
%     figure(201);
%     subplot(4,4,5);
%     plot(dis_grid , TE_density+TE_density_non,'g-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     title([' Unconnected + Connected '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,9);
%     plot(dis_grid , TEsim_density./TE_density,'b-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 1],'fontsize',12,'fontname','Arial');
%     title([' Bidirectional probability '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,1);
%     plot(dis_grid , TE_density,'b-', 'linewidth', 3);  hold on;
%     plot(dis_grid , TEsim_density,'r-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     title([' (p = .005) '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,13);
%     plot(dis_grid , TE_density./(TE_density+TE_density_non),'b-', 'linewidth', 3);  hold on;
%     xlabel([' Distance [ƒÊm]'],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 0.15],'fontsize',12,'fontname','Arial');
%     title([' Connection probability '],'fontsize', 12,'fontname','Arial');
%     
% %% Network visualization
%     pdf02           = (pdf+pdf' > 0.005);
%     G = double(sparse(pdf02));
%     X = fruchterman_reingold_force_directed_layout(G);
% 
%     figure(181);
%     subplot(2,2,1);
%      plot(X(:,1),X(:,2),'ko', 'MarkerSize', 3); hold on;
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
%             X2(i,:) = X(i,:) + 2.5*(X(j,:)-X(i,:))./norm(X(j,:)-X(i,:));
%             X2(j,:) = X(j,:) + 2.5*(X(i,:)-X(j,:))./norm(X(i,:)-X(j,:));
%             if pdf02(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  ); hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%                 plot_arrow( X2(j,1)-0, X2(j,2)-0, ...
%                     X2(i,1)+0, X2(i,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%             end
%         end
%     end
% 
%     set(gca,'YLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     set(gca,'XLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     title([' (p = .005) '],'fontsize', 12,'fontname','Arial');
% 
% %% Spatial map
%     figure(202);
%     subplot(2,2,1);
%     plot(location_all(:,1),location_all(:,2),'ko', 'MarkerSize', 3); hold on;
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
% 
%             location2(i,:) = location_all(i,:) + 30*(location_all(j,:) - location_all(i,:))./norm(location_all(j,:) - location_all(i,:));
%             location2(j,:) = location_all(j,:) + 30*(location_all(i,:) - location_all(j,:))./norm(location_all(i,:) - location_all(j,:));
%             if pdf(i,j) > 0.005
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );  hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );
%                 plot_arrow( location2(j,1)-0, location2(j,2)-0, ...
%                     location2(i,1)+0, location2(i,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );
%             end
%         end
%     end
% 
%     title([' (p = .005) '],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[-1000 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[-500 500],'fontsize',12,'fontname','Arial');
%    
% %%  Degree distribution
%     figure(111);
%     subplot(2,2,1);
%     pdf_sim_num = sum(pdf_sim, 2);
%     pdf02_num = sum(pdf02,2);
%     pdf_sim_sum_num = sum(sum(pdf_sim, 2)./2);
%     pdf02_sum_num     = sum(sum(pdf02-pdf_sim,2) + sum(pdf_sim, 2)./2);
%     semilogx(pdf02_sum_num*hist(pdf02_num(:,1),10),'b-','linewidth',2);  hold on;
%     semilogx(pdf_sim_sum_num*hist(pdf_sim_num(:,1),10),'r-','linewidth',2);  hold on;
% 
%     title([' (p = .005) '],'fontsize', 12,'fontname','Arial');
%     xlabel(['Link #'],'fontsize',12,'fontname','Arial');
%     ylabel(['Frequency'],'fontsize',12,'fontname','Arial');
%      
% %%  Thr = 0.01
%     clear pdf2  pdf_sim2  pdf2_non dis_grid dis_grid TE_density TEsim_density TE_density_non
%     
%     [pdf, errorrate, imageregion] = TECIConnectivity(TEpeak0, CI, TEpeak0_rand_all, CI_rand, [25 25],  [-7 -3; 0 1], 0.01);
% 
%     % TEpeak0 ‚Æ TE at time0 ‚ð”äŠr?B ŒãŽÒ‚ª‘å‚«‚¯‚ê‚Î?Apdf = 0 ‚Æ‚·‚é?B
%      pdf(TEdelays00>TEpeak0) = 0;
%      
%     pdf2            = reshape(pdf, 1,size(pdf,1)*size(pdf,2));
%     [dis_grid TE_density] = f_gridmake( pdf2, distance2 , distance_max, grid_num);
%     
%     pdf_sim = ((pdf + pdf')==2);
%     pdf_sim2  = reshape(pdf_sim, 1,size(pdf_sim,1)*size(pdf_sim,2));
%     [dis_grid TEsim_density] = f_gridmake( pdf_sim2, distance2 , distance_max, grid_num);
% 
%     pdf2_non  = ((pdf+pdf') == 0);
%     pdf2_non2  = reshape(pdf2_non, 1,size(pdf2_non,1)*size(pdf2_non,2));
%     [dis_grid_non  TE_density_non] = f_gridmake( pdf2_non2 , distance2, distance_max, grid_num);
%     
%     addpath ./lapper/
% %     save_pdf001_19mjitt
%     str = ['mkdir  ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr001;'];           eval(str);
%     str = ['save ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr001/pdf pdf;'];     eval(str);
%     
% %% Plot CI-TE map
% 
% figure(178);
%   filter2 = (errorrate>0.01);
%   
% % subplot(3,4,1);
% subplot(2,2,1);
% semilogx(nonzeros(TEpeak0), nonzeros(CI), 'b.', 'MarkerSize', 3);
% set(gca,'XLim',[10.^(-7) 10.^(-3)] ,'fontsize',18,'fontname','Arial');
% set(gca,'YLim',[0 1] ,'fontsize',18,'fontname','Arial');
% title(['TE-CI map (raw data)'],'fontsize', 18,'fontname','Arial');
% 
% figure(179);
% subplot(2,2,2);
% semilogx(reshape(TEpeak0_rand_all,1,numel(TEpeak0_rand_all)), reshape(CI_rand,1,numel(CI_rand)), 'b.', 'MarkerSize', 3); hold on;
% set(gca,'XLim',[10.^(-7) 10.^(-3)] ,'fontsize',18,'fontname','Arial');
% set(gca,'YLim',[0 1] ,'fontsize',18,'fontname','Arial');
% title(['TE-CI map (afer jitttering)'],'fontsize', 18,'fontname','Arial');
% 
% figure(180);
% subplot(2,2,3);
% semilogx(TEpeak0(pdf == 1), CI(pdf == 1), 'bo', 'MarkerSize', 5); hold on;
% semilogx(TEpeak0(pdf == 0), CI(pdf == 0), 'r+', 'MarkerSize', 2, 'LineWidth',1);
% set(gca,'XLim',[10.^(-7) 10.^(-3)] ,'fontsize',18,'fontname','Arial');
% set(gca,'YLim',[0 1] ,'fontsize',18,'fontname','Arial');
% title(['Remain (p = .01)'],'fontsize', 18,'fontname','Arial');
%   
%   filter2(x_range <= logte_thr, :) = 1/2;
%   filter2(:, y_range <= CI_rand_thr) = 1/2;
%               
% 
% figure(177);
% subplot(3,4,10);
% imagesc([-7: 4/25: -3-4/15], [1:-1/25:1/25], filter2'); 
% title(['Filter (p = .01)'],'fontsize', 10,'fontname','Arial');
% 
% subplot(3,4,6)
% imagesc(pdf); colormap gray;
% title(['Connectivity map'],'fontsize', 10,'fontname','Arial');
%     
% %% DIstance dependency
%     figure(201);
%     subplot(4,4,6);
%     plot(dis_grid , TE_density+TE_density_non,'g-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
% 
%     subplot(4,4,10);
%     plot(dis_grid , TEsim_density./TE_density,'b-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 1],'fontsize',12,'fontname','Arial');
% 
%     subplot(4,4,2);
%     plot(dis_grid , TE_density,'b-', 'linewidth', 3);  hold on;
%     plot(dis_grid , TEsim_density,'r-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     title([' p = .01 '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,14);
%     plot(dis_grid , TE_density./(TE_density+TE_density_non),'b-', 'linewidth', 3);  hold on;
%     xlabel([' Distance [ƒÊm]'],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 0.15],'fontsize',12,'fontname','Arial');
% 
% %% Network visualization
%     pdf02           = (pdf+pdf' > 0.01);
%     G = double(sparse(pdf02));
%     X = fruchterman_reingold_force_directed_layout(G);
% 
%     figure(181);
%     subplot(2,2,2);
%      plot(X(:,1),X(:,2),'ko', 'MarkerSize', 3); hold on;
% 
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
%             X2(i,:) = X(i,:) + 2.5*(X(j,:)-X(i,:))./norm(X(j,:)-X(i,:));
%             X2(j,:) = X(j,:) + 2.5*(X(i,:)-X(j,:))./norm(X(i,:)-X(j,:));
%             if pdf02(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );  hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%                 plot_arrow( X2(j,1)-0, X2(j,2)-0, ...
%                     X2(i,1)+0, X2(i,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%             end
%         end
%     end
% 
%     set(gca,'YLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     set(gca,'XLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     title([' (p = .01) '],'fontsize', 12,'fontname','Arial');
% 
% %% Spatial map
%     figure(202);
%     subplot(2,2,2);
%     plot(location_all(:,1),location_all(:,2),'ko', 'MarkerSize', 3); hold on;
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
%             location2(i,:) = location_all(i,:) + 30*(location_all(j,:) - location_all(i,:))./norm(location_all(j,:) - location_all(i,:));
%             location2(j,:) = location_all(j,:) + 30*(location_all(i,:) - location_all(j,:))./norm(location_all(i,:) - location_all(j,:));
%             if pdf(i,j) > 0.01
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 ); hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );
%                 plot_arrow( location2(j,1)-0, location2(j,2)-0, ...
%                     location2(i,1)+0, location2(i,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );
%             end
%         end
%     end
% 
%     title([' (p = .01) '],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[-1000 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[-500 500],'fontsize',12,'fontname','Arial');
% 
% %%  Degree distribution
%     figure(111);
%     subplot(2,2,2);
%     pdf_sim_num = sum(pdf_sim, 2);
%     pdf02_num = sum(pdf02,2);
%     pdf_sim_sum_num = sum(sum(pdf_sim, 2)./2);
%     pdf02_sum_num     = sum(sum(pdf02-pdf_sim,2) + sum(pdf_sim, 2)./2);
%     semilogx(pdf02_sum_num*hist(pdf02_num(:,1),10),'b-','linewidth',2);  hold on;
%     semilogx(pdf_sim_sum_num*hist(pdf_sim_num(:,1),10),'r-','linewidth',2);  hold on;
% 
%     title([' (p = .01) '],'fontsize', 12,'fontname','Arial');
%     xlabel(['Link #'],'fontsize',12,'fontname','Arial');
%     ylabel(['Frequency'],'fontsize',12,'fontname','Arial');
% 
% %%  Thr = 0.05
%     clear pdf2  pdf_sim2  pdf2_non dis_grid dis_grid TE_density TEsim_density TE_density_non
%     [pdf, errorrate, imageregion] = TECIConnectivity(TEpeak0, CI, TEpeak0_rand_all, CI_rand, [25 25],  [-7 -3; 0 1], 0.05);
% 
%     % comparison between TEpeak0 and TE at time0.
%      pdf(TEdelays00>TEpeak0) = 0;
%      
%     pdf2            = reshape(pdf, 1,size(pdf,1)*size(pdf,2));
%     [dis_grid TE_density] = f_gridmake( pdf2, distance2 , distance_max, grid_num);
%     
%     pdf_sim = ((pdf + pdf')==2);
%     pdf_sim2  = reshape(pdf_sim, 1,size(pdf_sim,1)*size(pdf_sim,2));
%     [dis_grid TEsim_density] = f_gridmake( pdf_sim2, distance2 , distance_max, grid_num);
% 
%     pdf2_non  = ((pdf+pdf') == 0);
%     pdf2_non2  = reshape(pdf2_non, 1,size(pdf2_non,1)*size(pdf2_non,2));
%     [dis_grid_non  TE_density_non] = f_gridmake( pdf2_non2 , distance2, distance_max, grid_num);
% 
%     addpath ./lapper/
% %     save_pdf005_19mjitt
%     str = ['mkdir  ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr005;'];           eval(str);
%     str = ['save ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr005/pdf pdf;'];     eval(str);
%     
% %% Plot CI-TE map
% 
% figure(177);
%   filter2 = (errorrate>0.05);
%   
%   filter2(x_range <= logte_thr, :) = 1/2;
%   filter2(:, y_range <= CI_rand_thr) = 1/2;
%   
% subplot(3,4,11);
% imagesc([-7: 4/25: -3-4/15], [1:-1/25:1/25], filter2'); 
% title(['Filter (p = .05)'],'fontsize', 10,'fontname','Arial');
% 
% subplot(3,4,7)
% imagesc(pdf); colormap gray;
% title(['Connectivity map'],'fontsize', 10,'fontname','Arial');
% 
% %% DIstance dependency
%     figure(201);
%     subplot(4,4,7);
%     plot(dis_grid , TE_density+TE_density_non,'g-', 'linewidth', 3);  hold on;
%     %     xlabel([' Distance [ƒÊm]'],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     %     title([' p = .005 '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,11);
%     plot(dis_grid , TEsim_density./TE_density,'b-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 1],'fontsize',12,'fontname','Arial');
%     %     title([' p = .005 '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,3);
%     plot(dis_grid , TE_density,'b-', 'linewidth', 3);  hold on;
%     plot(dis_grid , TEsim_density,'r-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     %     xlabel([' Distance [ƒÊm]'],'fontsize', 12,'fontname','Arial');
%     title([' p = .05 '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,15);
%     plot(dis_grid , TE_density./(TE_density+TE_density_non),'b-', 'linewidth', 3);  hold on;
%     xlabel([' Distance [ƒÊm]'],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 0.15],'fontsize',12,'fontname','Arial');
%     %     title([' p = .005 '],'fontsize', 12,'fontname','Arial');
% 
% %% Network visualization
% 
%     pdf02           = (pdf+pdf' > 0.05);
%     G = double(sparse(pdf02));
%     X = fruchterman_reingold_force_directed_layout(G);
% 
%     figure(181);
%     subplot(2,2,3);
%      plot(X(:,1),X(:,2),'ko', 'MarkerSize', 3); hold on;
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
%             X2(i,:) = X(i,:) + 2.5*(X(j,:)-X(i,:))./norm(X(j,:)-X(i,:));
%             X2(j,:) = X(j,:) + 2.5*(X(i,:)-X(j,:))./norm(X(i,:)-X(j,:));
%             if pdf02(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  ); hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%                 plot_arrow( X2(j,1)-0, X2(j,2)-0, ...
%                     X2(i,1)+0, X2(i,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%             end
%         end
%     end
% 
%     set(gca,'YLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     set(gca,'XLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     title([' (p = .05) '],'fontsize', 12,'fontname','Arial');
% 
% %% Spatial map
%     figure(202);
%     subplot(2,2,3);
%     plot(location_all(:,1),location_all(:,2),'ko', 'MarkerSize', 3); hold on;
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
%             location2(i,:) = location_all(i,:) + 30*(location_all(j,:) - location_all(i,:))./norm(location_all(j,:) - location_all(i,:));
%             location2(j,:) = location_all(j,:) + 30*(location_all(i,:) - location_all(j,:))./norm(location_all(i,:) - location_all(j,:));
% 
%             if pdf(i,j) > 0.05
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 ); hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );
%                 plot_arrow( location2(j,1)-0, location2(j,2)-0, ...
%                     location2(i,1)+0, location2(i,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );
%             end
%         end
%     end
% 
%     title([' (p = .05) '],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[-1000 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[-500 500],'fontsize',12,'fontname','Arial');
% 
% %%  Degree distribution
%     figure(111);
%     subplot(2,2,3);
%     pdf_sim_num = sum(pdf_sim, 2);
%     pdf02_num = sum(pdf02,2);
%     pdf_sim_sum_num = sum(sum(pdf_sim, 2)./2);
%     pdf02_sum_num     = sum(sum(pdf02-pdf_sim,2) + sum(pdf_sim, 2)./2);
%     semilogx(pdf02_sum_num*hist(pdf02_num(:,1),10),'b-','linewidth',2);  hold on;
%     semilogx(pdf_sim_sum_num*hist(pdf_sim_num(:,1),10),'r-','linewidth',2);  hold on;
% 
%     title([' (p = .05) '],'fontsize', 12,'fontname','Arial');
%     xlabel(['Link #'],'fontsize',12,'fontname','Arial');
%     ylabel(['Frequency'],'fontsize',12,'fontname','Arial');
%     
% %%  Thr = 0.3
%     clear pdf2  pdf_sim2  pdf2_non dis_grid dis_grid TE_density TEsim_density TE_density_non
%     [pdf, errorrate, imageregion] = TECIConnectivity(TEpeak0, CI, TEpeak0_rand_all, CI_rand, [25 25],  [-7 -3; 0 1], 0.3);
%     
% %     save_pdf01_19mjitt
%     str = ['mkdir  ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr03;'];           eval(str);
%     str = ['save ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr03/pdf pdf;'];     eval(str);
%     
% %%  Thr = 0.2
%     clear pdf2  pdf_sim2  pdf2_non dis_grid dis_grid TE_density TEsim_density TE_density_non
%     [pdf, errorrate, imageregion] = TECIConnectivity(TEpeak0, CI, TEpeak0_rand_all, CI_rand, [25 25],  [-7 -3; 0 1], 0.2);
%     
% %     save_pdf01_19mjitt
%     str = ['mkdir  ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr02;'];           eval(str);
%     str = ['save ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr02/pdf pdf;'];     eval(str);
%     
% %%  Thr = 0.1
%     clear pdf2  pdf_sim2  pdf2_non dis_grid dis_grid TE_density TEsim_density TE_density_non
%     [pdf, errorrate, imageregion] = TECIConnectivity(TEpeak0, CI, TEpeak0_rand_all, CI_rand, [25 25],  [-7 -3; 0 1], 0.1);
% 
%     % comparison between TEpeak0 and TE at time0.
%      pdf(TEdelays00>TEpeak0) = 0;
%      
%     pdf2            = reshape(pdf, 1,size(pdf,1)*size(pdf,2));
%     [dis_grid TE_density] = f_gridmake( pdf2, distance2 , distance_max, grid_num);
%     
%     pdf_sim = ((pdf + pdf')==2);
%     pdf_sim2  = reshape(pdf_sim, 1,size(pdf_sim,1)*size(pdf_sim,2));
%     [dis_grid TEsim_density] = f_gridmake( pdf_sim2, distance2 , distance_max, grid_num);
% 
%     pdf2_non  = ((pdf+pdf') == 0);
% %     [dis_grid_non  TE_density_non] = f_gridmake( pdf2_non , distance2, distance_max, grid_num);
%     pdf2_non2  = reshape(pdf2_non, 1,size(pdf2_non,1)*size(pdf2_non,2));
%     [dis_grid_non  TE_density_non] = f_gridmake( pdf2_non2 , distance2, distance_max, grid_num);
%     
%     addpath ./lapper/
% %     save_pdf01_19mjitt
%     str = ['mkdir  ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr01;'];           eval(str);
%     str = ['save ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr01/pdf pdf;'];     eval(str);
%     
% %% Plot CI-TE map
% 
% figure(177);
%   filter2 = (errorrate>0.1);
%   
%   filter2(x_range <= logte_thr, :) = 1/2;
%   filter2(:, y_range <= CI_rand_thr) = 1/2;
%               
% 
% subplot(3,4,12);
% imagesc([-7: 4/25: -3-4/15], [1:-1/25:1/25], filter2'); 
% title(['Filter (p = .1)'],'fontsize', 10,'fontname','Arial');
% 
% subplot(3,4,8)
% imagesc(pdf); colormap gray;
% title(['Connectivity map'],'fontsize', 10,'fontname','Arial');
%     
% %% DIstance dependency
%     figure(201);
%     subplot(4,4,8);
%     plot(dis_grid , TE_density+TE_density_non,'g-', 'linewidth', 3);  hold on;
%     %     xlabel([' Distance [ƒÊm]'],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     %     title([' p = .005 '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,12);
%     plot(dis_grid , TEsim_density./TE_density,'b-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 1],'fontsize',12,'fontname','Arial');
%     %     title([' p = .005 '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,4);
%     plot(dis_grid , TE_density,'b-', 'linewidth', 3);  hold on;
%     plot(dis_grid , TEsim_density,'r-', 'linewidth', 3);  hold on;
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     title([' p = .1 '],'fontsize', 12,'fontname','Arial');
% 
%     subplot(4,4,16);
%     plot(dis_grid , TE_density./(TE_density+TE_density_non),'b-', 'linewidth', 3);  hold on;
%     xlabel([' Distance [ƒÊm]'],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[0 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[0 0.15],'fontsize',12,'fontname','Arial');
%     %     title([' p = .005 '],'fontsize', 12,'fontname','Arial');
%     
% %% Network visualization
% 
%     pdf02           = (pdf+pdf' > 0.1);
%     G = double(sparse(pdf02));
%     X = fruchterman_reingold_force_directed_layout(G);
% 
%     figure(181);
%     subplot(2,2,4);
%      plot(X(:,1),X(:,2),'ko', 'MarkerSize', 3); hold on;
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
%             X2(i,:) = X(i,:) + 2.5*(X(j,:)-X(i,:))./norm(X(j,:)-X(i,:));
%             X2(j,:) = X(j,:) + 2.5*(X(i,:)-X(j,:))./norm(X(i,:)-X(j,:));
%             if pdf02(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  ); hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 plot_arrow( X2(i,1)-0, X2(i,2)-0, ...
%                     X2(j,1)+0, X2(j,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%                 plot_arrow( X2(j,1)-0, X2(j,2)-0, ...
%                     X2(i,1)+0, X2(i,2)+0,'linewidth',0.5,'headwidth',0.16,'headheight',0.16  );
%             end
%         end
%     end
% 
%     set(gca,'YLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     set(gca,'XLim',[-min(size(pdf))./2.4 min(size(pdf))./2.4],'fontsize',12,'fontname','Arial');
%     title([' (p = .1) '],'fontsize', 12,'fontname','Arial');
% 
% %% Spatial map
%     figure(202);
%     subplot(2,2,4);
%     plot(location_all(:,1),location_all(:,2),'ko', 'MarkerSize', 3); hold on;
% 
%     for i = 1: size(pdf,1)
%         for j = 1:size(pdf,2)
%             location2(i,:) = location_all(i,:) + 30*(location_all(j,:) - location_all(i,:))./norm(location_all(j,:) - location_all(i,:));
%             location2(j,:) = location_all(j,:) + 30*(location_all(i,:) - location_all(j,:))./norm(location_all(i,:) - location_all(j,:));
%             if pdf(i,j) > 0.01
%                 %                 plot( [location_all(i,1)-0.01: location_all(j,1)- location_all(i,1)+0.02 :location_all(j,1)+0.01] ,...
%                 %                          [location_all(i,2)-0.01: location_all(j,2)- location_all(i,2)+0.02 : location_all(j,2)+0.01],'b-','linewidth',0.5,'headwidth', 3,'headheight',3 );
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 ); hold on;
%             end
%             if pdf_sim(i,j) == 1
%                 %                 plot( [location_all(i,1)-0.01: location_all(j,1)- location_all(i,1)+0.02 :location_all(j,1)+0.01] ,...
%                 %                          [location_all(i,2)-0.01: location_all(j,2)- location_all(i,2)+0.02 : location_all(j,2)+0.01],'r-','linewidth',0.5,'headwidth', 3,'headheight',3 );
%                 plot_arrow( location2(i,1)-0, location2(i,2)-0, ...
%                     location2(j,1)+0, location2(j,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3);
%                 plot_arrow( location2(j,1)-0, location2(j,2)-0, ...
%                     location2(i,1)+0, location2(i,2)+0,'linewidth',0.5,'headwidth', 3,'headheight',3 );
%             end
%         end
%     end
%     
% %     plot(location_all(cell_label_all==1 & sum(pdf02+pdf02')>= 1,1), location_all(cell_label_all==1 & sum(pdf02+pdf02')>= 1,2), 'ro', 'MarkerSize', 5,'linewidth', 2); hold on; %  'r^', 'MarkerSize', 5,'linewidth', 2); hold on;
% %     plot(location_all(cell_label_all==2 & sum(pdf02+pdf02')>= 1,1), location_all(cell_label_all==2 & sum(pdf02+pdf02')>= 1,2), 'bo', 'MarkerSize', 5,'linewidth', 2); hold on; % , 'bv', 'MarkerSize', 5,'linewidth', 2); hold on;
% 
%     title([' (p = .1) '],'fontsize', 12,'fontname','Arial');
%     set(gca,'XLim',[-1000 1000],'fontsize',12,'fontname','Arial');
%     set(gca,'YLim',[-500 500],'fontsize',12,'fontname','Arial');
% 
% %%  Degree distribution
%     figure(111);
%     subplot(2,2,4);
%     pdf_sim_num = sum(pdf_sim, 2);
%     pdf02_num = sum(pdf02,2);
%     pdf_sim_sum_num = sum(sum(pdf_sim, 2)./2);
%     pdf02_sum_num     = sum(sum(pdf02-pdf_sim,2) + sum(pdf_sim, 2)./2);
%     semilogx(pdf02_sum_num*hist(pdf02_num(:,1),10),'b-','linewidth',2);  hold on;
%     semilogx(pdf_sim_sum_num*hist(pdf_sim_num(:,1),10),'r-','linewidth',2);  hold on;
%     %     loglog(pdf02_sum_num*hist(pdf02_num(:,1),10),'b-','linewidth',2);  hold on;
%     %     loglog(pdf_sim_sum_num*hist(pdf_sim_num(:,1),10),'r-','linewidth',2);  hold on;
% 
%     title([' (p = .1) '],'fontsize', 12,'fontname','Arial');
%     xlabel(['Link #'],'fontsize',12,'fontname','Arial');
%     ylabel(['Frequency'],'fontsize',12,'fontname','Arial');
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%  Thr = 0.03
%     clear pdf2  pdf_sim2  pdf2_non dis_grid dis_grid TE_density TEsim_density TE_density_non
%     [pdf, errorrate, imageregion] = TECIConnectivity(TEpeak0, CI, TEpeak0_rand_all, CI_rand, [25 25],  [-7 -3; 0 1], 0.03);
% 
%     % TEpeak0 ‚Æ TE at time0 ‚ð”äŠr?B ŒãŽÒ‚ª‘å‚«‚¯‚ê‚Î?Apdf = 0 ‚Æ‚·‚é?B
%      pdf(TEdelays00>TEpeak0) = 0;
%      
%     pdf2            = reshape(pdf, 1,size(pdf,1)*size(pdf,2));
%     [dis_grid TE_density] = f_gridmake( pdf2, distance2 , distance_max, grid_num);
%     
%     pdf_sim = ((pdf + pdf')==2);
%     pdf_sim2  = reshape(pdf_sim, 1,size(pdf_sim,1)*size(pdf_sim,2));
%     [dis_grid TEsim_density] = f_gridmake( pdf_sim2, distance2 , distance_max, grid_num);
% 
%     pdf2_non  = ((pdf+pdf') == 0);
%     pdf2_non2  = reshape(pdf2_non, 1,size(pdf2_non,1)*size(pdf2_non,2));
%     [dis_grid_non  TE_density_non] = f_gridmake( pdf2_non2 , distance2, distance_max, grid_num);
% 
%     addpath ./lapper/
% %     save_pdf005_19mjitt
%     str = ['mkdir  ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr003;'];           eval(str);
%     str = ['save ./data_post/th8_19mszit_data30ms/',num2str(data_name),'_graph/thr003/pdf pdf;'];     eval(str);
%     
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     str = ['mkdir ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph']; eval(str);
%     str = ['save ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph/TEpeak0 TEpeak0;']; eval(str);
%     str = ['save ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph/TEpeak0_sur TEpeak0_sur;']; eval(str);
%     str = ['save ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph/TEpeak0_cont TEpeak0_cont;']; eval(str);
%     str = ['save ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph/TEpeak0_rand_all TEpeak0_rand_all;']; eval(str);
%     str = ['save ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph/TEpeak0_rand_all_sur TEpeak0_rand_all_sur;']; eval(str);
%     str = ['save ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph/TEpeak0_rand_all_cont TEpeak0_rand_all_cont;']; eval(str);
%     str = ['save ./data_post/th8TE_s19ms_d20ms3/',num2str(data_name),'_graph/distance distance;']; eval(str);
%     
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     str = ['mkdir ./figure_post/th8TE_final2/',num2str(data_name),'_graph;']; eval(str);
%     str = ['print -f179 -djpeg ./figure_post/th8TE_final2/',num2str(data_name),'_graph/TE_shuffle_19m_fig09;']; eval(str);
%     str = ['print -f179 -depsc ./figure_post/th8TE_final2/',num2str(data_name),'_graph/TE_shuffle_19m_fig09;']; eval(str);
%     str = ['print -f180 -djpeg ./figure_post/th8TE_final2/',num2str(data_name),'_graph/TE_shuffle_19m_fig010;']; eval(str);
%     str = ['print -f180 -depsc ./figure_post/th8TE_final2/',num2str(data_name),'_graph/TE_shuffle_19m_fig010;']; eval(str);
%     
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;
    





