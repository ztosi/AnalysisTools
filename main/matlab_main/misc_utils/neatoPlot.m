f=figure
set(f, 'Position', [100, 100, 1400, 700]);
set(f, 'color', [0 0 0]);
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 70 70];
vid = VideoWriter('~/Desktop/Spikes_A', 'Motion JPEG AVI');
open(vid);
f.InvertHardcopy = 'off';
rastE=rastFV;
%rastI=rast(~ei,:);
positions = [x', y', zeros(length(y),1)];%isiMap11, log2(mean(cliqueDirAll+1))'];
posE= positions;%isiMaps{11}%positions(ei,:);
%posI=positions(~ei,:);
set(gcf, 'Renderer', 'painters');
%axis 'vis3d';
ax =gca;
%isa(ax,'matlab.graphics.axis.Axes')

xticks([]);
%set(ax,'xticklabel',[])
yticks([]);
%set(ax,'yticklabel',[])
zticks([]); 
%set(ax,'zticklabel',[])
box off;
%grid off;

%set(ax, 'BoxStyle', 'none');
ax.LineWidth =2;
ax.XColor = [.6 .6 .6];
ax.YColor = [.6 .6 .6];
ax.ZColor = [.6 .6 .6];
%ax.Projection = 'Perspective';
ax.Visible = 'on';
grid on;
ax.GridColor = [0.6 0.6 0.6];
%isa(ax,'matlab.graphics.axis.Axes')

ax.Color = [0 0 0];
xmin = -1000;
xmax = 1000;
ymin = -650;
ymax = 600;
zmin = 0;%-0.01;
zmax = 0;%1.01;
% xlim(ax,[0 2200]);
% ylim(ax, [0 2200]);
% zlim(ax, [0 6100]);
xlim(ax,[xmin xmax]);
ylim(ax, [ymin ymax]);
%zlim(ax, [zmin zmax]);


adj = wtMat ~= 0;
sigDotP = zeros(size(wtMat));
for i = 2000:6000 %4900
    %view(ax, [(i-100)/10 + 100, 20]);
    hold(ax, 'on');
    %ax=ax;
    %ax.XLimMode = 'manual';
    %ax.YLimMode = 'manual';
    %ax.ZLimMode = 'manual';
    caxis([-1.2 1.2]);
    %colormap(teal2Yel);
    lite = bsxfun(@times, double(wtMat~=0), double(rastFV(:,i)>0.2));
    scatter3(ax,positions(rastFV(:,i)<=0.05, 1), positions(rastFV(:,i)<=0.05, 2),...
        positions(rastFV(:,i)<=0.05, 3), 40, zeros(1, sum(rastFV(:,i)<=0.05)), 'filled');
    %%lite = lite ./max(max(lite));
    %lite=lite>0.1;
    %     [src, tar, val] = find((lite .* (wtMat*100)));
    firingNow = find(rastFV(:,i) > 0.2);
    for j=1:length(firingNow)
        %         if val(j)>-4 && val(j)<4
        %             continue;
        %         end
        % if val(j) < 0
        %     val(j) = -(val(j)+1.5);
        % end
        %     val(j) = (9*(val(j)-1)/19)+0.5;
        %     if val(j) > 10
        %         val(j) = 10;
        %     end
        %         if ei(src(j))
        clur = 'm';
        %         else
        %             clur = 'c';
        %         end
        
        outs = find(adj(:, firingNow(j)));
        inp = firingNow(j);
        for kk = 1:length(outs)
            alph = rastFV(outs(kk), i + round(Delays(inp, outs(kk))));
            alpha = alph * rastFV(inp,i);
           % alph = alph^3;
            if alph > 1
                alph = 1;
            end
            if alph > 0.05
                xyzSrc = [positions(inp,1); positions(inp,2);positions(inp,3)];
                xyzTar = [positions(outs(kk),1); positions(outs(kk),2);positions(outs(kk),3)];
                xyzPt = sigDotP(inp, outs(kk)) .* (xyzTar-xyzSrc) + xyzSrc;
                
                sigDotP(inp, outs(kk)) =  sigDotP(inp, outs(kk)) + 1/(20*(1+Delays(inp, outs(kk))));
                if  sigDotP(inp, outs(kk))> 1
                     sigDotP(inp, outs(kk)) = 0;  
                end

                patchline([xyzSrc(1), xyzTar(1)], ...
                    [xyzSrc(2), xyzTar(2)], ...
                    [xyzSrc(3), xyzTar(3)]-0.1,...
                    'edgecolor', clur, 'LineWidth', ...
                    abs(val(j))*20, 'edgealpha', alph, 'Parent', ax);
                if alph > 0.2
                    scatter3(xyzPt(1), xyzPt(2), xyzPt(3)-0.1, 30, 'm', 'filled', 'MarkerFaceAlpha', alph);
                end
            end
            
        end
        
%         alph = rast(src(j),i)/6;
%         alph = alph + (1-alph)*(rast(tar(j),i)*rast(src(j),i-round((Delays(src(j),tar(j))))))^3;
%         alph = alph * 0.8;
%         if alph > 1
%             alph = 1;
%         end
%         if alph < .1
%             alph = 0;
%         end
%         
%         patchline([positions(src(j),1) positions(tar(j),1)], ...
%             [positions(src(j),2) positions(tar(j),2)], ...
%             [positions(src(j),3) positions(tar(j),3)],...
%             'edgecolor', clur, 'LineWidth', ...
%             abs(val(j)), 'edgealpha', alph, 'Parent', ax);
    end
    disp(i);
    
    scatter3(ax,positions(rast(:,i)<=0.05, 1), positions(rast(:,i)<=0.05, 2),...
        positions(rast(:,i)<=0.05, 3), 45, 'k', 'LineWidth', 1);
    clrs = (rastE(rastE(:,i)>0.05, i));
    clrs(clrs>1) = 1;
    scatter3(ax, posE(rastE(:,i)>0.05, 1), posE(rastE(:,i)>0.05, 2), ...
        posE(rastE(:,i)>0.05, 3), 40 + ((rastE(rastE(:,i)>0.05, i).^2).*300),...
        clrs, 'filled');
    % scatter3(ax, posI(rastI(:,i)>0.05, 1), posI(rastI(:,i)>0.05, 2), ...
    %     posI(rastI(:,i)>0.05, 3), 40 + ((rastI(rastI(:,i)>0.05, i).^2).*200),...
    %     -sqrt((rastI(rastI(:,i)>0.05, i).^2)), 'filled');
    scatter3(ax,posE(rastE(:,i)>0.05, 1), posE(rastE(:,i)>0.05, 2),...
        posE(rastE(:,i)>0.05, 3), 40 + ((rastE(rastE(:,i)>0.05, i).^2).*300), ...
        'm', 'LineWidth', 1);
    % scatter3(ax,posI(rastI(:,i)>0.05, 1), posI(rastI(:,i)>0.05, 2),...
    %     posI(rastI(:,i)>0.05, 3), 40 + ((rastI(rastI(:,i)>0.05, i).^2).*200), ...
    %     'c', 'LineWidth', 1);
    
    % for k=0:15
    %    ps = patch([-1 2201 2201 -1], [-1 -1 2201 2201], k*6100/15 * ones(1,4), [.5 .5 .5], 'Parent', ax);
    %    ps.FaceAlpha=0.15;
    % end
    
    %pxy = patch([-1 -1 1501 1501], [-1 1501 1501 -1], [2500 2500 2500 2500], [.5 .5 .5]);
    %pxz = patch([-1 -1 1501 1501], [750 750 750 750], [-1 5001 5001 -1], [.5 .5 .5]);
    %pyz = patch([750 750 750 750], [-1 -1 1501 1501], [-1 5001 5001 -1], [.5 .5 .5]);
    %pxy.FaceAlpha = 0.3;
    %pxz.FaceAlpha = 0.3;
    %pyz.FaceAlpha = 0.3;
    
    %set(ax,'children',flipud(get(ax,'children')))
    drawnow;
    writeVideo(vid,getframe(gcf));
%     writeVideo(vid,getframe(gcf));
%     writeVideo(vid,getframe(gcf));

    if mod(i,100) == 0
        saveas(gcf, ['~/Desktop/Spks' num2str(i)], 'svg');
    end
    
    %     frame = getframe(f);
    %     im = frame2im(frame);
    %     [imind, cm] = rgb2ind(im, 256);
    %
    %     if i==2
    %         imwrite(imind,cm,'~/Network.gif','gif', 'DelayTime',0, ...
    %             'Loopcount', inf);
    %     else
    %         imwrite(imind,cm,'~/Network.gif','gif', 'DelayTime', 0, ...
    %             'WriteMode', 'append');
    %     end
    cla;
    hold(ax, 'off');
end


close(vid);
