f=figure;
set(f, 'Position', [100, 100, 700, 950]);
set(f, 'color', [0 0 0]);
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 70 95];
vid = VideoWriter('~/Desktop/Spikes.avi');
open(vid);
f.InvertHardcopy = 'off';
rastE=rast(ei,:);
rastI=rast(~ei,:);
posE=positions(ei,:);
posI=positions(~ei,:);
set(gcf, 'Renderer', 'painters');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])
box off;
grid off;

%set(gca, 'BoxStyle', 'none');
set(gca, 'LineWidth', 2);
set(gca,'Xcolor',[0.5 0.5 0.5]);
set(gca,'Ycolor',[0.5 0.5 0.5]);
set(gca,'Zcolor',[0.5 0.5 0.5]);
set(gca, 'Projection', 'Perspective');
set(gca, 'Visible', 'off');


ax = gca;
axis(ax, 'vis3d')
axes(ax)
set(ax, 'color', [0 0 0]);
xlim(ax,[0 2200]);
ylim(ax, [0 2200]);
zlim(ax, [0 6100]);


for i = 100:5000
view(ax, [(i-100)/10 -7.5]);
hold(ax, 'on');
%ax=gca;
%ax.XLimMode = 'manual';
%ax.YLimMode = 'manual';
%ax.ZLimMode = 'manual';
caxis([-1.2 1.2]);
colormap(teal2Yel);
lite = bsxfun(@times, double(wtMat~=0), double(rast(:,i)>0.3));
scatter3(ax,positions(rast(:,i)<=0.05, 1), positions(rast(:,i)<=0.05, 2),...
    positions(rast(:,i)<=0.05, 3), 40, zeros(1, sum(rast(:,i)<=0.05)), 'filled');
%%lite = lite ./max(max(lite));
%lite=lite>0.1;
[src, tar, val] = find((lite .* wtMat));
for j=1:length(src)
    if val(j)>-4 && val(j)<4
        continue;
    end
   % if val(j) < 0
   %     val(j) = -(val(j)+1.5);
   % end
    val(j) = (9*(val(j)-1)/19)+0.5;
%     if val(j) > 10
%         val(j) = 10;
%     end
if ei(src(j))
    clur = 'm';
else 
    clur = 'c';
end

alph = rast(src(j),i)/6;
alph = alph + (1-alph)*(rast(tar(j),i)*rast(src(j),i-round((Delays(src(j),tar(j))/4))))^3;
if alph > 1
    alph = 1;
end
if alph < .05
    alph = 0;
end

patchline([positions(src(j),1) positions(tar(j),1)], ...
 [positions(src(j),2) positions(tar(j),2)], ...
 [positions(src(j),3) positions(tar(j),3)],...
 'edgecolor', clur, 'LineWidth', ...
 abs(val(j)), 'edgealpha', alph, 'Parent', ax);
end
disp(i-100);

scatter3(ax,positions(rast(:,i)<=0.05, 1), positions(rast(:,i)<=0.05, 2),...
    positions(rast(:,i)<=0.05, 3), 45, 'k', 'LineWidth', 2);

scatter3(ax, posE(rastE(:,i)>0.05, 1), posE(rastE(:,i)>0.05, 2), ...
    posE(rastE(:,i)>0.05, 3), 40 + ((rastE(rastE(:,i)>0.05, i).^2).*200),...
    sqrt((rastE(rastE(:,i)>0.05, i).^2)), 'filled');
scatter3(ax, posI(rastI(:,i)>0.05, 1), posI(rastI(:,i)>0.05, 2), ...
    posI(rastI(:,i)>0.05, 3), 40 + ((rastI(rastI(:,i)>0.05, i).^2).*200),...
    -sqrt((rastI(rastI(:,i)>0.05, i).^2)), 'filled');
scatter3(ax,posE(rastE(:,i)>0.05, 1), posE(rastE(:,i)>0.05, 2),...
    posE(rastE(:,i)>0.05, 3), 40 + ((rastE(rastE(:,i)>0.05, i).^2).*200), ...
    'm', 'LineWidth', 1);
scatter3(ax,posI(rastI(:,i)>0.05, 1), posI(rastI(:,i)>0.05, 2),...
    posI(rastI(:,i)>0.05, 3), 40 + ((rastI(rastI(:,i)>0.05, i).^2).*200), ...
    'c', 'LineWidth', 1);

for k=0:15
   ps = patch([-1 2201 2201 -1], [-1 -1 2201 2201], k*6100/15 * ones(1,4), [.5 .5 .5], 'Parent', ax);
   ps.FaceAlpha=0.15;
end

%pxy = patch([-1 -1 1501 1501], [-1 1501 1501 -1], [2500 2500 2500 2500], [.5 .5 .5]);
%pxz = patch([-1 -1 1501 1501], [750 750 750 750], [-1 5001 5001 -1], [.5 .5 .5]);
%pyz = patch([750 750 750 750], [-1 -1 1501 1501], [-1 5001 5001 -1], [.5 .5 .5]);
%pxy.FaceAlpha = 0.3;
%pxz.FaceAlpha = 0.3;
%pyz.FaceAlpha = 0.3;

%set(gca,'children',flipud(get(gca,'children')))
drawnow;
writeVideo(vid,getframe(gcf));
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
