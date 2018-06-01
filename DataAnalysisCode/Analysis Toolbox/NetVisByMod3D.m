figure; hold on;
axis vis3d;
shading interp;
colo = get(gca, 'colororder');
[X,Y,Z] = getSphere(30, 50);
for ii=1:length(groups)
    for jj=1:length(groups{ii})
        S=surf(X + positions(groups{ii}(jj),1),  ...
            Y + positions(groups{ii}(jj),2), ...
            Z + positions(groups{ii}(jj),3));
        S.FaceColor=colo(mod(ii-1,7)+1,:);
        S.FaceLighting='gouraud';
        S.LineStyle='none';
    end
end
warning('off', 'all');
wtCut = wtMat .* (wtMat > 4 | wtMat < -9 & Delays < 12);
[src, tar] = find(wtCut);
plot3([positions(src,1), positions(tar,1)]', ...
    [positions(src,2), positions(tar,2)]', ...
    [positions(src,3), positions(tar,3)]', ...
    'Color', [.5 .5 .5 .3], 'LineWidth', 1, 'LineSmoothing', 'on');
warning('on', 'all');
camlight(20, 70);
set(gca, 'Visible', 'off');
%set(gca, 'BoxStyle', 'none');
camproj('perspective');