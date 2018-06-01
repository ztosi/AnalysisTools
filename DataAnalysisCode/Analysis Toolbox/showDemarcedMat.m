function showDemarcedMat(mMem, wtMat)

noMem = length(mMem);
modOrd = cell2mat(cellfun(@(x) x',mMem, 'UniformOutput', 0));
mMem = cellfun(@(x) ones(length(x),1), mMem, 'UniformOutput', 0);
for ii=1:noMem
    mMem{ii} = mMem{ii} .* ii;
end
mMem=cell2mat(mMem);

[mMemsUnq,~, ord] = unique(mMem);
mMemsUnq = 1:length(mMemsUnq);
mMem = mMemsUnq(ord);
offset = 0.5;


figure;
hold on;
if (any(sum(wtMat,2)<0))
    teal2Yel = zeros(65536,3);
    teal2Yel(1:end/2, 1) = (1-(1/32768)):-1/32768:0;
    teal2Yel(end/2:end, 2) = 0:1/32768:(1);
    teal2Yel(end/2:end, 3) = 0:1/32768:(1);
    teal2Yel(1:end/2, 2) = (1-(1/32768)):-1/32768:0;
    colormap(flipud(teal2Yel));
else
%     colormap(flipud(gray(256)));
end
[m, n] = size(wtMat(modOrd,modOrd));
imagesc(wtMat(modOrd, modOrd));

for ii = 1:length(mMemsUnq)-1
    mems = sum(mMem==mMemsUnq(ii));
    plot([0.5 n+0.5], [mems mems]+offset, 'm-', 'LineWidth',2);
    plot([mems mems]+offset, [0.5 m+0.5], 'm-', 'LineWidth',2);
    offset = offset+mems;
end
xlim([1 m]);
ylim([1,n]);
view([90 90]);
drawnow;

% [J, I] = find(wtMat');
% layOff = 0;
% figure; hold on;
% colorOrd = get(gca, 'ColorOrder');
% xycoors = zeros(length(I),2);
% for ii = 1:length(mMemsUnq)
%     mems = sum(mMem==mMemsUnq(ii));
%     xycoors(mMem==mMemsUnq(ii),1) = (1:mems)-mems/2;
%     xycoors(mMem==mMemsUnq(ii),2) = layOff;
%     layOff = layOff + 10;
% end
% 
% for ii = length(mMemsUnq):-1:1
%    color = colorOrd(mod(mMemsUnq(ii)-1,7)+1,:);
%    src = find(mMem == mMemsUnq(ii));
%    srcs = ismember(I, src);
%    tars = J(srcs);
%    srcs = I(srcs);
%    plot([xycoors(srcs,1), xycoors(tars,1)]', ...
%        [xycoors(srcs,2), xycoors(tars,2)]', '-', 'Color', [color 0.2]);
% end
% 
% for ii = 1:length(mMemsUnq)
%    mems = mMem==mMemsUnq(ii);
%    color = colorOrd(mod(mMemsUnq(ii)-1,7)+1,:);
%    scatter(xycoors(mems,1), xycoors(mems,2), 20, color);
% end

hold off;
