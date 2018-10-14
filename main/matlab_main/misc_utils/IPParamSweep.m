noi = 909;
%ef = EstFiringRates(noi);
%incEfs = EstFiringRates(wtMat(:,noi)~=0);
incEfs = PrefFRs(wtMat(:,noi)~=0);
incEfs=sort(incEfs,'ascend');
%incEfs=incEfs(2);
range = 0.01:0.01:70;
[pfs, efs] = meshgrid(range, range);
%pf = PrefFRs(noi);
fleeb = zeros(length(range),length(range));
%func = @(efInc, efLoc, pf) exp(-abs(efInc-efLoc)/pf)*sign(efLoc-efInc); 
raw = zeros(size(fleeb));
pfDfac = zeros(size(pfs));
pfDfac(:, range<=2) = pfs(:, range<=2)/2;
pfDfac(:, range>2) = ((1+(log(1+ (2.5*((pfs(:, range>2)/2)-1)))/2.5))); 
pfPfac = exp(-pfs/20);
for i=1:length(incEfs)
    blw=exp(-abs(incEfs(i)-efs)./pfs).*(sign(efs-incEfs(i)));
    abv = blw;
    blw = blw .* (blw>0);
    abv = abv .* (abv<0);
    abv = abv .* pfDfac;
    blw = blw .* pfPfac;
    raw = raw + abv + blw;
    i
end




% for i=1:length(incEfs)
%     i
%     fasterEfs = range>incEfs(i);
%         slowerEfs = ~fasterEfs;
%     fleeb(slowerEfs, range>2) = fleeb(slowerEfs,range>2) + ...
%         (exp(-pfs(fasterEfs,range>2)/20) .* ...
%         exp(-abs(efs(fasterEfs,range>2)-incEfs(i))./pfs(fasterEfs, range>2)));
%     fleeb(fasterEfs, range>2) = fleeb(fasterEfs,range>2) - ...
%         ((1+(log(1+ (4*((pfs(fasterEfs, range>2)/2)-1)))/4)) .* ...
%         exp(-abs(efs(slowerEfs, range>2)-incEfs(i))./pfs(slowerEfs, range>2)));
% end
% for i=1:length(incEfs)
%     i
%     fasterEfs = range>incEfs(i);
%         slowerEfs = ~fasterEfs;
%     fleeb(slowerEfs, range<=2) = fleeb(slowerEfs, range<=2) + ...
%         (exp(-pfs(fasterEfs, range<=2)/20) .* ...
%         exp(-abs(efs(slowerEfs, range<=2)-incEfs(i))./pfs(slowerEfs,range<=2)));
% 
%     fleeb(fasterEfs, range<=2) = fleeb(fasterEfs,range<=2) - ...
%         (pfs(fasterEfs, range<=2)./2 .* ...
%         exp(-abs(efs(fasterEfs, range<=2)-incEfs(i))./pfs(fasterEfs, range<=2)));
% end
% 

