function [conmat, errorrate, imageregion] = TECIConnectivity(mte, ci, mte_j, ci_j, resol, edge, thr)
% function [conmat, errorrate] = TECIConnectivity(mte, ci, jmte, jci, resol, edge)
%
%    mte - (c,1)
%    ci - {d,2}
%    mte_j - 
%    ci_j - 
%    resol - (1,2) resolution of each axis (default: [25 25])
%    edge - (2,2) edge of each axis (default: [-7 -3; 0 1])
%    thr - (1,1) threshold for error rate (default: 0.05)
%
% Returns:
%
%
% Description :
%
%
% Example :
%    
%
% Author   : Shinya Ito
%            Indiana University
%
% Last modified on 6/25/2011

% remove diagonal elements
mte = RemoveDiag(mte);
mte_j = RemoveDiag(mte_j);
ci = RemoveDiag(ci);
ci_i = RemoveDiag(ci_j);

if nargin < 5
	resol = [25 25]; % 25 by 25 by default
end

if nargin < 6
	edge = [-7 -3; 0 1]; % some reasonable thing
end

if nargin < 7
	thr = 0.05;
end

xincr = (edge(1,2) - edge(1,1))/resol(1);
xregions = 10.^(edge(1,1):xincr:edge(1,2));

yincr = (edge(2,2) - edge(2,1))/resol(2);
yregions = edge(2,1):yincr:edge(2,2);

js = size(mte_j);

if length(js) < 3
	jitnum = 1;
else
	jitnum = js(3); % 3rd dimension is the number of jittered data.
end

imageregion = zeros(resol(1), resol(2));
errorrate = zeros(resol(1), resol(2));
connected = [];


for i = 1:resol(1)
	isec = find(mte > xregions(i) & mte <= xregions(i+1));
	isecj = find(mte_j > xregions(i) & mte_j <= xregions(i+1));

	for j = 1:resol(2)
		jsec = find(ci > yregions(j) & ci <= yregions(j+1));
		jsecj = find(ci_j > yregions(j) & ci_j <= yregions(j+1));

		n = length(intersect(isec, jsec));
		if n ~= 0
			nj = length(intersect(isecj, jsecj)) / jitnum;
			errorrate(i,j) = nj/(n+nj);
			if nj/(n+nj) < thr
				%if i>10 % lower bound, needed?
					connected = [connected; intersect(isec, jsec)];
				%end
				imageregion(i,j) = 1;
			end
		end
	end
end

nNeu = length(mte);
conmat = zeros(nNeu);
conmat(connected) = 1;
