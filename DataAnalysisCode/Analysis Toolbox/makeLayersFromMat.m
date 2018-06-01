function [layerMem, ioScore] = makeLayersFromMat(adj_mat, varargin)

%narginchk(1, 3);
%nargoutchk(1,2);

[m, n] = size(adj_mat);
if m~=n
    % For now, just to simplify life...
    error('Matrix must be square');
end
%adj_mat = double(adj_mat~=0);
% Replace the above with the below for a weighted version
%adj_mat = abs(adj_mat);

indegs = sum(adj_mat)';
outdegs = sum(adj_mat,2);
recNess = (indegs - outdegs) ./ (indegs + outdegs);
adj_matT = adj_mat';
[~, recNI] = sort(recNess, 'ascend'); % TODO make sure to include all with 0 in-deg

if nargin == 1
    initL1 = recNI(1:int32(size(adj_mat,1)/5));
    %initL1 = randperm(size(adj_mat,1), int32(size(adj_mat,1)/5));
end

lscores = zeros(n, 2); % Layer, layer score
lscores(initL1,1) = 1;
lscores(initL1,2) = indegs(initL1)./2;
visited = zeros(n,1);
visited(initL1) = 1;
currLay = 1;
maxit = 100;
layers = -20:20;

for ii=1:200
    if ii > 1
        visited = zeros(n,1);
        visited(lscores(:,1)==currLay) = 1;
        visited(indegs==0) = 1;
    end
    %lscores(:,2) = 0;
    stLay = lscores(:,1);
    itcnt = 0;
    iter = 11;
    visitedPrev = [];
    
    while any(~visited) && ~isequal(visitedPrev, visited)
        currLay = layers(iter);
        layIni = double(lscores(:,1) == currLay);
        %size(layIni)
        LNextScore = adj_matT * layIni;
        visitedPrev = visited;
        visited = visited | (double(adj_matT~=0)*layIni);
        mvL = (lscores(:,2) < LNextScore) | (rand(size(LNextScore)) < exp(-ii/10)) ;
        nextLay = lscores(:,1) == (floor(curLay+1)) | mvL;
        penalty = sum(adj_mat(mvL, nextLay),2);
        temp = LNextScore;
        temp(mvL) = LNextScore(mvL)-penalty;
        mvL1 = mvL & (lscores(:,2) >= temp);
        mvL2 = (lscores(:,2) < (temp));
        %visited = visited & ~mvL;
        if ~isempty(mvL1)
            inBetween =  (currLay + floor(currLay+1))/2;
            if ~any(layers == inBetween)
                layers = [layers(1:iter) inBetween layers(iter+1:end)];
            end
            lscores(mvL1,2) = sum(adj_mat(lscores(:,1)==floor(inBetween), mvL1));
            lscores(mvL1,1) = inBetween;
            iter = iter+2;
        else
            iter = iter+1;
        end
        lscores(mvL2,2) = LNextScore(mvL2);
        lscores(mvL2,1) = floor(currLay)+1;
        %lscores(mvL,2) = lscores(mvL,2)-sum(adj_mat(mvL, nextLay),2);
        
        %any(~visited)
        %sum(mvL) ~= 0
        itcnt = itcnt+1;
    end
    [lys, ~, lord] = unique(lscores(:,1));
    lys = 1:length(lys);
    lscores(:,1) = lys(lord);
    layers = [(min(lys)-5) : (max(lys)+5)];
    iter = find(layers == max(lscores(:,1)));
    disp(sum(stLay ~= lscores(:,1)));
    lscores(:,1) = lscores(:,1) - min(lscores(:,1))+1;
    currLay = max(lscores(:,1));
    visited = zeros(n,1);
    visited(lscores(:,1)==currLay) = 1;
    visited(indegs==0) = 1;
    %maxLay =  max(lscores(:,1));
    %lscores(:,2)=0;
    itcnt=0;
    visitedPrev = [];
    
    while any(~visited) && ~isequal(visitedPrev, visited)
        currLay = layers(iter);
        layIni = double(lscores(:,1) == currLay);
        %size(layIni)
        LNextScore = adj_mat * layIni;
        visitedPrev = visited;
        visited = visited | (double(adj_mat~=0)*layIni);
        mvL = (lscores(:,2) < LNextScore) | (rand(size(LNextScore)) < exp(-ii/10)) ;
        %visited = visited & ~mvL;
        nextLay = lscores(:,1) == (ceil(currLay-1)) | mvL;
        penalty = sum(adj_mat(nextLay, mvL))';
        temp = LNextScore;
        temp(mvL) = LNextScore(mvL)-penalty;
        mvL1 = mvL & (lscores(:,2) >= temp);
        mvL2 = (lscores(:,2) < (temp));
        %         if ~any(mvL) && any(mvL1) %no movers after penalty
        %             lscores(lscores(:,1)<currLay) = lscores(lscores(:,1)<currLay);
        %         end
        if ~isempty(mvL1)
            inBetween = (currLay + ceil(currLay-1)) ./ 2;
            if ~any(layers == inBetween)
                layers = [layers(1:iter) inBetween layers(iter+1:end)];
            end
            lscores(mvL1,2) = sum(adj_mat(mvL1,lscores(:,1)==ceil(inBetween)),2);
            lscores(mvL1,1) = inBetween;
            iter = iter-2;
        else
            iter = iter-1;
        end
        lscores(mvL2,2) = LNextScore(mvL2);
        lscores(mvL2,1) = ceil(currLay)-1;
        
        %any(~visited)
        %sum(mvL) ~= 0
        itcnt = itcnt+1;
    end
    %lscores(:,2)=0;
    [lys, ~, lord] = unique(lscores(:,1));
    lys = 1:length(lys);
    lscores(:,1) = lys(lord);
    layers = [(min(lys)-5) : (max(lys)+5)];
    lscores(:,1) = lscores(:,1) - min(lscores(:,1))+1;
    ioScore = zeros(length(lscores),3);
    for jj=1:length(lscores)
        ioScore(jj,2) = mode(lscores(adj_mat(jj,:),1));
        ioScore(jj,1) = mode(lscores(adj_mat(:,jj),1));
    end
    ioScore(:,3) = (ioScore(:,1)+ioScore(:,2)) ./ 2;
    
end



layerMem = lscores;