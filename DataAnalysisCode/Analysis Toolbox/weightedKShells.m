function [k_parts, nodeId] = weightedKShells(wtMat, iob)

wtMat = abs(wtMat);
nodeId = 1:size(wtMat,2);
k_parts = zeros(1, size(wtMat,2));
nzrs = (sum(wtMat) .* sum(wtMat,2)') ~= 0;
wtMat = wtMat(nzrs, nzrs);
nodeId = nodeId(nzrs);

wtMat = wtMat ./ mean(nonzeros(wtMat(:)));

wtMat = wtMat ./ min(nonzeros(wtMat(:)));

wtMat = ceil(wtMat);

if strcmp(iob, 'both')
    
    kd = ceil(sqrt((sum(wtMat~=0) + sum(wtMat~=0,2)') ...
        .* (sum(wtMat) + sum(wtMat,2)')));
    level = min(kd);
    
    while any(k_parts == 0) && ~isempty(nodeId)
        
        while any(kd<=level)
            k_parts(nodeId(kd<=level)) = level;
            nodeId = nodeId(kd>level);
            wtMat = wtMat(kd>level, kd>level);
            kd = ceil(sqrt((sum(wtMat~=0) + sum(wtMat~=0,2)') ...
                .* (sum(wtMat) + sum(wtMat,2)')));
        end
        level = level + 1;
        
    end
    
elseif strcmp(iob, 'in')
    kd = ceil(sqrt((sum(wtMat~=0)) ...
        .* (sum(wtMat))));
    level = min(kd);
    
    while any(k_parts == 0) && ~isempty(nodeId)
        
        while any(kd<=level)&& ~isempty(nodeId)
            k_parts(nodeId(kd<=level)) = level;
            nodeId = nodeId(kd>level);
            wtMat = wtMat(kd>level, kd>level);
            kd = ceil(sqrt((sum(wtMat~=0) ) ...
                .* (sum(wtMat))));
        end
        level = level + 1;
        
    end
    
elseif strcmp(iob, 'out')
    kd = ceil(sqrt(( sum(wtMat~=0,2)') ...
        .* (sum(wtMat,2)')));
    level = min(kd);
    while any(k_parts == 0) && ~isempty(nodeId)
        
        while any(kd<=level)&& ~isempty(nodeId)
            k_parts(nodeId(kd<=level)) = level;
            nodeId = nodeId(kd>level);
            wtMat = wtMat(kd>level, kd>level);
            kd = ceil(sqrt((sum(wtMat~=0,2)') ...
                .* ( sum(wtMat,2)')));
        end
        level = level + 1;
        
    end
    
    
    
end