function [ optima, flags ] = findOptima( vec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    d1 = vec(2:end)-vec(1:end-1);
    negDer = d1 < 0;
    ndd = diff(negDer);
    optima=[];
    flags=[];
    n=0;
    if sum(ndd<0)~=0
        optima = find(ndd<0);
        n=length(optima);
        flags = zeros(1,n); % zeros represent attractors
            
    end
    if sum(ndd>0)~=0
        if ~isempty(optima)
        optima = [optima, find(ndd>0)];
        flags = [flags, ones(1, length(optima)-n)]; % ones represent repellors
        else
            optima = find(ndd>0);
            flags = ones(1,length(optima)-n); 
        end
        
        n=length(optima);
    end
    
    if ~isempty(optima)
        optima = [1 optima];
        flags = [negDer(1)>0 flags];
    else
        optima = 1;
        flags = negDer(1) > 0;
    end
    
%     if sum(ndd==0)~=0
%         if ~isempty(optima)
%             optima = [optima, find(ndd==0)];
%             flags = [flags, 0.5.* ones(1,length(optima)-n)];
%         else
%             optima = find(ndd==0);
%             flags = 0.5.* ones(1,length(optima)-n);
%         end
%     end
end

