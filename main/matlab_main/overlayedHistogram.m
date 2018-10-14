function [ h ] = overlayedHistogram( values, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    data = [];
    dcell = {};
    if iscell(values)
        numSets = length(values);
        for i = 1:length(values)
            if iscolumn(values{i})
                values{i} = values{i}';
            end
            data = [data values{i}];
            dcell(i,1) = values(i);
        end
        [~, ed] = histcounts(data, varargin{:});
        prevHeld = 1;
        if ~ishold
            figure; hold on;
            prevHeld = 0;
        end  
        h=[];
        for i = 1:length(dcell)
           h = [h histogram(dcell{i}, ed, varargin{:})];
        end
        if ~prevHeld
            hold off;
        end
    else
        numSets = 1;
        data = values;
        if iscolumn(data)
            data=data';
        end
        dcell{1, 1} = values; 
        binSpec = 0;
        edSpec = 0;
        ed = [];
        for i = 1:length(varargin)
            if isstring(varargin{i}) || isscalar(varargin{i}) ...
                    || ischar(varargin{i})
               binSpec = isscalar(varargin{i});
               numSets = i;
               edSpec = strcmp('Edges', varargin{i});
               if edSpec
                   ed = varargin{i+1};
               end
               
               break;
            end
            if iscolumn(varargin{i})
               data = [data varargin{i}'];
            else
               data = [data varargin{i}];
            end
            dcell(i+1, 1) = varargin(i);
            numSets = numSets + 1;
        end
        if ~edSpec
            [~, ed] = histcounts(data, varargin{numSets:end});
        end
        numSets = numSets + binSpec + (2*edSpec);
%         prevHeld = 1;
%         if ~ishold
%             figure; hold on;
%             prevHeld = 0;
%         end
        h= [];
        hold on;
        for i = 1:length(dcell)
           h = [h histogram(dcell{i}, ed, varargin{numSets:end})];
        end
        hold off;
%         if ~prevHeld
%            hold off;
%         end
    end
    
    alph = 1/length(h);
    if alph < 0.25
        alph = 0.25;
    end
    for i=1:length(h)
        h(i).FaceAlpha = alph;
        h(i).EdgeColor = 'none';
    end
end

