function [rNet] = createRichNet(spr, prefness, kappa, varargin)

narginchk(5, 5);
isDeg = 1;
for ii=1:2:length(varargin)
    switch varargin{ii}
        case 'RichParameter'
            isDeg = 0;
            rparam = varargin{ii+1};
            [m,n] = size(rparam);
            N=max([m,n]);
            M=N;
        case 'NumNodes'
            M = varargin{ii+1};
            N = M;
        otherwise
            error('Unknown argument');
    end
end

% TODO Create code for sparse for N > 5000 or something
rNet = zeros(M, N);

nCons = uint32(N*(N-1) * spr);

prefP = uint32(prefness * nCons);

nRCons = nCons-prefP;

conInds = randperm(N*M, nRCons);

rNet(conInds) = 1;
for ii=1:N
    rNet(ii,ii) = 0;
end

if isDeg
   rparam = double(sum(rNet,2) + sum(rNet)');
end

rParamP = rparam.^kappa ./ sum(rparam.^kappa);
figure; imagesc(rParamP);
rPPCSum = [0; cumsum(rParamP)];
figure; 
plot(rPPCSum);
figure;
hold on;
plot(rPPCSum);

counter = 0;
while prefP > 0
    
    tselect = discretize(rand, rPPCSum);
    if sum(rNet(:,tselect)) == N-1
        continue;
    end
    sselect = rParamP .* double(~(rNet(:,tselect)));
    sselect = sselect.^(2*kappa) ./ sum(sselect.^(2*kappa));
    if sum(sselect)>0
        sselect = discretize(rand, [0; cumsum(sselect)]);
%         if rNet(sselect, tselect) == 1
%             continue;
%         end
    else
        continue;
    end
    
    rNet(sselect, tselect) = 1;
    prefP = prefP - 1;
    
    if isDeg
        rparam(sselect) = rparam(sselect) + 1;
        rparam(tselect) = rparam(tselect) + 1;
        rParamP = rparam.^kappa ./ sum(rparam.^kappa);
        rPPCSum = [0; cumsum(rParamP)];
        if mod(counter, 100) == 0
            clf;
            plot(rPPCSum);
            ylim([0 1]);
            drawnow;
        end
    end
    counter = counter + 1;
    
end
% 
% switch pdFunc
%     
%     case 'lognormal'
%         
%     case 'normal'
%         
%     otherwise
        
end



