function [nn, selection, perCorrect, trainingData, trOut, rast] = GASearch(neuNo, N, asdf, varargin)
%%
% Searches for the optimal set of neurons in the network to use to predict
% whether or not a given neuron "neuNo" will fire.
%
%
%%
% 1-N: 1 where a given input neuron is selected
% N+1 : # of hidden units
% N+2 : # of time frames to look back (5-20)
% TODO: num agents is argument
numAgents = 200;
numSurvivors = 20;
numGenerations = 250;
noSamples = 20000;
mutationP = 0.01;
mxTB = 25;%25;
mnTB = 5;
mxHid = 100;%100;
mnHid = 3;
agentGenes = ones(N+2, numAgents);
agentGenes(1:N,:) = rand(N, numAgents) < repmat(rand(1, numAgents)*0.5 + 0.01 , N, 1);
agentGenes(neuNo,:) = 0; % not allowed to use your own activation
agentGenes(N+1,:) = randi([mnHid mxHid], 1, numAgents); % 10-100 hidden units
agentGenes(N+2,:) = randi([mnTB mxTB], 1, numAgents); % looking back between 5 and 20 time steps
agents = cell(numAgents,1);
fitness = ones(numAgents,1);
accuracy = zeros(numAgents,1);
best = cell(3,1);
algorithm = 'MLP';
useGPU = false;
verbose = 2;
showFig = 1;
tinyBias = 0.1;
dsN = '';

for ii=1:2:length(varargin)
    switch varargin{ii}
        case 'NumAgents'
            numAgents = varargin{ii+1};
        case 'Generations'
            numGenerations = varargin{ii+1};
        case 'MutationProb'
            mutationP = varargin{ii+1};
        case 'NumSamples'
            noSamples = varargin{ii+1};
        case 'Survivors'
            numSurvivors = varargin{ii+1};
        case 'Algorithm'
            algorithm = varargin{ii+1};
            if ~(strcmp(algorithm, 'MLP') || strcmp(algorithm, 'LSCOV'))
                error('Uknown Algorithm option. Options are MLP and LSCOV');
            end
        case 'UseGPU'
            useGPU = varargin{ii+1};
        case 'Verbose'
            verbose = varargin{ii+1};
            if verbose > 2 || verbose < 0
                error('Verbose levels are 0 (off), 1 or 2');
            end
        case 'ShowFigure'
            showFig = varargin{ii+1};
        case 'SmallNetBias'
            tinyBias = varargin{ii+1};
        case 'DataSetName'
            dsN = varargin{ii+1};
        otherwise
            error('Unknown input.');
    end
end

if useGPU
    if strcmp(algorithm, 'LSCOV')
        warning('UseGPU is only a valid option with MLP. Ignoring...');
    end
end

[trainingDataFull, outputs, noSamples] = getTrainAndOutData(neuNo, asdf, noSamples, exp(-[0:18]/3)', mxTB, 10, -0.1);
disp(['Number of Samples: ', num2str(noSamples)]);
if showFig
    figure('position',[0, 0, 1200, 800]);
end

outCat = zeros(noSamples,2,'single');
for ll=1:noSamples
   ind = outputs(ll)+1;
   outCat(ll,ind) = 1;
end

for ii=1:numGenerations
    if verbose >= 1
        disp(['GENERATION ' num2str(ii)]);
    end
    for jj=1:numAgents
        tic;
        ce = 0; % cross entropy
        acc = 0; % accuracy
        if(sum(agentGenes(1:N,jj)) == 0) % accidentally evolved an agent w/ no genes
            agentGenes(randperm(N, 10),jj) = 1;
        end
        
        t_back = agentGenes(N+2, jj);
        noHid = agentGenes(N+1, jj);
        inps = agentGenes(1:N, jj)==1;
        
        
        
        if strcmp(algorithm, 'MLP') %Use a multi-layer perceptron classifier
            outputLoc = single(outputs');
            if useGPU
                outputLoc = categorical(outputLoc);
                trainDatLoc = input4SelectAndType(inps, t_back, trainingDataFull, 'image');
                %plorp = [size(trainDatLoc,1), size(trainDatLoc,2), 1];
                %if ((size(plorp,2)<2) || (size(plorp,2)>3) || (nnz(plorp)<3))
                %    disp('WAT');
                %end
                layerSpec =  [ ...
                    imageInputLayer([size(trainDatLoc,1), size(trainDatLoc,2), 1]), ...
                    fullyConnectedLayer(noHid),...
                    tanhLayer,...
                    fullyConnectedLayer(2),...
                    softmaxLayer, ...
                    classificationLayer];
                                   % regressionLayer];
                options = trainingOptions('adam', 'ExecutionEnvironment', ...
                    'gpu', 'Verbose', 0, 'MaxEpochs', 30, 'Shuffle', 'every-epoch',...
                    'InitialLearnRate', 5E-3, 'MiniBatchSize', noSamples/2);
                
                net = trainNetwork( ...
                    trainDatLoc(:,:,:, 1:end-int32(noSamples/10)), ...
                    outputLoc(1:end-int32(noSamples/10)), layerSpec, options);
                
            else
                trainDatLoc = input4SelectAndType(inps, t_back, trainingDataFull, 'vector');
                net = train(net, trainDatLoc(:, 1:end-int32(noSamples/10)), ...
                    outCat(1:end-int32(noSamples/10)), 'UseParallel', 'yes', 'Reduction', 1);
                
            end
            
            % Get a test MSE
            y =  predict(net,trainDatLoc(:,:,1,(noSamples-int32(noSamples/10)+1):noSamples));
            testOut = outputs((end-int32(noSamples/10)+1):end)+1;
            for kk=1:length(testOut)
               ce = ce - log(y(kk,testOut(kk)));
               [~,sl] = max(y(kk,:));
               acc = acc + (sl == testOut(kk));
            end
            acc = acc / length(testOut);
            %mseLoc = mean(((single(y)')-outputs((end-int32(noSamples/10)+1):end)).^2);
            accuracy(jj) = acc;%mseLoc;
            fitness(jj) = ce * (1+(tinyBias*(sum(inps)/N)));%mses(jj) + (sum(inps)/N) * tinyBias;
            agents{jj} = net;
            
        else
            outputLoc = outCat;
            trainDatLoc = input4SelectAndType(inps, t_back, trainingDataFull, 'vector');
            net = lscov(trainDatLoc(:, 1:end-int32(noSamples/10))', outputLoc(1:end-int32(noSamples/10), :));
            y = trainDatLoc(:, (end-int32(noSamples/10)+1):end)'*net;
            y = exp(y)./sum(exp(y),2); % softmax
            testOut = outputs((end-int32(noSamples/10)+1):end)+1;
            for kk=1:length(testOut)
               ce = ce - log(y(kk,testOut(kk)));
               [~,sl] = max(y(kk,:));
               acc = acc + (sl == testOut(kk));
            end
            acc = acc / length(testOut);
            %mseLoc = mean((y-outputs((end-int32(noSamples/10)+1):end)').^2);
            accuracy(jj) = acc;%mseLoc;
            fitness(jj) = ce * (1+(tinyBias*(sum(inps)/N)));%mseLoc + (sum(inps)/N) * tinyBias;
            agents{jj} = net;
        end
        if isempty(best{1}) || best{1} > acc
            best{1} = acc;
            best{2} = agentGenes(:, jj);
            best{3} = net;
            save(['Neu',num2str(neuNo),'_', dsN,'umweltSelect.mat'], 'best', 'trainDatLoc', 'outputs', 'ce', 'acc');
        end
        if verbose == 2
            trainDatLoc(trainDatLoc<0) = 0;
            if ndims(trainDatLoc) == 4
                trainDatLoc = reshape(trainDatLoc, size(trainDatLoc,1) * ...
                    size(trainDatLoc,2), size(trainDatLoc,4));
            end    
            st = sum(trainDatLoc);
            perZero = sum(st==0)/length(st);
            st = find(st == 0);
            spo = sum(outputs(st))/length(st);
            disp(['Agent: ', num2str(jj), '  Fitness: ', num2str(fitness(jj)),...
                '  Accuracy: ', num2str(acc*100),'%', ...
                '  Zero%: ', num2str(100*perZero), ...
                '  Spont%: ', num2str(100*spo)]);
        end
        toc;
    end
    
    % Begin the "genetic" part of the GA
    [~,I] = sort(fitness, 'ascend');
    slcts = agentGenes(1:N,I(1:numSurvivors)); % Genes of the elites
    toDie = I(numSurvivors+1:end); % the indices of all to die
    
    % recombined (mate), mutate, and generate random agents
    genesTemp = recombAndMut(agentGenes, I(1:numSurvivors), I(1:2*numSurvivors), length(toDie), ...
        [mnTB mxTB], [mnHid mxHid], .7, .2, mutationP);
    agentGenes(:, toDie) = genesTemp;
    %agentGenes(N+1,agentGenes(N+1,:)<mnHid) = mnHid;
    
    if showFig
        clf;
        subplot(2, 3, [1 4]);
        hold on;
        title(['Gen: ', num2str(ii), ' Elite Agent Inputs']);
        imagesc(slcts);
        xlim([0.5 numSurvivors+.5]);
        ylim([0.5 N+.5]);
        hold off;
        subplot(2,3,3);
        hold on;
        title('Accuracy');
        histogram(accuracy(I(1:numSurvivors)));
        hold off;
        subplot(2,3, [2 5]);
        hold on;
        title('Frequency of input appearance');
        plot( [zeros(1,N); sum(slcts,2)'], [1:N;1:N], '-k');
        scatter(sum(slcts, 2), 1:N, 30, 'filled');
        ylim([.5, N+.5]);
        hold off;
        subplot(2,3,6);
        if strcmp(algorithm, 'MLP')
            hold on;
            title('Input vs Hidden Units');
            scatter(sum(slcts), agentGenes(N+1, I(1:numSurvivors)), 30, 'filled');
            ylim([mnHid-1 mxHid+1]);
            unityline(gca);
            hold off;
        else
            hold on;
            title('# of input neurons');
            histogram(sum(slcts));
            hold off;
        end
        drawnow;
    end
    
    if verbose >= 1
        disp(['Best Accuracy: ', num2str(max(accuracy))]);
        disp(['Best Fitness: ', num2str(min(fitness))]);
        [~,miFit] = min(fitness);
        disp(['Hidden Units: ', num2str(agentGenes(N+1, miFit))]);
        disp(['Looks Back: ', num2str(agentGenes(N+2, miFit))]);
        disp(['Selected Neurons, ', num2str(sum(agentGenes(1:N, miFit))), ':']);
        disp(find(agentGenes(1:N,miFit)));
        disp(' ');
    end
    
end

[~,miFit] = min(accuracy);
perCorrect = best{1};%mses(miFit);
selection = find(best{2}(1:N)); %find(agentGenes(1:N, miFit));
nn = best{3};%agents{miFit};
t_back = agentGenes(N+2,miFit);
trainingData = trainingDataFull(end-t_back+1:end, selection, :);
trOut = outputs;

end