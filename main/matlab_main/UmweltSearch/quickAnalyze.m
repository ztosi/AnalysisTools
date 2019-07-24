function [dat, out, randS, hidActs, map, mapHid, classNet, ce, acc] = quickAnalyze(neuNo, best, trainDatLoc, outputs, varargin)

showFig = 1;
if ~isempty(varargin)
   showFig = varargin{2};
end

if size(trainDatLoc,4) > 5000
    randS = randperm(size(trainDatLoc,4), 5000);
else
    randS = 1:size(trainDatLoc,4);
end
tb = size(trainDatLoc,1);
hid = best{2}(end-1);
inp = size(trainDatLoc,2);
samples = size(trainDatLoc,4);

dat = reshape(trainDatLoc, tb*inp, samples);
out = outputs;

map = tsne(dat(:, randS)', 'Verbose', 2, 'NumDimensions', 3, ...
    'Perplexity', 40, 'Exaggeration', 65, ...
    'Algorithm', 'barneshut', 'NumPCAComponents', int32(tb*hid/10));

hidActs = reshape(activations(best{3}, trainDatLoc, 'layer'), hid, samples);

mapHid = tsne(hidActs(:, randS)', 'Verbose', 2, 'NumDimensions', 3, ...
    'Perplexity', 30, 'Exaggeration', 65, ...
    'Algorithm', 'barneshut', 'NumPCAComponents', 10);
layerSpec =  [ ...
    imageInputLayer([size(trainDatLoc,1), size(trainDatLoc,2), 1]), ...
    fullyConnectedLayer(hid),...
    tanhLayer,...
    fullyConnectedLayer(2),...
    softmaxLayer, ...
    classificationLayer];
options = trainingOptions('adam', 'ExecutionEnvironment', ...
    'gpu', 'Verbose', 0, 'MaxEpochs', 30, 'Shuffle', 'every-epoch',...
    'InitialLearnRate', 5E-3, 'MiniBatchSize', samples/2);
classNet = trainNetwork( ...
    trainDatLoc(:,:,:, 1:end-int32(samples/10)), ...
    categorical(outputs(1:end-int32(samples/10))'==1), layerSpec, options);
y =  predict(classNet,trainDatLoc(:,:,1,(samples-int32(samples/10)+1):samples));

testSamples = int32(samples/10);
ce = zeros(testSamples,1);

for ii=1:testSamples
    ce(ii) = -log(y(ii,outputs(((samples-int32(samples/10)+ii)))+1));
end
acc = 0;
for ii=1:testSamples
    gss = round(y(ii,:));
    acc = acc + gss(outputs(((samples-int32(samples/10)+ii)))+1);
end
acc = acc / single(testSamples);

if showFig
    quickPlot(neuNo, map, mapHid, outputs, randS, ce, acc);
end
end