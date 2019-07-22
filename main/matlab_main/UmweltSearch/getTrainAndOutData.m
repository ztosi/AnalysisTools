function [trainingDataFull, outputs, noSamples] = getTrainAndOutData(neuNo, asdf, noSamples, kernel, ...
    mxTB, smearWindow, offset)

rast=ASDFToRaster(asdf(1:end-2), 'BinSize', 1, '-nometa', 'row');
N = length(asdf)-2;
fireTimes = find(rast(neuNo,:));
fireTimes = fireTimes(fireTimes>mxTB);
if length(fireTimes) < noSamples
    noSamples = length(fireTimes);
end
fireTimes = fireTimes(sort(randperm(length(fireTimes), noSamples), 'ascend'));
fireTimeSmear = zeros(2*smearWindow + 1, length(fireTimes));
for ii=1:length(fireTimes)
    fireTimeSmear(:, ii) = (fireTimes(ii)-smearWindow):(fireTimes(ii)+smearWindow);
end
fireTimeSmear = fireTimeSmear(:);
fireTimeSmear = unique(fireTimeSmear(fireTimeSmear>=1));
nF = zeros(size(rast, 2), 1,'logical');
nF(fireTimeSmear) = 1;
noFireTimes = find(~nF);
noFireTimes = noFireTimes(noFireTimes > mxTB);
noFireTimes = noFireTimes(sort(randperm(length(noFireTimes), noSamples), 'ascend'));

trainingDataFull = zeros(mxTB, N, 1, 2*noSamples, 'single');
figure;
for ii=1:noSamples
    kernel = (1-offset)*kernel/max(kernel);
    onInd = ii*2 - 1;
    offInd = ii*2;
    onArr = full(rast(:, (fireTimes(ii)-mxTB):(fireTimes(ii)-1)))';
    offArr = full(rast(:, (noFireTimes(ii)-mxTB):(noFireTimes(ii)-1)))';
    for jj = 1:N
        tmp = conv(onArr(:,jj), kernel, 'full');
        onArr(:, jj) = tmp(1:end-length(kernel)+1);
        tmp = conv(offArr(:,jj), kernel, 'full');
        offArr(:, jj) = tmp(1:end-length(kernel)+1);
    end
    trainingDataFull(:, :, 1, onInd) = single(onArr);
    trainingDataFull(:, :, 1, offInd) = single(offArr);
end

outputs = zeros(1,2*noSamples);
outputs(1:2:end) = 1;
trainingDataFull = trainingDataFull + offset;
noSamples = noSamples * 2;

end