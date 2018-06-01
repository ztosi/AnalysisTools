tic;
wtMatO = wtMat; %Save original
wtMat = wtMat .* (abs(wtMat)>0.01);
[N, m] = size(wtMat);
nullMods = zeros(N,m,100, 'uint8');
nullMods10 = zeros(N,m,100, 'uint8');
nullModsee = zeros(sum(ei),sum(ei),100, 'uint8');
nullModsee10 = zeros(sum(ei),sum(ei),100, 'uint8');
wtValsEx = nonzeros(wtMat(ei,ei));
wtValsEx = sort(wtValsEx, 'descend');
cutoff10 = wtValsEx(uint32(length(wtValsEx)/10));
wtVals=nonzeros(abs(wtMat));
wtVals=sort(wtVals, 'descend');
locCut=wtVals(uint32(length(wtVals)/10));
inpKins = sum(inpmat~=0);
inpKins10 = sum(inpmat>locCut);
inpKinsee10 = sum(inpmat>cutoff10);
wtMEx10 = wtMat(ei,ei)>cutoff10;
wtMatee = wtMat(ei,ei);
wtMat10 = abs(wtMat)>locCut;
wtMatee10 = wtMEx10;
disp('Generating null models...');
tic;
parfor i=1:100
    nullMods(:,:,i) = uint8(dir_generate_srand_bid_prev(wtMat~=0));
    nullModsee(:,:,i) = uint8(dir_generate_srand_bid_prev(wtMat(ei,ei)~=0));
    nullModsee10(:,:,i) = uint8(dir_generate_srand_bid_prev(wtMEx10));
    nullMods10(:,:,i) = uint8(dir_generate_srand_bid_prev(wtMat10));
end
wtMat = wtMat ~= 0;
wtMatee = wtMatee ~= 0;
toc;
disp('DONE');
