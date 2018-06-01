names = {'02-0',  '03-1',  '04-1', '05-1',  '06-1', '07-1',  '08-1', ... % '09-1', ...
'03-0' , '04-0',  '05-0',  '06-0',  '07-0',  '08-0'  '09-0' };
names = sort(names);
cd '/home/zoe/Documents/Sunny Data/';
asdfPrefix = 'asdf_data';
wvfrmPrefix = 'WaveformCellLabel';
% sumSqDists = cell(length(names),1);
% wvjsds = cell(length(names),1);
%isijsds = cell(length(names),1);
frs = cell(length(names),1);
excInds = cell(length(names),1);
tsneWVJsd = cell(length(names),1);
%tsneRawWV = cell(length(names),1);
tsneComb = cell(length(names),1);
allwvf = [];
szs = zeros(1,14);
counter = 0;
for kk=1:length(names)
   disp(['Dataset ' num2str(kk) '__________________________________' ]);
   load(['./' asdfPrefix '/' num2str(kk) '.mat']);
   load(['../Spike Waveforms/' wvfrmPrefix names{kk} '.mat']);
   excInds{kk} = excind; 
   frs{kk} = cellfun(@(x) 1000*length(x)/asdf_raw{end}(2), asdf_raw(1:end-2));
   %isijsd = ssJSD_DistMatrix(asdf_raw, [.05 .2 .5 .2 .05]);
 %  isijsds{kk} = isijsd ./ max(max(isijsd));
   N = length(asdf_raw)-2;
   szs(ii)=N;
   wvjsd = zeros(N);
   ssd = zeros(N);
   nwfN = nwf - min(min(nwf));
   if length(nwf(1,:))==41
    allwvf = [allwvf; nwfN];
    counter=counter+1;
   end
%    for ii=1:N
%        for jj=ii+1:N
%            [wvjsd(ii,jj), ssd(ii,jj)] = curveCopare(nwfN(ii,:), nwfN(jj,:));
%        end
%    end
%    wvjsds{kk} = (wvjsd+wvjsd') ./ max(max(wvjsd));
%    sumSqDists{kk} = (ssd + ssd') ./ max(max(ssd));
   %tsneRawWV{kk} = tsne(nwfN, [], 3, 15, 25);
   tsneWVJsd{kk} = tsne_d(wvjsds{kk}, [], tsneISI{kk}, 35);
   %tsneComb{kk} = tsne_d(wvjsds{kk}.*0.5 + isijsds{kk}.*0.5, [], 3, 15);
end