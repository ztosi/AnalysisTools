dnames = {'02-0' '04-0'  '05-1'  '07-0'  '08-1'...  
'03-0'  '04-1'  '06-0'  '07-1'  '09-0' ...
'03-1'  '05-0'  '06-1'  '08-0'  '09-1' };

frs = cell(15,1);
wts = cell(15,1);
xys = cell(15,1);
isiEnts = cell(15,1);
outDs = cell(15,1);

for ii=1:15
   load(string(strcat(dnames(ii), '/asdf.mat')));
   load(strcat(dnames{ii}, '/PDF_1_16_30ms', '.mat'));
   %, dnames{ii}, '.mat'));
   load(strcat(dnames{ii}, '/wgts_1_16ms.mat'));
   load(strcat(dnames{ii}, '/xy.mat'));
   
   frs{ii} = cellfun(@(x) 1000*length(x)/max(x), asdf_raw(1:end-2));
   isiEnts{ii} = cellfun(@quickDirtyEntropy, cellfun(@diff, asdf_raw(1:end-2), 'UniformOutput', 0));
   
   wts{ii} = wgt .* PDF(:,:, 48);
   outDs{ii} = sum(wts{ii}~=0,2);
   xys{ii} = [x', y'];
    
end
