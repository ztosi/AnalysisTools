                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                function [diffEnt] = timingEntropy(asdf, adj_mat)

diffEnt = zeros(length(asdf)-2);

asdf = cellfun(@(x) reshape(x, [length(x), 1]), asdf, 'UniformOutput', 0);

for ii=1:size(diffEnt,1)
    inpInds = find(adj_mat(:,ii)~=0);
    ii
   for jp=1:length(inpInds)
       jj = inpInds(jp);
       indices = zeros(length(asdf{ii}) + length(asdf{jj}), 1);
       indices(1:length(asdf{ii})) = 1;
      [srtDiffs, srtInds] = sort([asdf{ii}; asdf{jj}], 'ascend');
      
      indices = indices(srtInds);
      start = 0;
      for kk=1:length(srtDiffs)
            if(indices(kk)==1)
                start = kk;
                break;
            end
      end
      
      if start == 0
          break;
      end
       tDiffs = [];
       
      chgs = diff([1; indices]);
      acu = [];
      srcFlag = 1;
      tarFlag = 0;
      tind = 1;
      for kk=1:length(chgs)
          if chgs(kk) == -1
              srcFlag = 1;
              tarFlag = 0;
              acu = [];
          end
          
          if srcFlag
              acu = [acu srtDiffs(kk)];
          end
          
          if chgs(kk) == 1
              tarFlag = 1;
              srcFlag = 0;
          end
          
          if tarFlag
              if srtDiffs(kk)-acu(end) > 50
                  continue;
              end
             tDiffs = [tDiffs (srtDiffs(kk)-acu)];
          end
          
      end
      
%       curTime = srtDiffs(start);
%       tDiffs = zeros(length(asdf{jj} - start + 1), 1);
%       tind = 1;
%       for kk=start+1:length(srtDiffs)
%             if ~indices(kk)
%                 tDiffs(tind) = curTime - srtDiffs(kk);
%                 tind = tind + 1;
%             else
%                 curTime = srtDiffs(kk);
%             end
%       end
      
      %Quick and dirty...      
       diffEnt(jj,ii) = quickDirtyEntropy(tDiffs);
       
   end
    
    
end


end