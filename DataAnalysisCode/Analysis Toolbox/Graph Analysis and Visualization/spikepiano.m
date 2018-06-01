[keyfrs, keyno] = frs2piano(PrefFRs);
kernel = exp(-(0:(1/1470):50)/10);
audbL = zeros(bperframe*framest, 88);
audbR = zeros(bperframe*framest, 88);
sinebase = (0:(1/44100):((size(audb,1)/44100)-(1/44100)))';
keyno=keyno+1;

for i=1:size(rast,1)
   % kernel = exp(-(0:(1/1470):(5*(89-keyno(i))))/(89-keyno(i)));
    train = full(rast(i,1:framest))';
    inds = find(train);
    inds = inds * bperframe;
    fulltrain = zeros(size(audb,1),1);
    szor = numel(fulltrain);
    for j=0:49
        fulltrain(inds+j)=1;
    end
    fulltrain = conv(fulltrain, kernel, 'full');
    fulltrain=fulltrain./max(fulltrain);
    fulltrain = fulltrain(1:(end-(length(fulltrain)-szor)));
    %fulltrain = fulltrain(1:szor) .* (89-keyno(i))/88;%((45-abs(44-keyno(i)))/44);
    
    %fulltrain(fulltrain>1) = 1;
    fulltrain(isnan(fulltrain))=0;
    if i<186
        audbL(:,keyno(i)) = audbL(:,keyno(i)) + fulltrain;
        audbL(audbL(:,keyno(i))>1,keyno(i)) = 1;
    else
        audbR(:,keyno(i)) = audbR(:,keyno(i)) + fulltrain;
        audbR(audbR(:,keyno(i))>1,keyno(i)) = 1;
    end
    
    disp(i);
end

kfs = (2.^(((1:88 )- 49)/12)) .* 440;

for i=1:88
    audbL(:,i) =  audbL(:,i) .* sin(sinebase * kfs(i));
    audbR(:,i) =  audbR(:,i) .* sin(sinebase * kfs(i));
end

audbL=sum(audbL,2);
audbR=sum(audbR,2);

% 
%     if i<186
%         audbL(:,keyno(i)) = audbL(:,keyno(i)) + (sin(sinebase * keyfrs(i)) .* fulltrain);
%     else
%         audbR(:,keyno(i)) = audbR(:,keyno(i)) + (sin(sinebase * keyfrs(i)) .* fulltrain);
%     end

