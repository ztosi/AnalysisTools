patK = [patternKey patternKey];
% patK(patK<7)=1;
% patK(patK<11 & patK>6)=2;
% patK(patK==11 | patK==12)=3;
% patK(patK==13 | patK==14)=4;
% patK(patK==15 | patK==16)=5;
% patK(patK==17 | patK==18)=6;
% patK(patK==19 | patK==20)=7;
% patK(patK==21 | patK==22)=8;
% patK(patK==23)=9;
% patK(patK==24)=10;

patK(patK<7)=1;
patK(patK<11 & patK>6)=2;
patK(patK==11 | patK==12)=3;
patK(patK==13 | patK==14)=3;
patK(patK==15 | patK==16)=3;
patK(patK==17 | patK==18)=4;
patK(patK==19 | patK==20)=4;
patK(patK==21 | patK==22)=4;
patK(patK==23)=5;
patK(patK==24)=6;