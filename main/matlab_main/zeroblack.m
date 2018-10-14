cm = colormap;
cl=caxis;
if cl(1) < 0 && cl(2) > 0 
    
range = cl(2)-cl(1);

divs = range / size(cm,1);

dr = ceil((cl(1):divs:cl(2)) * (1/divs))*divs;
    
cm(dr==0,:) = [0 0 0];

colormap(cm);

clear cm cl range divs dr;

end