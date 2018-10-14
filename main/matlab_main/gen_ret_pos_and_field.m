function [ basex, basey, prob, hexposx, hexposy ] = ...
    gen_ret_pos_and_field( nn, fx, fy, sdm, sc, ts )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    radi = int32(fy/nn);
    
    [basex, basey] = meshgrid(-radi:radi);
    basex=basex(:);
    basey=basey(:);
    bxt = basex;
    basex=basex(sqrt(single(basey.^2 + basex.^2))<=radi);
    basey=basey(sqrt(single(basey.^2+bxt.^2))<=radi);
    % Normal distribution centered on foveation area
    % scaled according to time-step (ts) so we get a desired p of firing
    % relative to each ms. Also scaled by just some scalar sc ~0.01?
    sd=double(radi/sdm);
    dis = sqrt(double(basex.^2 + basey.^2));
    prob=-sc*ts*2/(sqrt(3*sd)*pi^0.25).*(1-(dis./sd).^2).*exp(-0.5.*(dis./sd).^2);
    %prob=normpdf(dis, 0, sd)*sc*ts;
    %prob=prob-min(prob);
    
    
    [hexposx, hexposy] = meshgrid(1:radi:fx, 1:radi:fy);
    hexposy(:, mod(1:size(hexposy,2),2)==0) = ...
        hexposy(:, mod(1:size(hexposy,2),2)==0) + radi/2;
    hexposx=hexposx(:);
    hexposy=hexposy(:);

end

