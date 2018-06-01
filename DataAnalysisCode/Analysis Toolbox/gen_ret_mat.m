function [ bx, by, p, hx, hy, out, im] = gen_ret_mat(  nn, fx, fy, sdm, sc, ts )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
[bx, by, p, hx, hy] = gen_ret_pos_and_field(nn, fx, fy, sdm, sc, ts);

inds=cell(1, length(hx));
vals=cell(1, length(hx));
im=zeros(144,180);
for i=1:length(hx)
    [inds{i}, okv] = pos_subs([fy fx], hy(i)+by, hx(i)+bx);
    inds{i}=inds{i}+((i-1)*fy*fx);
    vals{i}=p(okv)';   
    im(inds{i}-((i-1)*fy*fx))=im(inds{i}-((i-1)*fy*fx))+vals{i};
end
%figure; imagesc(im);
[I,J] = ind2sub([fx*fy length(hx)], cell2mat(inds));
out=sparse(J,I,cell2mat(vals), length(hx), fx*fy);
end

