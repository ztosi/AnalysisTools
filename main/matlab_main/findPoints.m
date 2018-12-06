function [inds] = findPoints(query, pts, tol)

n = length(query);
[M,N] = size(pts);

if M ~= n
   if N ~= n
       error('Query point must have the same dimension as points in pts.');
   else
       pts = pts';
       tmp = M;
       M = N;
       N = tmp;
   end
end

if isrow(query)
    query = query';
end

inds = [];
for ii=1:N
   if ~any(abs(pts(:,ii)-query) > tol)
      inds = [inds; ii]; 
   end
end

end