function [ out ] = convM( M,v,shape )
%[ out ] = convM( M,v ) convolves every line of M by vector v

if nargin < 3,
    shape = 'full';
end

out = [];
for i = 1:size(M,1)
     out(i,:) = conv(M(i,:),v,shape);
end

end

