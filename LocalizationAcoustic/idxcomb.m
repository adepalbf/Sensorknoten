function A = idxcomb(phi1, phi2)

beta = [phi1, phi2];
args = {beta(1,:), beta(2,:), beta(3,:)};

NC = size(beta,1);
ii = NC:-1:1;

[A{ii}] = ndgrid(args{ii});
% concatenate
A = reshape(cat(NC+1,A{:}), [], NC);

end
