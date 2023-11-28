function y = temp20230511(x,theta,n)

% a localGenerator that overwrites symbol positions uniformly sampled (with
% replacement) in a suffix determined by the realized position of an
% initial symbol overwrite. Cf. temp20230418.m

y = x;
ind = randi([1,numel(y)]);
y(ind) = randi(n);
for j = 2:ceil(theta)
    ind_j = randi([ind,numel(y)]);
    y(ind_j) = randi(n);
end