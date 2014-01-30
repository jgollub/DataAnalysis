function m = maxall(A)

dims = length(size(A));

m = A;
for n=1:dims
    m = max(m);
end