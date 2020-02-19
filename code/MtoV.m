% Converting a (nx X ny) matrix to a column vector
% Counting by rows from the top
function V=MtoV(M,nx,ny)
k=1;
V=zeros(nx*ny,1);
for i=1:nx
    for j=1:ny
        V(k,1)=M(i,j);
        k=k+1;
    end
end
end