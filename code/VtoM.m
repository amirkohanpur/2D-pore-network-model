% Creating a (nx X ny) matrix from a column vector
% Notice that selection of nx and ny should satisfies n=nx.ny
% Matrix will be filled by row from the top
function M=VtoM(V,n,nx,ny)
    if n~=nx*ny
        disp('Selection of nx and ny is not consistent with n.');
    end
    M=zeros(nx,ny);
    function r=res(a,b)
    r=a-b*floor(a/b);
    end
for i=1:n
    row=ceil(i/ny);
    if res(i,ny)==0
        col=ny;
    else
        col=res(i,ny);
    end
    M(row,col)=V(i);
end
end