% START
%

clear;

% physical properties
rho=1; mu=1; Pmax=11000; Pmin=0;

% number of pore bodies in the network
nx=20; ny=20;
nB=nx*ny; nT=nx*(ny+1)+(nx-1)*ny;

% required data for random distribution
Bmin=0.12e-3; Bmax=0.24e-3; Bmean=0.188e-3; Bsd=0.25;
Tmin=0.010e-3; Tmax=0.102e-3; Tmean=0.024e-3; Tsd=0.25;

% a sample size for network
lx=0.3e-3; ly=0.3e-3;
Lx=(nx-1)*lx; Ly=(ny-1)*ly;

% indtroducing the matrices
Bxp=zeros(ny,nx); Byp=zeros(ny,nx); 
TVd=zeros(ny+1,nx); TVl=zeros(ny+1,nx); TVxp=zeros(ny+1,nx); TVyp=zeros(ny+1,nx);
THd=zeros(ny,nx-1); THl=zeros(ny,nx-1); THxp=zeros(ny,nx-1); THyp=zeros(ny,nx-1);

% variables for the linear system of the problem
% AX=b with initial guess X0 ; A=BA, X=BP, b=Bb,
BA=zeros(nB);
Bb=zeros(nB,1);

% random truncated log-normal distrubtion of pore bodies and pore throats
% random data are diameters
Bd=lognm(ny,nx,Bmin,Bmax,Bmean,Bsd);
Td=lognv(nT,Tmin,Tmax,Tmean,Tsd);

% splitting throats vector into two parts (vertical & horizontal)
% vertical throats
for i=1:(nx*(ny+1))
    row=ceil(i/nx);
    if resdiv(i,nx)==0
        col=nx;
    else
        col=resdiv(i,nx);
    end
    TVd(row,col)=Td(i);
end
% horizontal throats
for j=((nx*(ny+1))+1):nT
    i=j-nx*(ny+1);
    row=ceil(i/(nx-1));
    if resdiv(i,(nx-1))==0
        col=nx-1;
    else
        col=resdiv(i,(nx-1));
    end
    THd(row,col)=Td(j);
end

% computing the lengths of pore throats
% vertical throats
for i=1:ny+1
    for j=1:nx
        if i==1
            TVl(i,j)=ly-0.5*Bmean-0.5*Bd(i,j);
        elseif i==(ny+1)
            TVl(i,j)=ly-0.5*Bmean-0.5*Bd(i-1,j);
        else
            TVl(i,j)=ly-0.5*(Bd(i,j)+Bd(i-1,j));
        end
    end
end
% horizontal throats
for i=1:ny
    for j=1:nx-1
        THl(i,j)=lx-0.5*(Bd(i,j)+Bd(i,j+1));
    end
end

% computing positions
% pore bodies
for i=1:ny
    for j=1:nx
        Bxp(i,j)=(j-1)*lx;
        Byp(i,j)=(ny-i+1)*ly;
    end
end
% horizontal pore throats
for i=1:ny
    for j=1:(nx-1)
        THxp(i,j)=((Bxp(i,j)+0.5*Bd(i,j))+(Bxp(i,j+1)-0.5*Bd(i,j+1)))/2;
        THyp(i,j)=Byp(i,j);
    end
end
% vertical pore throats
for i=1:(ny+1)
    for j=1:nx
        if i==1
            % top boundary
            TVxp(i,j)=Bxp(i,j);
            TVyp(i,j)=Byp(i,j)+0.5*TVl(i,j)+0.5*Bd(i,j);
        elseif i==(ny+1)
            % bottom boundary
            TVxp(i,j)=Bxp(i-1,j);
            TVyp(i,j)=Byp(i-1,j)-0.5*TVl(i,j)-0.5*Bd(i-1,j);
        else
            % inner
            TVxp(i,j)=Bxp(i,j);
            TVyp(i,j)=((Byp(i,j)+0.5*Bd(i,j))+(Byp(i-1,j)-0.5*Bd(i-1,j)))/2;
        end
    end
end

% plotting the network
clf;
subplot(2,2,1);
title('Geometry of Pore-Network System');
hold on;
axis equal;
% pore bodies
for i=1:ny
    for j=1:nx
        square(Bxp(i,j),Byp(i,j),Bd(i,j));
    end
end
% horizontal pore throats
for i=1:ny
    for j=1:(nx-1)
        rect(THxp(i,j),THyp(i,j),THl(i,j),THd(i,j))
    end
end
% vertical pore throats
for i=1:(ny+1)
    for j=1:nx
        rect(TVxp(i,j),TVyp(i,j),TVd(i,j),TVl(i,j))
    end
end
% pressure reservoirs
rect(Bxp(1,1)+0.5*Lx,Byp(1,1)+1.5*ly-0.5*Bmean,Lx+2*Bmax,ly);
text(Bxp(1,1)+0.5*Lx,Byp(1,1)+1.5*ly-0.5*Bmean,'Pmax');
rect(Bxp(1,1)+0.5*Lx,Byp(ny,nx)-1.5*ly+0.5*Bmean,Lx+2*Bmax,ly);
text(Bxp(1,1)+0.5*Lx,Byp(ny,nx)-1.5*ly+0.5*Bmean,' Pmin');

% main part
% constructing A matrix and b vector for the linear system AX=b
% A=BA (coefficients vector), X=Bp (unknown vector), b=Bb (left side vector)

% inner points
for i=2:(ny-1)
    for j=2:(nx-1)
      k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=(THd(i,j-1)^4)/THl(i,j-1);
      ar=(THd(i,j)^4)/THl(i,j);
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      BA(k,k+1)=ar;
      BA(k,k-nx)=at;
      BA(k,k+nx)=ad;
      Bb(k)=0;
      
    end
end

% top edge (including corners)
i=1;
for j=1:nx
    % top left corner
    if j==1
      k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=0;
      ar=(THd(i,j)^4)/THl(i,j);
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      % (no k-1 index)
      BA(k,k+1)=ar;
      % (no k-nx index)
      BA(k,k+nx)=ad;
      Bb(k)=-at*Pmax;
      
    % top right corner
    elseif j==nx
         k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=(THd(i,j-1)^4)/THl(i,j-1);
      ar=0;
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      BA(k,k+1)=ar;
      % (no k-nx index)
      BA(k,k+nx)=ad;
      Bb(k)=-at*Pmax;
      
    % top edge
    else
      k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=(THd(i,j-1)^4)/THl(i,j-1);
      ar=(THd(i,j)^4)/THl(i,j);
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      BA(k,k+1)=ar;
      % (no k-nx index)
      BA(k,k+nx)=ad;
      Bb(k)=-at*Pmax;
      
    end
end

% bottom edge (including corners)
i=ny;
for j=1:nx
    % bottom left corner
    if j==1
      k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=0;
      ar=(THd(i,j)^4)/THl(i,j);
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      BA(k,k+1)=ar;
      BA(k,k-nx)=at;
      % (no k+nx index)
      Bb(k)=-ad*Pmin;
      
    % bottom right corner
    elseif j==nx
      k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=(THd(i,j-1)^4)/THl(i,j-1);
      ar=0;
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      % (no k+1 index)
      BA(k,k-nx)=at;
      % (no k+nx index)
      Bb(k)=-ad*Pmin;
      
    % bottom edge
    else
      k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=(THd(i,j-1)^4)/THl(i,j-1);
      ar=(THd(i,j)^4)/THl(i,j);
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      BA(k,k+1)=ar;
      BA(k,k-nx)=at;
      % (no k+nx index)
      Bb(k)=-ad*Pmin;
      
    end
end

% left edge
j=1;
for i=2:(ny-1)
      k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=0;
      ar=(THd(i,j)^4)/THl(i,j);
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      BA(k,k+1)=ar;
      BA(k,k-nx)=at;
      BA(k,k+nx)=ad;
      Bb(k)=0;
      
end

% right edge
j=nx;
for i=2:(ny-1)
    k=(i-1)*nx+j;
      
      at=(TVd(i,j)^4)/TVl(i,j);
      ad=(TVd(i+1,j)^4)/TVl(i+1,j);
      al=(THd(i,j-1)^4)/THl(i,j-1);
      ar=0;
      ap=-(at+ad+al+ar);
      
      BA(k,k)=ap;
      BA(k,k-1)=al;
      BA(k,k+1)=ar;
      BA(k,k-nx)=at;
      BA(k,k+nx)=ad;
      Bb(k)=0;
      
end

% solving the linear system (BA.BP=Bb)
% determining the pore body pressures
% method is BiCGStab
tol = 1e-10; maxit = 10000;
BAsp=sparse(BA);
setup.type='ilutp';
setup.milu='row';
[M1,M2]=ilu(BAsp,setup);
BP=bicgstab(BAsp,Bb,tol,maxit,M1,M2);
P=VtoM(BP,nB,ny,nx);

% post-process
axis equal;
cn=20;
subplot(2,2,[2,4]);
contourf(Bxp,Byp,P,cn,'--');
title('Pore Bodies Pressure Field');
cbl = colorbar;
cbl.Label.String = 'Pressure (Pa)';

Pavg=zeros(ny,1);
lev=zeros(ny,1);
for i=1:ny
    lev(i,1)=ny+1-i;
    sumPV=0;
    sumV=0;
    for j=1:nx
        sumPV=sumPV+P(i,j)*(Bd(i,j)^3);
        sumV=sumV+(Bd(i,j)^3);
    end
    Pavg(i)=sumPV/sumV;
end
subplot(2,2,3);
scatter(Pavg,lev,'filled');
axis([Pmin,Pmax,1,ny]);
title('Average Pore Bodies Pressure in Each Level');

%
% END

% plotting the network
figure(2)
title('Geometry of Pore-Network System');
hold on;
axis equal;
% pore bodies
for i=1:ny
    for j=1:nx
        circle(Bxp(i,j),Byp(i,j),Bd(i,j));
    end
end
% horizontal pore throats
for i=1:ny
    for j=1:(nx-1)
        rect(THxp(i,j),THyp(i,j),THl(i,j),THd(i,j))
    end
end
% vertical pore throats
for i=1:(ny+1)
    for j=1:nx
        rect(TVxp(i,j),TVyp(i,j),TVd(i,j),TVl(i,j))
    end
end
% pressure reservoirs
rect(Bxp(1,1)+0.5*Lx,Byp(1,1)+1.5*ly-0.5*Bmean,Lx+2*Bmax,ly);
text(Bxp(1,1)+0.5*Lx,Byp(1,1)+1.5*ly-0.5*Bmean,'Pmax');
rect(Bxp(1,1)+0.5*Lx,Byp(ny,nx)-1.5*ly+0.5*Bmean,Lx+2*Bmax,ly);
text(Bxp(1,1)+0.5*Lx,Byp(ny,nx)-1.5*ly+0.5*Bmean,' Pmin');
