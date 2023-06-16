%% Supersonic flow over a flat plate
clear all
close all

%%%%%%%%%%%%%%%%%%%%%Define parameters%%%%%%%%%%%%%%%%%%%%%
mach=10; %mach number
Pr=0.71; %prandtl number
cp=1005; %cp
cv=718; %cv
R=cp-cv; %gas constant
s1=110.4; %sutherland temperature
soundspeed0=340.18; %m/s
u0=mach*soundspeed0; %free stream velocity(x direction)
v0=0; %free stream velocity(y direction)
p0=101300; %free stream pressure
T0=288.16; %free stream temperature
rho0=p0/(R*T0); %free stream density from p=rho*R*T
gamma=1.4; %specific heat ratio
mu0=1.735*10^-5; %initial viscosity
k0=(cp/Pr)*mu0; %initial thermal conductivity
lengthhorzi=10^-5; %horizontal length of domain(m)
lengthverti=8*10^-6; %vertical length of domain(m)


%%%%%%%%%%%%%%%%%%%%%Setting up grids%%%%%%%%%%%%%%%%%%%%%
nx=75; 
ny=80;
x= linspace(0,lengthhorzi,nx);
y= linspace(0,lengthverti,ny);
dx= x(2)-x(1);
dy= y(2)-y(1);
[xx,yy]=ndgrid(x,y);

%Initializing
tstep= 2.35*10^-11; %timestep
totaltime = 3525*10^-11;  %total time(this value=1500*tstep)
tvec=0:tstep:totaltime;%time vector
U=prim2cons(rho0,u0,v0,T0,cv);%Initializing U (all the conservative variables) this is a 4x75x80 matrix
[rho,u,v,T,p,e,Et]=cons2prim(U,R,cv);%Initializing primitive variables (rho,u,v,T,p,e,Et)

E=zeros(size(U));%Preallocate E
Ebar=zeros(size(U));%Preallocate Ebar
F=zeros(size(U));%Preallocate F
Fbar=zeros(size(U));%Preallocate Fbar
Ubar=zeros(size(U));%Preallocate Ubar
uconvervec=zeros(1,length(tvec));%Preallocating the convergence vector

%Preallocating gradients in the shear stresses
%In E
dudxinE=zeros(size(U,[2 3]));
dudyinE=zeros(size(U,[2 3]));
dvdxinE=zeros(size(U,[2 3]));
dvdyinE=zeros(size(U,[2 3]));
dTdx=zeros(size(U,[2 3]));
%In F
dudxinF=zeros(size(U,[2 3]));
dudyinF=zeros(size(U,[2 3]));
dvdxinF=zeros(size(U,[2 3]));
dvdyinF=zeros(size(U,[2 3]));
dTdy=zeros(size(U,[2 3]));
%Preallocating shear stresses and heat q
tauxx=zeros(size(U,[2 3]));
tauxyinE=zeros(size(U,[2 3]));
tauxyinF=zeros(size(U,[2 3]));
tauyy=zeros(size(U,[2 3]));
qx=zeros(size(U,[2 3]));
qy=zeros(size(U,[2 3]));

%Calculate the initial physical parameters
mu=calmu(T);
k=(cp/Pr)*mu;


%%%%%%%%%%%%%%%%%%%%%%%%%LOOP%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:length(tvec)-1

%%%%%%%%%%%%PREDICTOR STEP%%%%%%%%%%%%(forward in x and y)
%Update mu and k
mu=calmu(T);
k=(cp/Pr)*mu;

%presaving u for convergence check
uold=u;

%Compute the gradients:dudx,dudy,dvdx,dvdy,dTdx and dTdy(6 gradients)
%note that all dds need to be opposite biased
%in E flux(E flux is forward in x, so the ddxs need to be bwd and the ddys need to be central) 
dudxinE=ddx_bwd(u,dx);
dudyinE=ddy_central(u,dy);
dvdxinE=ddx_bwd(v,dx);
dvdyinE=ddy_central(v,dy);
dTdx=ddx_bwd(T,dx);
%in F flux(F flux is forward in y, so the ddys need to be bwd and the ddxs need to be central)
dudxinF=ddx_central(u,dx);
dudyinF=ddy_bwd(u,dy);
dvdxinF=ddx_central(v,dx);
dvdyinF=ddy_bwd(v,dx);
dTdy=ddy_bwd(T,dy);

%Compute shear stresses and heat q
%in E flux 
tauxyinE=mu.*(dudyinE+dvdxinE);
tauxx=2*mu.*(dudxinE-(1/3)*(dudxinE+dvdyinE) );
qx=-k.*dTdx;
%in F flux
tauxyinF=mu.*(dudyinF+dvdxinF);
tauyy=2*mu.*(dvdyinF-(1/3)*(dudxinF+dvdyinF) );
qy=-k.*dTdy;


%Assemble E(note that prim variables are the following:rho,u,v,T,p,e,Et)
E(1,:,:)=rho.*u;
E(2,:,:)=rho.*u.*u+p-tauxx;
E(3,:,:)=rho.*u.*v-tauxyinE;
E(4,:,:)=(Et+p).*u-u.*tauxx-v.*tauxyinE+qx;
%Assemble F
F(1,:,:)=rho.*v;
F(2,:,:)=rho.*u.*v-tauxyinF;
F(3,:,:)=rho.*v.*v+p-tauyy;
F(4,:,:)=(Et+p).*v-v.*tauyy-u.*tauxyinF+qy;

%Compute Ubar
Ubar=U-tstep*ddx_fwdE(E,dx)-tstep*ddy_fwdF(F,dy); 

%Obtain primitive variables from Ubar
[rho,u,v,T,p,e,Et]=cons2prim(Ubar,R,cv); %rho,u,v,T,p,e,Et these variables are updated 

%Enforcing the boundary conditions(prim variables are: rho,u,v,T,p,e,Et)
%inflow
u(1,2:end)=u0; %setting u=u0 aka u freestream
v(1,2:end)=0; %setting v=0 
p(1,2:end)=p0; %setting p=p0 aka p freestream
T(1,2:end)=T0; %setting T=T0 aka T freestream
rho(1,2:end)=p0/(R*T0); %updating rho=pressure/(R*T)
e(1,2:end)=cv*T0;  %updating e=cv*T
%outflow
for o=2:size(U,3)-1
    u(end,o)= 2*u(end-1,o)-u(end-2,o); %setting u(nx,j)=2*u(nx-1,j)-u(nx-2,j)
end
for o=2:size(U,3)-1
    v(end,o)= 2*v(end-1,o)-v(end-2,o); %setting v(nx,j)=2*v(nx-1,j)-v(nx-2,j)
end
for o=2:size(U,3)-1
   T(end,o)=2*T(end-1,o)-T(end-2,o); %setting T(nx,j)=2*T(nx-1,j)-T(nx-2,j)
    e(end,o)=cv*T(end,o);%updating e=cv*T
end
for o=2:size(U,3)-1
    p(end,o)= 2*p(end-1,o)-p(end-2,o); %setting p(nx,j)=2*p(nx-1,j)-p(nx-2,j)
    rho(end,o)=p(end,o)./(R*T(end,o));%updating rho=pressure/(R*T)
end
%farfield
u(1:end,end)=u0; %setting u=u0 aka u freestream
v(1:end,end)=0; %setting v=0 
p(1:end,end)=p0; %setting p=p0 aka p freestream
T(1:end,end)=T0; %setting T=T0 aka T freestream
rho(1:end,end)=p0/(R*T0); %updating rho=pressure/(R*T)
e(1:end,end)=cv*T0; %updating e=cv*T
%wall
u(2:end,1)=0; %setting u=0
v(2:end,1)=0; %setting v=0 
for w=2:size(U,2)
    p(w,1)= 2*p(w,2)-p(w,3); %setting p(i,1)=2*p(i,2)-p(i,3)
    rho(w,1)=p(w,1)/(R*T0);%updating rho=pressure/(R*T)
end
T(2:end,1)=T0; %setting T=T0=Twall 
e(2:end,1)=cv*T0; %updating e=cv*T 
%origin point(leading edge)
u(1,1)=0; %setting u=0
v(1,1)=0; %setting v=0 
p(1,1)=p0; %setting p=p0 aka p freestream
T(1,1)=T0; %setting T=T0 aka T freestream
rho(1,1)=p0/(R*T0); %updating rho=pressure/(R*T)
e(1,1)=cv*T0; %updating e=cv*T


%%%%%%%%%%%%CORRECTOR STEP%%%%%%%%%%%%(backward in x and y)
%Update mu and k
mu=calmu(T);
k=(cp/Pr)*mu;

%Compute the gradients:dudx,dudy,dvdx,dvdy,dTdx and dTdy(6 gradients)
%note that all dds need to be opposite biased
%in E flux, so the ddxs need to be fwd and the ddys need to be central
dudxinE=ddx_fwd(u,dx);
dudyinE=ddy_central(u,dy);
dvdxinE=ddx_fwd(v,dx);
dvdyinE=ddy_central(v,dy);
dTdx=ddx_fwd(T,dx);
%in F flux, so the ddys need to be bwd and the ddxs need to be central
dudxinF=ddx_central(u,dx);
dudyinF=ddy_fwd(u,dy);
dvdxinF=ddx_central(v,dx);
dvdyinF=ddy_fwd(v,dx);
dTdy=ddy_fwd(T,dy);

%Compute shear stresses and heat q
%in E flux 
tauxyinE=mu.*(dudyinE+dvdxinE);
tauxx=2*mu.*(dudxinE-(1/3)*(dudxinE+dvdyinE) );
qx=-k.*dTdx;
%in F flux
tauxyinF=mu.*(dudyinF+dvdxinF);
tauyy=2*mu.*(dvdyinF-(1/3)*(dudxinF+dvdyinF) );
qy=-k.*dTdy;


%Assemble E(note that prim variables are:rho,u,v,T,p,e,Et)
E(1,:,:)=rho.*u;
E(2,:,:)=rho.*u.*u+p-tauxx;
E(3,:,:)=rho.*u.*v-tauxyinE;
E(4,:,:)=(Et+p).*u-u.*tauxx-v.*tauxyinE+qx;
%Assemble F
F(1,:,:)=rho.*v;
F(2,:,:)=rho.*u.*v-tauxyinF;
F(3,:,:)=rho.*v.*v+p-tauyy;
F(4,:,:)=(Et+p).*v-v.*tauyy-u.*tauxyinF+qy;

%Compute U at the next time step
U=0.5*(U+Ubar-tstep*ddx_bwdE(E,dx)-tstep*ddy_bwdF(F,dy));

%Obtain primitive variables at the next time step(prim) from U
[rho,u,v,T,p,e,Et]=cons2prim(U,R,cv);

%Enforcing the boundary conditions(prim variables are: rho,u,v,T,p,e,Et)
%inflow
u(1,2:end)=u0; %setting u=u0 aka u freestream
v(1,2:end)=0; %setting v=0 
p(1,2:end)=p0; %setting p=p0 aka p freestream
T(1,2:end)=T0; %setting T=T0 aka T freestream
rho(1,2:end)=p0/(R*T0); %updating rho=pressure/(R*T)
e(1,2:end)=cv*T0;  %updating e=cv*T
%outflow
for o=2:size(U,3)-1
    u(end,o)= 2*u(end-1,o)-u(end-2,o); %setting u(nx,j)=2*u(nx-1,j)-u(nx-2,j)
end
for o=2:size(U,3)-1
    v(end,o)= 2*v(end-1,o)-v(end-2,o); %setting v(nx,j)=2*v(nx-1,j)-v(nx-2,j)
end
for o=2:size(U,3)-1
   T(end,o)=2*T(end-1,o)-T(end-2,o); %setting T(nx,j)=2*T(nx-1,j)-T(nx-2,j)
    e(end,o)=cv*T(end,o);%updating e=cv*T
end
for o=2:size(U,3)-1
    p(end,o)= 2*p(end-1,o)-p(end-2,o); %setting p(nx,j)=2*p(nx-1,j)-p(nx-2,j)
    rho(end,o)=p(end,o)./(R*T(end,o));%updating rho=pressure/(R*T)
end
%farfield
u(1:end,end)=u0; %setting u=u0 aka u freestream
v(1:end,end)=0; %setting v=0 
p(1:end,end)=p0; %setting p=p0 aka p freestream
T(1:end,end)=T0; %setting T=T0 aka T freestream
rho(1:end,end)=p0/(R*T0); %updating rho=pressure/(R*T)
e(1:end,end)=cv*T0; %updating e=cv*T
%wall
u(2:end,1)=0; %setting u=0
v(2:end,1)=0; %setting v=0 
for w=2:size(U,2)
    p(w,1)= 2*p(w,2)-p(w,3); %setting p(i,1)=2*p(i,2)-p(i,3)
    rho(w,1)=p(w,1)/(R*T0);%updating rho=pressure/(R*T)
end
T(2:end,1)=T0; %setting T=T0=Twall 
e(2:end,1)=cv*T0; %updating e=cv*T 
%origin point(leading edge)
u(1,1)=0; %setting u=0
v(1,1)=0; %setting v=0 
p(1,1)=p0; %setting p=p0 aka p freestream
T(1,1)=T0; %setting T=T0 aka T freestream
rho(1,1)=p0/(R*T0); %updating rho=pressure/(R*T)
e(1,1)=cv*T0; %updating e=cv*T

%compute U from the boundary condition "enforced" prim variables
%(note that prim is rho,u,v,T,p,e,Et)
U(1,:,:)=rho;
U(2,:,:)=rho.*u;
U(3,:,:)=rho.*v;
U(4,:,:)=rho.*(cv.*T+0.5*(u.^2+v.^2));

%checking for convergence 
%remove 4 boundaries of u and uold(to avoid dividing by zero)
unob=u(2:end-1,2:end-1);%u with no boundaries(no zeros)
uoldnob=uold(2:end-1,2:end-1);%u at previous time step with no boundaries(no zeros)
uconver=(mean((unob-uoldnob),"all"))./uoldnob; %take the mean of the difference and divide by uoldnob
uconvervec(1,t)=mean(abs(uconver),"all");%insert the mean value of the matrix into a vector, we will plot this vector later 

%saving T and p matrix for later use
Tcase1=T;
pcase1=p;
save ('midterm','Tcase1', 'pcase1')

%Numerical Schlieren
beta=0.8;
kappa=10;
gradrho=sqrt( (ddx_central(rho,dx)).^2+(ddy_central(rho,dy)).^2  );
S=beta*exp(-kappa*( gradrho ./ max(gradrho,[],"all") )  );

%Mach angle
mangle=asind(1/mach);%this is the analytical mach angle
pfixedaty1=p(:,29); %pressure vector as a function of x at a fixed y1 value
pfixedaty2=p(:,35); %pressure vector as a function of x at a fixed y2 value
dpfixedaty1=ddx_central(pfixedaty1,dx); % the derivative of pressure vector as a function of x at a fixed y1 value
dpfixedaty2=ddx_central(pfixedaty2,dx); % the derivative of pressure vector as a function of x at a fixed y2 value

t  %time count to keep track of runtime in the command window
    
%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%
if mod(t,100)==0 %plot when remainder of t/100 is 0
    figure(1)
    subplot(2,3,1)
    
    
    %rho field
    pcolor(xx,yy,rho)
    title('Density');  
    xlabel('X')
    ylabel('Y')
    colorbar
    shading interp
   
    %u field
    subplot(2,3,2)
    pcolor(xx,yy,u)
    title('u');  
    xlabel('X')
    ylabel('Y')
    colorbar
    shading interp

    %v field
    subplot(2,3,3)
    pcolor(xx,yy,v)
    title('v');  
    xlabel('X')
    ylabel('Y')
    colorbar
    shading interp

    %internal energy
    subplot(2,3,4)
    pcolor(xx,yy,e)
    title('Internal energy');  
    xlabel('X')
    ylabel('Y')
    colorbar
    shading interp

    %pressure 
    subplot(2,3,5)
    pcolor(xx,yy,p)
    title('Pressure');  
    xlabel('X')
    ylabel('Y')
    colorbar
    shading interp

    %temperature
    subplot(2,3,6)
    pcolor(xx,yy,T)
    title('Temperature');  
    xlabel('X')
    ylabel('Y')
    colorbar
    shading interp
    
    
    drawnow
end
end

%plotting the convergence plot
figure(3)
plot(tvec,uconvervec)
title('Convergence of u(velocity)');  
xlabel('time')
ylabel('convergence of u')

%% all functions


function r=ddx_fwdE(f,dx)%forward difference for E(x direction)
r=zeros(size(f));
for e=1:size(f,1)
r(e,:,:)=ddx_fwd(squeeze(f(e,:,:)),dx);
end
end

function r=ddy_fwdF(f,dy)%forward difference for F(y direction)
r=zeros(size(f));
for e=1:size(f,1)
r(e,:,:)=ddy_fwd(squeeze(f(e,:,:)),dy);
end
end

function r=ddx_bwdE(f,dx)%backward difference for E(x direction)
r=zeros(size(f));
for e=1:size(f,1)
r(e,:,:)=ddx_bwd(squeeze(f(e,:,:)),dx);
end
end

function r=ddy_bwdF(f,dy)%backward difference for F(y direction)
r=zeros(size(f));
for e=1:size(f,1)
r(e,:,:)=ddy_bwd(squeeze(f(e,:,:)),dy);
end
end

%below are functions from previous homeworks
function r=ddx_fwd(f,dx)
r=zeros(size(f));
  for i=1:size(f,2) %fixing the y direction(the columns)
    for j=1:size(f,1)-1 %do operation on the x direction(on every row)
      r(j,i)=( -f(j,i)+f(j+1,i) ) /dx;
    end
    %use bwd to calculate the final point, i.e., 50th point
    r(end,i)=( f(end,i)-f(end-1,i) ) /dx;
  end

end

function r=ddx_bwd(f,dx)
r=zeros(size(f));
  for i=1:size(f,2) %fixing the y direction(the columns)
    for j=2:size(f,1) %do operation on the x direction(on every row)
      r(j,i)=( f(j,i)-f(j-1,i) ) /dx;
    end
    %use fwd to calculate the first point
     r(1,i)=( -f(1,i)+f(2,i) ) /dx;
  end
end


function r=ddx_central(f,dx)
r=zeros(size(f));
  for i=1:size(f,2) %fixing the y direction(the columns)
    for j=2:size(f,1)-1 %do operation on the x direction(on every row)
      r(j,i)=(f(j+1,i)-f(j-1,i))/(2*dx);
    end
    %use fwd & bwd to calculate the first and final point
    %point
      r(1,i)=( -f(1,i)+f(2,i) ) /dx;
   r(end,i)=( f(end,i)-f(end-1,i) ) /dx;
  end

end
    
function r=ddy_bwd(f,dy)
r=zeros(size(f));
  for i=1:size(f,1) %fixing the x direction(the rows)
    for j=2:size(f,2) %do operation on the y direction(on every columns)
      r(i,j)=(  f(i,j)-f(i,j-1))/dy;
    end
    %use fwd to calculate first point
     r(i,1)=(  -f(i,1)+f(i,2) )/dy;
  end
end

function r=ddy_fwd(f,dy)
r=zeros(size(f));
  for i=1:size(f,1) %fixing the x direction(the rows)
    for j=1:size(f,2)-1 %do operation on the y direction(on every columns)
      r(i,j)=(  -f(i,j)+f(i,j+1))/dy;
    end
   %use bwd to calculate final point
     r(i,end)=(  f(i,end)-f(i,end-1))/dy;
  end
end

function r=ddy_central(f,dy)
r=zeros(size(f));
  for i=1:size(f,1) %fixing the x direction(the rows)
    for j=2:size(f,2)-1 %do operation on the y direction(on every columns)
      r(i,j)=(f(i,j+1)-f(i,j-1))/(2*dy);
    end
    %use fwd & bwd to calculate the first and final point
    %point
    r(i,1)=(  -f(i,1)+f(i,2) )/dy;
    r(i,end)=(  f(i,end)-f(i,end-1))/dy;
  end
end

%the conversion function from conservative variables to prim variables
function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv) %rho,u,v,T,p,e,Et

rho     = squeeze(U(1,:,:));
u       = squeeze(U(2,:,:))./rho;
v       = squeeze(U(3,:,:))./rho;
Et      = squeeze(U(4,:,:));
e       = Et./rho - (u.^2 + v.^2)/2;
T       = e/cv;
p       = rho*R.*T;
end

%the conversion function from prim variables to conservative variables
function r=prim2cons(rho,u,v,T,cv)
r=zeros(4,75,80);
r(1,:,:)=rho; 
r(2,:,:)=rho.*u;
r(3,:,:)=rho.*v;
r(4,:,:)=rho.*(cv.*T+0.5.*(u.^2+v.^2));
end

%sutherland function
function r=calmu(T)
T0=280.16;
vis0=0.00001735;
s1=110.4;
r=vis0.*((T./T0).^1.5).*(   (T0+s1)./(T+s1)   );
end