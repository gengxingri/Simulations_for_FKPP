%FKPP: u_t = D*(u_{xx}+u_{yy}) + gamma*q(u) where q(u)='u.*(1-u)';%
%-----------------------------------------------------------------%
clc; clear;
tic;
n =200; 
time=100; %time steps 100
D=1; gamma=1; 

%Grid
eps=1; delta_t=0.2;
% eps=0.3;gamma_y=1;kappa=-1.9718;delta_t=0.0005;

h = 200/n; %epsilon=eps*1i/2;
x=linspace(-100, 100, n);
dt = delta_t;

%x and y meshgrid
y=x';
[xx,yy]=meshgrid(x,y);


%initial conditions
exp_mat=exp(-(xx.^2+yy.^2)/(2*eps));
u=1/sqrt(pi*eps)*exp_mat;

%inital condition
%u=0.1*ones(100);
%u = padarray(u, [50 50]);

grad=u*0; 
% Vectorization/index for u(i,j) and the loop --------
I = 2:n-1; J = 2:n-1;    


% ---- Time stepping ---------------------------------
for step=1:200          
% Laplace gradient of the equation     
 grad(I,J)= u(I,J-1)+u(I,J+1)+u(I-1,J)+u(I+1,J);
 u =(1-0.8*D)*u+0.2*D*grad+gamma*dt*u.*(1-u);   
%Plots each timestep
 meshc(u); 
 title(['Time ',num2str(step)]); %axis([0 N 0 N 0 0.4]); 
 xlabel x; ylabel y; zlabel u;
 view(43,22); drawnow;  
end
% ----- Topology of the final surface ----------------
surf(u);                               
shading interp; colormap jet;          
view([-25 70]); 

toc;