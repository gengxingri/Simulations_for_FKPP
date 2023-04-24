%FKPP: u_t = D*(u_{xx}+u_{yy}) + gamma*q(u) where q(u)='u.*(1-u)';%
%     with the EBC: [u_y](x,0,t) = -2*a*u_{xx}(x,0,t)             %
%-----------------------------------------------------------------%
clc; clear;
tic;
n =200; 
d = n/2;
time=100; %time steps 100
D=1; gamma=1; 
a = 10;


%Grid
eps=0.2; delta_t=0.2; dt = delta_t;
% eps=0.3;gamma_y=1;kappa=-1.9718;delta_t=0.0005;
%eps = dt/h^2;
h = 200/n; 
x=linspace(-100, 100, n);


%x and y meshgrid
y=x';
[xx,yy]=meshgrid(x,y);


%initial conditions
exp_mat=exp(-(xx.^2+yy.^2)/(4*pi));
u=0.5*1/sqrt(4*pi)*exp_mat;

%initial grad
grad=u*0; 

% Vectorization/index for u(i,j) and the loop --------
% without road
%I = 2:n-1; J = 2:n-1;  

% ---- Time stepping ---------------------------------
for step=1:1500 
     % with a road
     Un = u;
     for i = 2 : n-1
         for j = 2: n-1
             if i ~= d
                 grad = Un(i-1,j) + Un(i+1,j)+ Un(i,j-1) + Un(i,j+1);

                 %FKPP
                 %u(i,j) = (1-4*eps)*Un(i,j) + eps*grad + dt*Un(i,j)*(1-Un(i,j));

                 %Fujita
                 u(i,j) = (1-4*eps)*Un(i,j) + eps*grad + dt*Un(i,j).^2.25;
             else
                 u(d, j) = (2*a*(Un(d,j+1) + Un(d,j-1)) + h*(Un(d+1,j) + Un(d-1,j)))/(4*a+2*h);
        
             end
         end
     end

% Laplace gradient of the equation 
% Without a road
%   grad(I,J)= u(I,J-1)+u(I,J+1)+u(I-1,J)+u(I+1,J);
%   u =(1-0.8*D)*u+0.2*D*grad+gamma*dt*u.^2;% Fujita equation


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