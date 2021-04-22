clear all;                  % clear all variables

c = 1;                      % c value of pde
xs=0; xe=1;                 % x domain [xs,xe]
I = 4;                     % I: number of division for x
dx = (xe-xs) / I;           % dx: mesh size
tf = 1;                   % final simulation time
J = 36;                    % J: number of time steps
dt = tf/J;
k = 0;
t = 0;

lambd = dt*c*c/(dx)^2;

if lambd > 0.5         % make sure dt satisy stability condition
    error('lambda should < 0.5!')
end
% Evaluate the initial conditions
x = xs : dx : xe;              % generate the grid point

f = sin(pi*x);    % Initial condition

% store the solution at all grid points for all time steps
u = zeros(I+1,J);
u_ex = zeros(I+1,J) ;
% Find the approximate solution at each time step
for j = 1:J-1
    t = j*dt;         % current time
    
    % boundary condition at left side
    gl = 0;
    % boundary condition at right side
    gr = 0;
    if j==1    % first time step
        for i=2:I    % interior nodes using schmidt method
            u(i,j) = lambd*f(i+1)-(1-2*lambd)*(f(i)+f(i))+lambd*f(i-1);
        end
        u(1,j) = gl;   % the left-end point
        u(I+1,j) = gr; % the right-end point
    else
        for i=2:I    % interior nodes
            u(i,j)=((1-2*lambd)/(1+2*lambd))*u(i,j-1)+lambd*2*(u(i+1,j)+u(i-1,j))/(1+2*lambd);
        end
        u(1,j) = gl;   % the left-end point
        u(I+1,j) = gr; % the right-end point
    end
    
    % calculate the analytic solution
    
    for i=1:I+1
        xj = xs + (i-1)*dx;
        u_ex(i,j)=sin(pi*i)*exp(-pi*pi*t);
    end
end

% Plot the results
tt = dt : dt : J*dt;
figure(1)
colormap(gray);     % draw gray figure
surf(x,tt, u');     % 3-D surface plot
xlabel('x')
ylabel('t')
zlabel('u')
title('Numerical solution of 1-D parabolic equation')

figure(2)
surf(x,tt, u_ex');     % 3-D surface plot
xlabel('x')
ylabel('t')
zlabel('u')
title('Analytic solution of 1-D parabolic equation')

maxerr=max(max(abs(u-u_ex)));