clc;
clear;

N = 30;
M = 300;
H = 1;
L=10;
dx = L/ (M - 1);
dy = H/ (N - 1);
psi = zeros(M, N);
omega = zeros(M, N);
beta = dx / dy;
Re = 50;
u=zeros(M,N);
v=zeros(M,N);

% Boundary conditions
iterations = 1000; % Adjust the number of iterations as needed

for iter = 1:iterations
    
    for j=10:29
        psi(1,j+1)=psi(1,j)+(3*dy)/2*H; %baundry 1 %since U* is 1
        omega(1,j)=-(1/dy^2)*(psi(1,j+1)-2*psi(1,j)+psi(1,j-1));
        u(1,j)=1;
        v(1,j)=0;

    end   
    for i=2:30
        psi(i,10)=0; %baundry 2
        omega(i,10)=(-2/dy^2)*(psi(i,11)-psi(i,10));
        u(i,10)=0;
        v(i,10)=0;
    end
    
    for i=30:299
        psi(i,1)=0;  %baundry 3
        omega(i,1)=(-2/dy^2)*(psi(i,2)-psi(i,1));
        u(i,1)=0;
        v(i,1)=0;
    end
    for j= 1:30
        psi(300,j)=2*psi(299,j)-psi(298,j); %baundry 4
        omega(300,j)=omega(299,j);
        v(300,j)=v(299,j);
        u(300,j)=u(299,j);
    end
    
    for i=1:299
        psi(i,30)=1; %boundry 5
        omega(i,30)=(-2/dy^2)*(psi(i,29)-psi(i,30));
        u(i,30)=0;
        v(i,30)=0;
    end 
    for j=2:9
        psi(30,j)=0; %baundry 6
        omega(30,j)=(-2/dx^.2)*(psi(31,j)-psi(30,j));
        u(30,j)=0;
        v(30,j)=0;
    end

      psi_old = psi(i, j);
      omega_old = omega(i, j);
      
    for j = 11:29
        for i = 2:30
            psi_old = psi(i, j);
            psi(i, j) = (1 / (2 * (beta^2 + 1))) * (psi(i + 1, j) + psi(i - 1, j) + beta^2 * (psi(i, j + 1) + psi(i, j - 1)) + dx^2 * omega(i, j));
     
            if abs(psi(i, j) - psi_old) < 1e-6
                break;
            end
        end
    end   


    for j = 2:29
        for i = 31:299
            psi_old = psi(i, j);
            psi(i, j) = (1 / (2 * (beta^2 + 1))) * (psi(i + 1, j) + psi(i - 1, j) + beta^2 * (psi(i, j + 1) + psi(i, j - 1)) + dx^2 * omega(i, j));
            % Check for convergence
            if abs(psi(i, j) - psi_old) < 1e-6
                break;
            end
        end
    end  
   
end
for iter = 1:iterations
    for j = 11:29
        for i = 2:30
            omega_old = omega(i, j);
            omega(i,j)=(.5/(1+beta^2))*((1-(psi(i,j+1)-psi(i,j-1))*beta*Re*.25)*omega(i+1,j) ...
            +(1+(psi(i,j+1)-psi(i,j-1))*beta*Re*.25)*omega(i-1,j)...
            +(1+(psi(i+1,j)-psi(i-1,j))*((Re*.25)/beta))*omega(i,j+1) ...
            +(1-(psi(i+1,j)-psi(i-1,j))*((Re*.25)/beta))*omega(i,j-1)*beta^2);
            
            % Check for convergence
            if abs(omega(i, j) - omega_old) < 1e-6
                break;
            end
        end
    end
     for j = 2:29
        for i = 31:299
            omega_old = omega(i, j);
            omega(i,j)=(.5/(1+beta^2))*((1-(psi(i,j+1)-psi(i,j-1))*beta*Re*.25)*omega(i+1,j) ...
            +(1+(psi(i,j+1)-psi(i,j-1))*beta*Re*.25)*omega(i-1,j)...
            +(1+(psi(i+1,j)-psi(i-1,j))*((Re*.25)/beta))*omega(i,j+1) ...
            +(1-(psi(i+1,j)-psi(i-1,j))*((Re*.25)/beta))*omega(i,j-1)*beta^2);
            
            % Check for convergence
            if abs(omega(i, j) - omega_old) < 1e-6
                break;
            end
        end
    end

end

for j = 2:29
    for i = 2:5:299
        u(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / (2 * dy); % Central difference approx.
        v(i, j) = -(psi(i + 1, j) - psi(i - 1, j)) / (2 * dx); % Central difference approx.
    end
end


% Make sure v at the last column (N) is the same as v at N-1
v(:, N) = v(:, N-1);

% Create a vector plo
figure(3)
quiver(u', v', 'AutoScaleFactor', 3);
title('Velocity Vectors');
xlabel('X-axis');
ylabel('Y-axis');
hold off

% Define grid points for plotting
[Y,X] = meshgrid(linspace(0, 1, N), linspace(0, 10, M));

% Create a figure for streamlines
figure;
contour(X,Y,psi, 20); % Adjust the number of contours as needed
title('Streamlines');
xlabel('X');
ylabel('Y');
axis equal;
colorbar;

% Create a figure for vorticity contours
figure;
contour(X, Y, omega, 20); % Adjust the number of contours as needed
title('Vorticity Contours');
xlabel('X');
ylabel('Y');
axis equal;
colorbar;


