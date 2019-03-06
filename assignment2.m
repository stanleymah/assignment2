%% Part 1 A
% Part 1 is used to explore the representation of the electrostatic
% potential in a 2-D rectangular region. Beucase the potential will not
% change with respect to the y-coordinate, the problem will be treated as a
% 1-D problem.

clear all;
close all;

% Initial conditions
Length = 30;
Width = 20;
v0 = 1;

% Create matrices
N = Length*Width;
F = zeros(N,N);
G = zeros(N,1);

for i=1:N
    if (i == N) % right 
        F(i,i) = 1;
        G(i) = 0;
    elseif i == 1 % left
        F(i,i) = 1;
        G(i) = 1;
    else % Middle
        F(i,i)   = -2;
        F(i,i+1) =  1;
        F(i,i-1) =  1;
    end
end
    
V = F\G;

% Plot graph
figure(1);
plot(V);
xlabel('x');
ylabel('Voltage');
title('V(x)');

% The change is linear. As moving from left to right in the x direction
% shows a decrease in voltage. Voltage is uniform in the y diredtion. If
% this was mapped out on a voltage map, there would be a gradient in the
% x-direction with the same colour in the y-direction.

clear all

%% Part 1 B
% The finite differences method will be used to implement a matrix
% calculation, $GV=F$. $V$ represents the voltages at different points,
% $F$ is the matrix used to set the boundary conditions and G represents 
% the relation of voltages throughout.
% The the equation used:
% 
% $$\frac{V_{x-1,y}-2V_{x,y}+V_{x+1,y}}{(\Delta x)^2} + \frac{V_{x,y-1}-2V_{x,y}+V_{x,y+1}}{(\Delta y)^2}=0$$

dx = 0.25; % spacing along x
dy = 0.25; % spacing along y

Const1 = -2*(1/dx^2 + 1/dy^2);
Const2 = 1/(dx^2);
Const3 = 1/(dy^2);

% Initial conditions
Length = 30;
Width = 20;
v0 = 1;

% Create Grid
x = linspace(0,Length);
y = linspace(0,Width);

% Create matrices
N = Length*Width;
F = zeros(N,N);
G = zeros(N,1);

% Analytical solution
% Middle
for i = 2:Length-1
    for j = 2:Width-1
        n = i + (j-1)*Length;
        F(n,n) = Const1;
        F(n,n-1) = Const2;
        F(n,n+1) = Const2;
        F(n,n-Length) = Const3;
        F(n,n+Length) = Const3;
        G(n,1) = 0;
    end
end

% Left BC
i = 1;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Right BC 
i = Length;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

% Bottom BC
j = 1;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
end

% Top BC
j = Width;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
end

% To find the solution
v = F\G;

% converting for plot
for i = 1 : Length
    for j = 1 : Width
        n = i + (j-1)*Length;
        Ph(i,j) = v(n);
    end
end

% Plot
figure(2);
mesh(Ph);
xlabel('x');
ylabel('y');
zlabel('Voltage (V)');
title('Surface Plot of V - Analytic Solution');

% It can be seen that the potential is 1 V on one end and 0 V on the other.

% Numerical soliution
% Middle
for i = 2:Length-1
    for j = 2:Width-1
        n = i + (j-1)*Length;
        F(n,n) = -4;
        F(n,n-1) = 1;
        F(n,n+1) = 1;
        F(n,n-Length) = 1;
        F(n,n+Length) = 1;
        G(n,1) = 0;
    end
end

% Left BC
i = 1;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Right BC 
i = Length;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

% Top BC
j = Width;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

v = F\G;

for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*Length;
        Ph(i,j) = v(n);
    end
end

% Plot
figure(3);
mesh(Ph);
xlabel('x');
ylabel('y');
zlabel('Voltage (V)');
title('Surface Plot of V - Numeric Solution');

%% 
% The numeric solution has boundry conditions that is elevated on both
% ends with a potential at 1 V. The middle of the region contains the
% lowest potential as it is furthest from the influence of both ends. 

ph2 = Ph;
a = 30;
b = 10;
anew = 0;
for i = 1:Length
    for j = 1:Width
        for n = 1:2:1000
            anew = anew + ((1/n)*(cosh(n*pi*i/a))*(sin(n*pi*j/a))*(1/(cosh(n*pi*b/a))));
        end
    ph2(i,j) = (4/pi) * anew;
    end
end
figure(4);
surf(ph2);
xlabel('x');
ylabel('y');
zlabel('Voltage (V)');
title('Graphical Plot of V - Anaylytic Solution');

%%
% The graph above is a alternate representation of the potential, V using 
% the equation in the problem sheet. The analytic solution is not nearly as 
% accurate as the numeric solution. The graph does not properly show the 
% raised voltages on both ends of the region making the numeric solution 
% more accurate. 

%% Part 2 A
% Part 2 will provide an understanding of current flow in a 2D rectangular
% region by looking at the current density, electric potential, and electric 
% field with a bottleneck region. 

% Middle
for i = 2:Length-1
    for j = 2:Width-1
        n = i + (j-1)*Length;
        F(n,n) = -4;
        F(n,n-1) = 1;
        F(n,n+1) = 1;
        F(n,n-Length) = 1;
        F(n,n+Length) = 1;
        G(n,1) = 0;
    end
end

% Left BC
i = 1;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Right BC 
i = Length;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

% Top BC
j = Width;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

v = F\G;

for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*Length;
        Ph(i,j) = v(n);
    end
end

% Plot
figure(5);
mesh(Ph);
xlabel('X');
ylabel('Y');
zlabel('Voltage (V)');
title('Voltage Plot');

%% 
% Electric Field 
Elecmap = zeros(Length,Width);
for i = 1:Length
    for j = 1:Width
        n = j + (i-1)*Width;
        Elecmap(i,j) = v(n);
    end
end

for i = 1:Length
    for j = 1:Width
        if i == 1
            Ex(i,j) = (Elecmap(i+1,j) - Elecmap(i,j));
        elseif i == Length
            Ex(i,j) = (Elecmap(i,j) - Elecmap(i-1,j));
        else
            Ex(i,j) = (Elecmap(i+1,j) - Elecmap(i-1,j))*0.5;
        end
        if j == 1
            Ey(i,j) = (Elecmap(i,j+1) - Elecmap(i,j));
        elseif j == Width
            Ey(i,j) = (Elecmap(i,j) - Elecmap(i,j-1));
        else
            Ey(i,j) = (Elecmap(i,j+1) - Elecmap(i,j-1))*0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Et = Ex + Ey; % electric field plot generated by adding the x and y components

figure(6);
mesh(Et);
xlabel('X');
ylabel('Y');
zlabel('E');
title('Electric Field Plot');

% Sigma
sigma = ones(Length,Width);
for i = 1:Length
    for j = 1:Width
        if j <= (Width/3) || j >= (Width*2/3)
            if i >= (Length/3) && i <= (Length*2/3)
                sigma(i,j) = 10^-12;
            end
            
        end
    end
end

figure(7);
mesh(sigma);
xlabel('X');
ylabel('Y');
zlabel('Conduction');
title('Conductivity Map');
%%
% Figure 7 shows that the bottleneck area contain a conductivity lower
% than that of the surrounding area. This means the 2 boxes are areas of
% high resisitivity which will resist current flow. 

% Current Density
J = sigma .* Et;
figure(8);
mesh(J);
xlabel('X');
ylabel('Y');
zlabel('J');
title('Current Density Map');

%% Part 2 C

% Inputs
Length = 30;
Width = 20;

%Grid
x = linspace(0,Length);
y = linspace(0,Width);
dx = x(2) - x(1);
dy = y(2) - y(1);

% Make matrices
N = Length*Width;
F = zeros(N,N);
G = zeros(N,1);

% Middle
for i = 2:Length-1
    for j = 2:Width-1
        n = i + (j-1)*Length;
        F(n,n) = -4;
        F(n,n-1) = 1;
        F(n,n+1) = 1;
        F(n,n-Length) = 1;
        F(n,n+Length) = 1;
        G(n,1) = 0;
    end
end

% Left BC
i = 1;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Right BC 
i = Length;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

% Top BC
j = Width;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

v = F\G;

for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*Length;
        Ph(i,j) = v(n);
    end
end

% Sigma
sigma1 = ones(Length,Width);
sigma2 = ones(Length,Width);
sigma3 = ones(Length,Width);
sigma4 = ones(Length,Width);

for i = 1:Length %Changing the lengths and widths
    for j = 1:Width
        if j <= (Width/4) || j >= (Width*3/4)
            if i >= (Length/3) && i <= (Length*2/3)
                sigma1(i,j) = 10^-2;
            end 
        end
        if j <= (Width/3.1) || j >= (Width - (Width/3.1))
            if i >= (Length/3) && i <= (Length*2/3)
                sigma2(i,j) = 10^-2;
            end 
        end 
        if j <= (Width/2.5) || j >= (Width - (Width/2.5))
            if i >= (Length/3) && i <= (Length*2/3)
                sigma3(i,j) = 10^-2;
            end 
        end
        if j <= (Width/2.1) || j >= (Width - (Width/2.1))
            if i >= (Length/3) && i <= (Length*2/3)
                sigma4(i,j) = 10^-2;
            end 
        end
    end
end

t1 = sigma1;
t2 = sigma2;
t3 = sigma3;
t4 = sigma4;

for i = 1:Length
    for j = 1:Width
        if sigma1(i,j) == (10^-2)
            t1(i,j) = 1 / (10^-2);
        end
    end
end
for i = 1:Length
    for j = 1:Width
        if sigma2(i,j) == (10^-2)
            t2(i,j) = 1 / (10^-2);
        end
    end
end
for i = 1:Length
    for j = 1:Width
        if sigma3(i,j) == (10^-2)
            t3(i,j) = 1 / (10^-2);
        end
    end
end
for i = 1:Length
    for j = 1:Width
        if sigma4(i,j) == (10^-2)
            t4(i,j) = 1 / (10^-2);
        end
    end
end

Current1 = Ph ./ t1;
C01 = sum(Current1(1,:));
CL1 = sum(Current1(Length,:));
c1 = (C01 + CL1) / 2;
figure(9);
subplot(2,2,1);
mesh(Current1);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title('I');

Current2 = Ph ./ t2;
C02 = sum(Current2(1,:));
CL2 = sum(Current2(Length,:));
c2 = (C02 + CL2) / 2;
subplot(2,2,2);
mesh(Current2);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title('I');

Current3 = Ph ./ t3;
C03 = sum(Current3(1,:));
CL3 = sum(Current3(Length,:));
c3 = (C03 + CL3) / 2;
subplot(2,2,3);
mesh(Current3);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title('I');

Current4 = Ph ./ t4;
C04 = sum(Current4(1,:));
CL4 = sum(Current4(Length,:));
c4 = (C04 + CL4) / 2;
subplot(2,2,4);
mesh(Current4);
xlabel('X');
ylabel('Y');
zlabel('I');
title('I');
%%
% Narrowing the boxes should have caused the conductivity to decreased.

%% Part 2 D

% Middle
for i = 2:Length-1
    for j = 2:Width-1
        n = i + (j-1)*Length;
        F(n,n) = -4;
        F(n,n-1) = 1;
        F(n,n+1) = 1;
        F(n,n-Length) = 1;
        F(n,n+Length) = 1;
        G(n,1) = 0;
    end
end

% Left BC
i = 1;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Right BC 
i = Length;
for j = 1:Width
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 1;
end

% Bottom BC
j = 1;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

% Top BC
j = Width;
for i = 1:Length
    n = i + (j-1)*Length;
    F(n,n) = 1;
    G(n,1) = 0;
end

v = F\G;

for i = 1 : 30
    for j = 1 : 20
        n = i + (j-1)*Length;
        Ph(i,j) = v(n);
    end
end

% Sigma
sigma = ones(Length,Width);
for i = 1:Length
    for j = 1:Width
        if j <= (Width/3) || j >= (Width*2/3)
            if i >= (Length/3) && i <= (Length*2/3)
                sigma(i,j) = 10^-2;
            end
            
        end
    end
end

% Current Flow I = V/R
t1 = sigma;
t2 = sigma;
t3 = sigma;
t4 = sigma;

for i = 1:Length
    for j = 1:Width
        if sigma(i,j) == (10^-2)
            t1(i,j) = 1 / (10^-5);
            t2(i,j) = 1 / (10^0);
            t3(i,j) = 1 / (10^1);
            t4(i,j) = 1 / (10^5);
        end
    end
end

Cur = Ph ./ t1;
Ce = sum(Cur(1,:));
CL = sum(Cur(Length,:));
c = (Ce + CL) / 2;
figure(10);
subplot(2,2,1);
mesh(Cur);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title('S = 10^-5');

Cur = Ph ./ t2;
Ce = sum(Cur(1,:));
CL = sum(Cur(Length,:));
c = (Ce + CL) / 2;
subplot(2,2,2);
mesh(Cur);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title('S = 10^0');

Cur = Ph ./ t3;
Ce = sum(Cur(1,:));
CL = sum(Cur(Length,:));
c = (Ce + CL) / 2;
subplot(2,2,3);
mesh(Cur);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title('S = 10^1');

Cur = Ph ./ t4;
Ce = sum(Cur(1,:));
CL = sum(Cur(Length,:));
c = (Ce + CL) / 2;
subplot(2,2,4);
mesh(Cur);
xlabel('X');
ylabel('Y');
zlabel('I(x,y)');
title('S = 10^5');
%%
% As the box conductivity increased, current increased. 