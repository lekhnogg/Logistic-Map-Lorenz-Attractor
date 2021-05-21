%Assignment 3
%Gleb Lekhno 
%7776511


%Uncomment 17,41,42,43 Individually

function par = setup
clc

%%% ---> Logistic bifurcation <--- %%%
par.pNum = 100;                 % Number of points for each r
par.iterT = 100;                % Number of times to reach limit 
par.xx = zeros(par.pNum,1);     % Initialize empty zero vector
par.intV = 0.001;               % Value spacing

getBif(par) %<-- Uncomment to get bifurcation


%%% ---> Lorenz Attractor <--- %%%

% Initial positions
par.xInt = [-1;-1;1];
par.pPoints = 100;
% Parameters 
% General -> (1) 10,8/3,28 abstract ->(2) 20,2/3,17 (3) 19,11/3,32
par.sig = 10;
par.beta = 8/3;
par.rho = 28;
% Iteration steps
par.dt = 0.001;
% Maximum iteration time
par.tMax = 50;
% Iteration pan
par.tSpan = (par.dt:par.dt:par.tMax);
% Number of perturbations to compare
par.pertN = 20;
% Diameter of Gaussian perturbation ball
par.eps = 1e-12;

%simLorenz(par);        %<--Uncomment to get Lorentz Attractor
%divLorenz(par);        %<--Uncomment to see Exponential Divergence
%simLorenzBif(par);     %<--Uncomment to get Lorentz Bifurcation Diagram
end

function getBif(par)

% Iterate over each r value between 0 and 4 

for rr = 0:par.intV:4
    %iterate over arbitrary number of times
    
        % Using initial "population" of 0.5
        par.xx(1)= 0.5;
       
        % Iterate to reach limit
        for i = 1:par.iterT
        par.xx(1) = par.xx(1)*rr*(1-par.xx(1));
        end
        % Once value has "settled", we save the values to plot against
        % current "rr" 
        for i = 1:(par.pNum-1)
            par.xx(1+i) = rr*par.xx(i)*(1-par.xx(i));
        end
        hold on;
        % Plot array of same r value against xx
        plot(rr*ones(par.pNum,1),par.xx,'.k','markersize', 1.5);
end
title('Logistic function bifurcation diagram');
xlabel('r value');
ylabel('x value');
end

function getDisp = divLorenz(par)
% Retrieve parameters
vals = [par.sig; par.beta; par.rho];
% Establish column length for perturbed arrays
colLen = par.tMax/par.dt;
% Initialize vectors
xp = zeros(colLen,par.pertN);
yp = zeros(colLen,par.pertN);
zp = zeros(colLen,par.pertN);
% Store perturbed positions
pertPos = par.xInt;
% Loop over pertN iterations
for i = 1:par.pertN
    % Solve ODE
    options = odeset('RelTol',1e-13);
    [~,yy] = ode45(@(t,x)getLorenz(t,x,vals),par.tSpan,pertPos,options);
    
    % Perturb the initial positions with tiny Gaussian ball
    pertPos = par.xInt + (par.xInt.*(par.eps.*randn(3,1)));
    % Collect Perturbed Solutions
    xp(:,i) = yy(:,1);
    yp(:,i) = yy(:,2);
    zp(:,i) = yy(:,3);
end
% Get variances in positions
vx = var(xp.').';
vy = var(yp.').';
vz = var(zp.').';
vt = sqrt(vx+vy+vz);
% Take the log of variances
lvt = log(vt);

% Establish interval of exponential dispersion
intT = par.tSpan(18587);
endT = par.tSpan(40069);
% Fit a line to interval, return slope
polvt = lvt(18587:40069);
slope = polyfit(intT:par.dt:endT,polvt,1)

% Return plot of dispersion over time
getDisp = plot(par.tSpan,lvt,'.k','markersize',2);
title('Dispersion of trajectories with initial perturbations');
xlabel('Time (ms)');
ylabel('Log of position dispersion');
end


function pos = getLorenz(t,x,vals)
% Initialize empty array - stores x,y,z in that order

% Introduce differential equations
dx = vals(1)*(x(2)-x(1));
dy = x(1)*(vals(3)-x(3))-x(2);
dz = x(1)*x(2)-vals(2)*x(3);

% Pass column vector
pos = [dx; dy; dz];

end

function simLorenz(par)

vals = [par.sig; par.beta; par.rho];
% Use ode45 solver
options = odeset('RelTol',1e-13);
[~,yy] = ode45(@(t,x)getLorenz(t,x,vals),par.tSpan,par.xInt,options);

plot3(yy(:,1),yy(:,2),yy(:,3),'k');

end


function simLorenzBif(par)


hold on
eps = 1e-1;
% Get parameters
vals = [par.sig; par.beta; par.rho];
for rho = 0:0.25:250
    % Update plane crossing
    zz = rho-1;
    % Update rho value
    vals(3) = rho;
    options = odeset('RelTol',1e-13);
    % Use ode45 solver
    [~,yy] = ode45(@(t,x)getLorenz(t,x,vals),par.tSpan,par.xInt,options);
    % Estimate vector sizes
    yCross = zeros(length(yy((abs(yy(:,2))<eps))),1);
    zCross = zeros(length(yy((abs(yy(:,3))<eps))),1);
    % Begin iteration after trajectory sets
    for i = 1001:(length(yy(:,2))-1)
       if ((yy(i-1,3)-zz)*((yy(i+1,3))-zz)) < 0 && yy(i,1) ~= 0 && length(zCross) < 300
        zCross(i) = yy(i,1);
       end
       
       if ((yy(i-1,2))*((yy(i+1,2)))) < 0 && yy(i,1) ~= 0 
           % Collect crossings
        yCross(i) = yy(i,1);
       end
    end
    % Get only non zero contributions
    yCross = yCross(yCross~=0);
    zCross = zCross(zCross~=0);
    
    % Plot crossings per each rho
    plot(rho*ones(length(yCross),1),yCross,'.k','markersize', 1.25);
    plot(rho*ones(length(zCross),1),zCross,'.k','markersize', 1.25);
end
title('Lorenz Bifurcation Diagram');
xlabel('rho value');
ylabel('x value');

end

