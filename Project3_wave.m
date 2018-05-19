
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ghaida El-Saied
%HW4 Finite Element Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
%%% subfunctions defined/from lecture notes


%% Define parameters and boundary conditions
tf = 6; 
dt = 0.005;
theta = 0;
Nt = tf/dt; 
ts = tf/dt/100 ; 
%mass lumping condition
masslumping = 1;
%% Mesh loading and data loading
coordinates = load('dsg-coordinates.dat');
elements =load('dsg-connectivity.dat');
elemD = load('dsg-dirichlet.dat'); 
elemN = load('dsg-neumann.dat');
ne = size(elements,1); 
nn = size(coordinates,1);
%% matrices
M = sparse(nn,nn); 
K = sparse(nn,nn); 
dof = setdiff(1:size(coordinates,1),unique(elemD)); % non-dirichlet nodes
DOF = [dof dof+nn];    
%% ASSEMBLY
if masslumping == 1
    M = diag(sum(M,2));
end

for e = 1:ne
    nodeIDs = elements(e,:);
    vert = coordinates(nodeIDs,:);
    M(nodeIDs,nodeIDs) = M(nodeIDs,nodeIDs) + mlocal(vert);
    K(nodeIDs,nodeIDs) = K(nodeIDs,nodeIDs) + klocal(vert);
end

%from lhsm AND rhsm
LHSM = [M, -theta*dt*M; theta*dt*K, M];
RHSM = [M, (1-theta)*dt*M; (theta-1)*dt*K, M];
dout = zeros(nn,round(Nt/ts)+1);
d = zeros(nn,1);
v = zeros(nn,1);

to = 1;
tic 
 for ti = 1:Nt
    F = zeros(2*nn,1);
    f_theta = zeros(nn,1);
    F = F + RHSM * [d;v];
    
    %specific consitions
    
    if ti >= Nt/2 && tf > 3
       u = zeros(nn,1);
       udot = zeros(nn,1);
       u2dot = zeros(nn,1);
    else
        u = zeros(nn,1);
        u(unique(elemD)) = theta* g(coordinates(unique(elemD),:),ti*dt) + ...
            (1-theta)*g(coordinates(unique(elemD),:),(ti-1)*dt);
       
        %first D
        
        udot=zeros(nn,1);
        udot(unique(elemD)) = theta * gdot(coordinates(unique(elemD),:),ti*dt) + ...
            (1-theta) * gdot(coordinates(unique(elemD),:),(ti-1)*dt);
        
        %second D
        
        u2dot = zeros(nn,1); 
        u2dot(unique(elemD)) = theta * g2dot(coordinates(unique(elemD),:),ti*dt) + ...
            (1-theta) * g2dot(coordinates(unique(elemD),:),(ti-1)*dt);
    
    end
    f_theta = - dt*(M*u2dot + K*u);
    F = F + [zeros(nn,1);f_theta];
    solution = [u;udot]; 
    
 %%% SOLVE 
  solution(DOF) = LHSM(DOF,DOF)\F(DOF);
  d = solution(1:nn);
  v = solution(nn+1:2*nn);
 
    if ~mod(ti,ts)
        to = to + 1;
        dout(:,to) = d;
    end
    if ~mod(ti,Nt/100)
        string = sprintf('%d percent complete',100*ti/Nt); 
        disp(string); 
    end
end
toc 
 %% results
trisurf(elements,coordinates(:,1),coordinates(:,2),dout(:,end))
xlim([-0.5 1]);ylim([0 1]);zlim([-0.2 0.2])
l = axis;
Frame(Nt) = struct('cdata',[],'colormap',[]);
vid = VideoWriter('Video_Name');
vid.FrameRate = 15;
vid.Quality = 75;
open(vid);
for ti = 1:(round(Nt/ts)+1)
 
    trisurf(elements,coordinates(:,1),coordinates(:,2),dout(:,ti))
    view(30,60);
    axis(l);
    drawnow;
    Frame(ti) = getframe;
    writeVideo(vid,Frame(ti));
end
close(vid);
function k = klocal(X)
kappa = 1; 
B = [-1 1 0; -1 0 1]; 
detJ = det([1,1,1; X']); 
Jinv = (1/detJ)*[X(3,2)-X(1,2) X(1,1)-X(3,1); ...
    X(1,2)-X(2,2) X(2,1)-X(1,1)]; 
k = kappa * 0.5 * (B') * Jinv * (Jinv') * B * detJ;
end
function m = mlocal(X)

J = det([1 1 1; X']); 
m = J * ... 
    [2 1 1; 1 2 1; 1 1 2]/24; 
end

function u = g(x,t)
A = 0.05;
w = 3;
u = A*sin(2*pi*w*t);
end

function udot = gdot(x,t)
A = 0.05;
w=3;
udot = 2*pi*w*A*cos(2*pi*w*t);
end

function u2dot = g2dot(x,t)
A = 0.05;
w=3;
u2dot=-4*pi*pi*w*w*A*sin(2*pi*w*t);
end

