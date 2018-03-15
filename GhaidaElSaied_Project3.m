%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ghaida El-Saied
%Project 3
%Finite Element Analysis, UC Berkeley
%credit to lecture 3 notes for h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;


%Load Mesh data
%where ne is the number of elements and nn is number nodes
elems = load('naca5012_connectivity.dat');
coords = load('naca5012_coordinates.dat');
ne = size(elems,1); 
nn = size(coords,1); 


% Load Boundary Mesh
%nN= neuman nodes and dN = drichlet nodes

dN = load('naca5012_airfoil_dirichlet.dat');
nN = load('naca5012_airfoil_neumann.dat');
right_dirichletNodes = load('naca5012_right_dirichlet.dat');
left_neumannNodes = load('naca5012_left_neumann.dat');



% FEM
% Initilization
K = sparse(nn,nn);
F = zeros(nn,1);

e = 1;
while e <= ne
    nodeIds = elems(e,:);
    vertices = coords(nodeIds,:);
    K(nodeIds,nodeIds) = K(nodeIds,nodeIds) + klocal(vertices);
    e = e+1;
end

% Assembling Neumann BC
for side = 1:size(left_neumannNodes,1)
    elementId = left_neumannNodes(side,:);
    if elementId(2) == 0
        msg = 'error';
        msg;
    elseif elementId(2) == 1
        nnIDs = elems(elementId(1),[2,3]);
    elseif elementId(2) == 2
      nnIDs = elems(elementId(1),[1,3]);
    elseif elementId(2) == 3
      nnIDs = elems(elementId(1),[1,2]);
    end
    verts = coords(nnIDs,:);
    F(nnIDs) = F(nnIDs) + norm(verts(1,:)-verts(2,:))/2*2*0.5*h(sum(verts)/2);
end


% Dirichlet conditions
u = zeros(nn,1);
u(unique(right_dirichletNodes)) = g(coords(unique(right_dirichletNodes),:));
F = F - K*u;

% computation of the solution
dof = setdiff (1: size (coords ,1) ,unique(right_dirichletNodes));
u(dof) = K(dof ,dof) \ F(dof);


% Quiver Plot
v = zeros(2,ne);
final_center = zeros(2,ne); 

e= 1;
while e <= ne
  NI = elems(e,:);
  vertices = coords(NI,:); X = vertices;
  e_center = sum(vertices)/3; 
  final_center(:,e) = e_center; 
  
  %computing element stiffness matrix
  
  B = [-1 1 0; -1 0 1];
  det_J = det ( [ 1 , 1 , 1 ; X' ] ) ; 
  Jinv = (1/ det_J ) *[X(3 ,2)-X(1 ,2), X(1 ,1)-X(3 ,1) ; ... 
  X(1 ,2)-X(2 ,2), X(2 ,1)-X(1 ,1) ] ;

  grad_N = (Jinv') * B;
  nodalPHIs = u(NI);
  velocitygrad = grad_N*nodalPHIs; 
  v(:,e) = velocitygrad;
  e = e+1;

end
x_val = final_center(1,:);
y_val = final_center(2,:);
vxvalues = v(1,:);
vyvalues = v(2,:);

%quiver plot of velocity field

quiver(x_val,y_val,vxvalues,vyvalues,0.5)
title('Velocity Field Plot', 'FontSize', 24)
xlabel('X')
ylabel('Y')
axis([-0.4 1.4 -0.8 0.8])

% Calculating pressure
e=1;
while e <=ne
    velocityElements = v(:,e);
    pressureElements(e) = 0.5-(0.5*(velocityElements(1)^2 + velocityElements(2)^2));
    e = e+1;
end
figure
colormap(jet)
patch('Faces', elems, 'Vertices', coords, 'CData', pressureElements, 'FaceColor', 'flat', 'LineStyle', 'none')
title('Gage Pressure Field', 'FontSize', 24)
xlabel('x')
ylabel('y')
axis([-0.4 1.4 -0.8 0.8])
airfoil_x = coords(dN, 1);
airfoil_y = coords(dN, 2);
k = boundary(airfoil_x, airfoil_y, 0);
patch(airfoil_x(k), airfoil_y(k), 'white')
colorbar

% Lift:
forc_El = zeros(3,size(nN,1));
for e = 1:size(nN,1)  
   foil_element = nN(e,:);
    if foil_element(2) == 1
        nnIDs = elems(foil_element(1),[2,3]);  
        verts = coords(nnIDs,:);
        D = verts(2,:)-verts(1,:);
        D = [D, 0];
    elseif foil_element(2) == 2
      nnIDs = elems(foil_element(1),[1,3]); 
      verts = coords(nnIDs,:);
      D = verts(1,:)- verts(2,:);
      D = [D, 0];
    elseif foil_element(2) == 3
      nnIDs = elems(foil_element(1),[1,2]);
      verts = coords(nnIDs,:);
      D = verts(2,:)-verts(1,:);
      D = [D, 0];
    end   
    t = D/norm(D);
    z_axis = [0, 0, 1];
    n = cross(z_axis,t);

% calculate force
forc_El(:,e) = norm(verts(1,:) - verts(2,:))/2 * ...
 2 * ...
 pressureElements((nN(e,1)))*(n');
end
lift_force_vector = sum(forc_El,2)
lift = norm(lift_force_vector)

% computation of element stiffness 
function k = klocal(X)
    B = [-1 1 0; -1 0 1];
    detJ = det([1, 1, 1;X']);
    Jinv = (1/detJ)* [X(3,2)-X(1,2), X(1,1)-X(3,1);...
                      X(1,2)-X(2,2), X(2,1)-X(1,1)];
    k = 0.5 * (B') * Jinv * (Jinv') * B * detJ;
end

%computation of element force vector
function f = flocal(X)
J=det([1,1,1; X']) ;
J;
end

%functions from lecture 3. 

function out = g(X)
    out = zeros(size(X,1),1);
end

function out = h(X)
    out = -1.*ones(size(X,1),1);
end

function source = s(X)
    source = ones(size(X,1),1);
end