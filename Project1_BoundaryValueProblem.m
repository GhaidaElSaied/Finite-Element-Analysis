%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ghaida El-Saied
%Project 1
%Finite Element Analysis
%University of California at Berkeley

% A program to solve the fh = -(k^2*sin((pi*k*xh)/L))/A - (2*xh)/A;
% Boundary value problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Boundary conditions
%L=length, A = amplitude, u0 and uL are first and second boundary
%conditions respectively

A = 0.2; L= 1; u0 = 0; uL = 1; 

%BVP
%fh = -(k^2*sin((pi*k*xh)/L))/A - (2*xh)/A;
%true solution
%u_t = (x/L) + (x/(3.*A)).*(x.^2 - L.^2) - (L.^2/(A.*pi.^2)).*sin((pi.*k.*x)/L); 


%k and ne (number of elements) values:
k = [1 2 4 8 16 32];
ne = [2 4 8 16 32 64 128 256 512 1024 2048]; 
test_array = [1 1 1 1 1 1 1 1 1 1 1];



%setting up array for plotting

count = 1; 
label = cell(1,count);
bar_x = zeros(1,6);
bar_y = zeros(1,6); 
i=1;



%displaying figure for k values - is outputted below
for k = [1 2 4 8 16 32]
    figure;       
    
    for ne = [2 4 8 16 32 64 128 256 512 1024 2048] % Ne values to loop around
        %Global Matrices: stiffness = P, Force = F, Displacement = D
        P = zeros(ne + 1,ne + 1); F = zeros(ne + 1, 1); D = zeros(ne + 1, 1);

       % create a mesh
        he = 1/ne;
        xh = 0:he:L;
        
        fh = -(k^2*sin((pi*k*xh)/L))/A - (2*xh)/A;
            e =1;
              while e <= ne %loops over the ne matrix
                stiffness_matrix = 1/he * [1 -1; -1 1]; 
                %discrete forces at each node:
                f_left = fh(e); f_right = fh(e+1);
                %force for each element:
                f_e = he/6 * [2 1;1 2] * [f_left; f_right]; 
                %creating force vectors
                if e == 1 
                   f_e = f_e - u0 * stiffness_matrix(:,1);
                end
                if e == ne
                    f_e = f_e - uL * stiffness_matrix(:,2); 
                end
                %combining stiffness and force elemennts
                P(e:e+1, e: e+1) = P(e:e+1, e: e+1) + stiffness_matrix;
                F(e:e+1) = F(e:e+1) + f_e;
                e = e + 1;
              end
              
    %Displacement
    D(1) = 0;
    D(ne + 1) = 1;
    D(2:ne) = P(2:ne, 2:ne)\F(2:ne);
    
   
    plot(xh, D,'Linewidth',3)
    hold on;

  %calculating error for each # of elements and k value  
    x = 0:L/1000:L;
    u_h = interp1(xh,D,x);
    u_t = (x/L) + (x/(3.*A)).*(x.^2 - L.^2) - (L.^2/(A.*pi.^2)).*sin((pi.*k.*x)/L);
    diff = (u_t - u_h).^2; error = trapz(x,diff);
    label{1,count} = ['Ne = ' num2str(ne) ', '  'error = ' num2str(error)];
    if error <= 0.01
        k_val = k; ne_val = ne;
        tolerance = error;
        bar_x(1,i) = k;
        bar_y(1,i) = ne; 
        count = 1;
        i = i+1;
        break;
    end
    count = count +1;
    end

    %plot for true solution:
    
    hold on;
    plot(x,u_t,'Linewidth',2)
    legend(label,'u^t true solution');
    legend('Location', 'NorthWestOutside');
end

figure;
cat = categorical(bar_x);
bar(cat,bar_y); 
grid on;
xlabel('K');ylabel('number of elements ne(k)');
title('K vs. ne(k)');
