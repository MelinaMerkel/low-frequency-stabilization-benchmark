%-----------------------------------%
%------------ Settings -------------%
%-----------------------------------%
hmin = 1e-10;
hmax = 1e10;
N = 11; %Number of evaluation points

load 'low_frequency_stabilization_benchmark.mat'
%   K:       stiffness matrix
%   M:       mass matrix
%   idx_b:   indices of (Dirichlet) boundary nodes
%   idx_nb:  indices of non-boundary nodes
%   idx_c:   indices of conductor nodes (excluding boundary nodes)
%   idx_nc:  indices of non-conductor nodes (excluding boundary nodes)


% Material parameterss
eps0 = 8.8541878128e-12;
epsilon = 2*eps0;
kappa = 5.96e7;

% h-loop
h_range = logspace(log10(hmin), log10(hmax), N);

% initialize
conditions = zeros(N,7);

K11 = K(idx_c,idx_c); K12 = K(idx_c,idx_nc); K21 = K(idx_nc,idx_c); K22 = K(idx_nc,idx_nc);
M11 = M(idx_c,idx_c); M12 = M(idx_c,idx_nc); M21 = M(idx_nc,idx_c); M22 = M(idx_nc,idx_nc);


for i = 1:N
    h = h_range(i);
    
    for ii = 0:6
        
        switch ii
            case 0
                a1 = 1; a2 = 1; b1 = 1; b2 = 1;
                A = [a1*b1*K11, K12; K21, K22]+[1/h*a1*b1*M11, 1/h*a1*b2*M12; 1/h*a2*b1*M21, 1/h*a2*b2*M22];
            case 1
                a1 = 1; a2 = sqrt(h); b1 = 1; b2 = sqrt(h);
                A = [a1*b1*K11, K12; K21, K22]+[1/h*a1*b1*M11, 1/sqrt(h)*M12; 1/sqrt(h)*b1*M21, M22];
            case 2
                a1 = 1; a2 = h; b1 = 1; b2 = 1;
                A = [a1*b1*K11, K12; K21, K22]+[1/h*a1*b1*M11, 1/h*a1*b2*M12; b1*M21, b2*M22];
            case 3
                a1 = sqrt(1/(kappa+1/h*epsilon)); a2 = sqrt(1/epsilon); b1 = sqrt(1/kappa+1/h*epsilon); b2 = sqrt(1/epsilon);
                A = [1/(kappa+epsilon/h)*K11 K12; K21 K22] + [1/(kappa*h+epsilon)*M11, 1/sqrt(kappa*epsilon*h+epsilon^2)*M12; 1/sqrt(kappa*epsilon*h+epsilon^2)*M21, 1/epsilon*M22];
            case 4
                a1 = 1/(kappa+1/h*epsilon); a2 = 1/epsilon; b1 = 1; b2 = 1;
                A = [a1*b1*K11, K12; K21, K22]+[1/h*a1*b1*M11, 1/h*a1*b2*M12; 1/1*a2*b1*M21, 1/1*a2*b2*M22];
            case 5
                inva1 = (K11+1/h*M11); inva2 = (M22); b1 = 1; b2 = 1;
                A = [eye(size(K11)), K12; K21, K22]+[zeros(size(M11)), inva1\M12/h; inva2\M21,  eye(size(M22))];
            case 6
                inva1 = (K11); inva2 = (M22); b1 = 1; b2 = 1;
                A = [eye(size(K11)), K12; K21, K22]+[inva1\M11/h, inva1\M12/h; inva2\M21, eye(size(M22))];
        end
        
        
        A([idx_c,idx_nc],[idx_c,idx_nc]) = A;
        conditions(i,ii+1) = condest(transp(A(idx_nb, idx_nb)));
        
    end
end
close all
markerShapes = {'d', 'o', 's', '^', 'v','d', 'o', 's', '^', 'v','d', 'o', 's', '^', 'v'};
for k = 1:7
    loglog(h_range', conditions(:,k), 'Marker', markerShapes{k}, 'MarkerSize',6,'LineWidth',mod(7-k,3)+1)
    hold on
end
xlabel('Time step size'); ylabel('Condition number');
legend('original','Scaling (i)','Scaling (ii)','Scaling (iii)','Scaling (iv)','Scaling (v)','Scaling (vi)')
