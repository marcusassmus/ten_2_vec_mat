%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this file, the tensor notation of the symmetric stiffness of Hooke's 
% law T = C : E is transferred into a 6 x 6 matrix notation w.r.t an ONB, 
%
% i.e. from
%
% C = C_{ijkl} e_{i} dyad e_{j} dyad e_{k} dyad e_{l}  
%                                               for all i,j,k,l in {1,2,3}
%
% to
%
% C = [ C_{11} C_{12} C_{13} C_{14} C_{15} C_{16};
%       C_{12} C_{22} C_{23} C_{24} C_{25} C_{26}; 
%       C_{13} C_{23} C_{33} C_{34} C_{35} C_{36};
%       C_{14} C_{24} C_{34} C_{44} C_{45} C_{46};
%       C_{15} C_{25} C_{35} C_{45} C_{55} C_{56};
%       C_{16} C_{26} C_{36} C_{46} C_{56} C_{66} ];
%   = [C_{alpha beta}]                 for all alpha, beta in {1,2,3,4,5,6}
% . 
% We hereby employ two different methodologies to relate the material 
% parameters C_{ijkl} and C_{alpha beta}. 
%
% Please see ten_2_vec-mat.pdf for a comprehensive description of the 
% problem.
%
% Correspondence should be adressed to Marcus A{\ss}mus
%
% Work address :    Chair of Engineering Mechanics
%                   Institute of Mechanics
%                   Faculty of Mechanical Engineering
%                   Otto von Guericke University
%                   Universitaetsplatz 2                    
%                   39108 Magdeburg 
%                   Germany
% Website:          https://www.ifme.ovgu.de/ltm
% E-Mail:           marcus.assmus@ovgu.de
%                   
%
% June 2020; Last revision: dd-month-yyyy
%
% When refering to the script in publications please cite as follows:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

I = eye(3); 
% some abbreviations for normalization factors
v = (1./2); % simply symmetric part of shear terms [Voigt]
k = (sqrt(2)./2); % normalized symmetric part of shear terms [Kelvin]

%%
% test stuff

% [x,y,z,m]=ndgrid(1:3,1:3,1:3,1:3);
% C = sym( reshape( ...
%       cellfun(@sym, ...
%         strcat('C_', ...
%           cellstr(num2str(z(:))), ...
%           cellstr(num2str(m(:))), ...
%           cellstr(num2str(x(:))), ...
%           cellstr(num2str(y(:)))), ...
%       'Uniform', 0), ...
%       3, 3, 3, 3 ) );
%%
% stiffness tensor with major, left minor, and right minor symmetry
%
C = sym( zeros(3,3,3,3) );
%C2 = sym( zeros(3,3,3,3) );
% C(1,1,1,1) = C_1111;

for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                C(ii,jj,kk,ll) = [ 'C_',num2str(ii),num2str(jj),num2str(kk),num2str(ll) ];
%                 C2(ii,jj,kk,ll) = [ 'C_',num2str(ii),num2str(jj),num2str(kk),num2str(ll) ];
            end
        end
    end
end

for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                 if ii == jj 
                   if ll > kk % ij23
                      C(ii,jj,ll,kk) = [ 'C_',num2str(ii),num2str(jj),num2str(kk),num2str(ll) ];  
                   end
                 elseif ii < jj
                   if ll > kk
                      C(ii,jj,ll,kk) = [ 'C_',num2str(ii),num2str(jj),num2str(kk),num2str(ll) ];
                   elseif ll == kk
                      C(ii,jj,kk,ll) = C(kk,ll,ii,jj);  
                   end 
                 elseif ii > jj
                   C(ii,jj,kk,ll) = C(jj,ii,kk,ll);
                 end
            end
        end
    end
end
% 
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                if ( ii + jj ) >= ( kk + ll )
                   C(ii,jj,kk,ll) = C(kk,ll,ii,jj);
%                    if ii == jj
%                       C(ii,jj,kk,ll) = [ 'C_',num2str(kk),num2str(ll),num2str(ii),num2str(jj) ]; 
%                    end
                end   
            end
        end
    end
end

C11 = [ C(1,1,1,1) C(1,1,1,2) C(1,1,1,3);
        C(1,1,2,1) C(1,1,2,2) C(1,1,2,3); 
        C(1,1,3,1) C(1,1,3,2) C(1,1,3,3) ];
    
C12 = [ C(1,2,1,1) C(1,2,1,2) C(1,2,1,3);
        C(1,2,2,1) C(1,2,2,2) C(1,2,2,3); 
        C(1,2,3,1) C(1,2,3,2) C(1,2,3,3) ];
    
C13 = [ C(1,3,1,1) C(1,3,1,2) C(1,3,1,3);
        C(1,3,2,1) C(1,3,2,2) C(1,3,2,3); 
        C(1,3,3,1) C(1,3,3,2) C(1,3,3,3) ];

C21 = [ C(2,1,1,1) C(2,1,1,2) C(2,1,1,3);
        C(2,1,2,1) C(2,1,2,2) C(2,1,2,3); 
        C(2,1,3,1) C(2,1,3,2) C(2,1,3,3) ];
    
C22 = [ C(2,2,1,1) C(2,2,1,2) C(2,2,1,3);
        C(2,2,2,1) C(2,2,2,2) C(2,2,2,3); 
        C(2,2,3,1) C(2,2,3,2) C(2,2,3,3) ];
    
C23 = [ C(2,3,1,1) C(2,3,1,2) C(2,3,1,3);
        C(2,3,2,1) C(2,3,2,2) C(2,3,2,3); 
        C(2,3,3,1) C(2,3,3,2) C(2,3,3,3) ];   
    
C31 = [ C(3,1,1,1) C(3,1,1,2) C(3,1,1,3);
        C(3,1,2,1) C(3,1,2,2) C(3,1,2,3); 
        C(3,1,3,1) C(3,1,3,2) C(3,1,3,3) ];
    
C32 = [ C(3,2,1,1) C(3,2,1,2) C(3,2,1,3);
        C(3,2,2,1) C(3,2,2,2) C(3,2,2,3); 
        C(3,2,3,1) C(3,2,3,2) C(3,2,3,3) ];
    
C33 = [ C(3,3,1,1) C(3,3,1,2) C(3,3,1,3);
        C(3,3,2,1) C(3,3,2,2) C(3,3,2,3); 
        C(3,3,3,1) C(3,3,3,2) C(3,3,3,3) ];    
    
% write stiffness tensor as 9 by 9 hypermatrix    
Chyper = [ C11 C12 C13;
           C21 C22 C23;
           C31 C32 C33 ];   

%%
% definition of Kelvin bases
K1 = sym(I(:,1)*I(:,1)'); % e1 dyad e1
K2 = sym(I(:,2)*I(:,2)'); % e2 dyad e2
K3 = sym(I(:,3)*I(:,3)'); % e3 dyad e3
K4 = sym(k*( I(:,2)*I(:,3)' + I(:,3)*I(:,2)' )); % e2 dyad e3
K5 = sym(k*( I(:,1)*I(:,3)' + I(:,3)*I(:,1)' )); % e1 dyad e3
K6 = sym(k*( I(:,1)*I(:,2)' + I(:,2)*I(:,1)' )); % e1 dyad e2

KK = zeros(3,3,6);

for bb = 1:6
    for ii = 1:3
        for jj = 1:3
            if bb == 1
               K = K1;
            elseif bb == 2
               K = K2; 
            elseif bb == 3
               K = K3; 
            elseif bb == 4
               K = K4; 
            elseif bb == 5
               K = K5;    
            elseif bb == 6
               K = K6;   
            end
            KK(ii,jj,bb) = K(ii,jj);
        end
    end
end

%%
% definition of Voigt bases
V1 = sym(I(:,1)*I(:,1)'); % e1 dyad e1
V2 = sym(I(:,2)*I(:,2)'); % e2 dyad e2
V3 = sym(I(:,3)*I(:,3)'); % e3 dyad e3
V4 = sym(v*( I(:,2)*I(:,3)' + I(:,3)*I(:,2)' )); % e2 dyad e3
V5 = sym(v*( I(:,1)*I(:,3)' + I(:,3)*I(:,1)' )); % e1 dyad e3
V6 = sym(v*( I(:,1)*I(:,2)' + I(:,2)*I(:,1)' )); % e1 dyad e2

VV = zeros(3,3,6);

for aa = 1:6
    for ii = 1:3
        for jj = 1:3
            if aa == 1
               V = V1;
            elseif aa == 2
               V = V2; 
            elseif aa == 3
               V = V3; 
            elseif aa == 4
               V = V4; 
            elseif aa == 5
               V = V5;    
            elseif aa == 6
               V = V6;   
            end
            VV(ii,jj,aa) = V(ii,jj);
        end
    end
end


%%
% 6 by 6 Matrix representation of C with Kelvin bases
CK = sym( zeros(6,6));

for aa = 1:6
    for bb = 1:6
        for ii = 1:3
            for jj = 1:3
                for kk = 1:3
                    for ll = 1:3
                        CK(aa,bb) = CK(aa,bb) + KK(ii,jj,aa)*C(ii,jj,kk,ll)*KK(kk,ll,bb); 
                    end
                end
            end
        end
    end
end

% 6 by 6 Matrix representation of C with Voigt bases
CV = sym( zeros(6,6));

for aa = 1:6
    for bb = 1:6
        for ii = 1:3
            for jj = 1:3
                for kk = 1:3
                    for ll = 1:3
                        CV(aa,bb) = CV(aa,bb) + VV(ii,jj,aa)*C(ii,jj,kk,ll)*VV(kk,ll,bb); 
                    end
                end
            end
        end
    end
end
%% definition of symmetric second-order tensors, i.e. strain E and stress T
%
E = sym( zeros(3,3) );

for ii = 1:3
    for jj = 1:3
         %if ( ii ) >= ( jj )
                E(ii,jj) = [ 'E_',num2str(ii),num2str(jj)];
         %end
    end
end

for ii = 1:3
    for jj = 1:3
         if ( ii ) >= ( jj )
                E(ii,jj) = E(jj,ii);
         end
    end
end

T = sym( zeros(3,3) );

for ii = 1:3
    for jj = 1:3
         %if ( ii ) >= ( jj )
                T(ii,jj) = [ 'T_',num2str(ii),num2str(jj)];
         %end
    end
end

for ii = 1:3
    for jj = 1:3
         if ( ii ) >= ( jj )
                T(ii,jj) = T(jj,ii);
         end
    end
end

%%
% 6 by 1 vector representation of E and T with Voigt bases
EV = sym( zeros(6,1));

for aa = 1:6
        for ii = 1:3
            for jj = 1:3
                        EV(aa) = EV(aa) + E(ii,jj)*VV(ii,jj,aa); 
            end   
        end
end
for aa = 4:6
       EV(aa) = EV(aa) * 2; 
end

TV = sym( zeros(6,1));

for aa = 1:6
        for ii = 1:3
            for jj = 1:3
                        TV(aa) = TV(aa) + T(ii,jj)*VV(ii,jj,aa); 
            end   
        end
end



% 6 by 1 vector representation of E with Kelvin bases
EK = sym( zeros(6,1));

for aa = 1:6
        for ii = 1:3
            for jj = 1:3
                        EK(aa) = EK(aa) + E(ii,jj)*KK(ii,jj,aa); 
            end   
        end
end

TK = sym( zeros(6,1));

for aa = 1:6
        for ii = 1:3
            for jj = 1:3
                        TK(aa) = TK(aa) + T(ii,jj)*KK(ii,jj,aa); 
            end   
        end
end

% %% results
% 
% % matrix representations
% 
% Chyper
 CV
 CK
% 
% % compute inverse of matrix
% iCV = inv(CV);
% iCK = inv(CK);
% 
% % compute eigenvalues and eigenvectors
% %[V,D] = eig(CV) % slow!
% [W, lambda] = eig(CK) % slow!
% %eig(C)
% %eig(CV)
% %eig(CK,'matrix')
% 
% %% some computations
% % scalar valued distance d of stiffnesses
% D  = (CK - CV);
% d = norm(D,'fro')./norm(CK,'fro');
% simplify(d,'steps',3)
% %simplify(D,'steps',3)
% 
% IV = simplify(CV * iCV);
% IK = simplify(CK * iCK);
% 
% % Frobenius norm of elasticity matrices
% nCV = simplify(norm(CV,'fro'));
% nCK = simplify(norm(CK,'fro'));
% 
%%
EK
EV

TK
TV

