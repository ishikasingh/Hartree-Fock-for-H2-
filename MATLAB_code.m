%Ishika Singh and Shivali Agrawal
a1=[0.168856 0.623913 3.42525]; %exponents of the STO-3G basis functions 
a=[a1 a1]; % for 2 1s orbitals
%guessed coefficients, we won't require to guess beacuse we don't need
%density matrix which is required for 2-electron integrals

%initiating other matrices
T=zeros(6,6);
S=zeros(6,6);
F=zeros(6,6);
X=zeros(6,6);
V=zeros(6,6);

% function required for evaluating Vext matrix elements
F0= @(t) 0.5*(pi/t)^0.5*erf(t^0.5);

% starting value of nuclear distance, since 0.5 a.u. is the minimum atomic
% radius;
% the nuclear distance is varied till 5 a.u., the energy is constant beyond
% this point as there isn't any interaction between the atoms
ii=1;
for R=0.5:0.05:5
R
%evaluating the integrals for different blocks of the Fock matrix

%top-left block: 1s_a:1s_a type integrals
for i=1:3
    for j=1:3
    % all the integrals are reduced analytically to less complicated forms
    % so that the computations expenses could be cut down
    %(ref: ?Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory?, by Attila
    % Szabo and Neil S. Ostlund, Dover Publications, ISBN-10: 0486691861.)
    % p, m, K are some short hand notations for certain constant terms used
    % in the integrals 
    p=a(i)+a(j);
    m=a(i)*a(j);
    K=1;
    %overlap term
    S(i,j)=K*(pi/p)^1.5;
    %kinetic energy term
    %%T(i,j)=c(i)*c(j)*(m/p)*3*(pi/p)^1.5;
    T(i,j)=0.5*(m/p)*(6-4*(m/p)*(0)^2)*S(i,j);
    k=F0(p*R^2);
    %nuclear attraction energy term
    V(i,j)=(-2*pi/p)*(1+k);
    %fock matrix
    F(i,j)=T(i,j)+V(i,j);
    end
end
%similar computations are performed for other blocks as well
%top-right block: 1s_a:1s_b type integrals
for i=1:3
    for j=4:6
    p=a(i)+a(j);
    m=a(i)*a(j);
    K=exp(-(m*R^2)/p);
    S(i,j)=K*(pi/p)^1.5;  
    %%T(i,j)=c(i)*c(j)*(m/p)*(3-2*m*R^2/p)*(pi/p)^1.5*exp(-m*R^2/p);
    T(i,j)=0.5*(m/p)*(6-4*(m/p)*(R)^2)*S(i,j);
    Rp=a(j)*R/(p);
    V(i,j)=(-2*pi/p)*exp(-(m/p)*R^2)*(F0(p*Rp^2)+F0(p*(Rp-R)^2));
    F(i,j)=T(i,j)+V(i,j);
    end
end
%bottom-left block: 1s_b:1s_a type integrals
for i=4:6
    for j=1:3
    p=a(i)+a(j);
    m=a(i)*a(j);
    K=exp(-(m*R^2)/p);
    S(i,j)=K*(pi/p)^1.5;   
    %%T(i,j)=c(i)*c(j)*(m/p)*(3-2*m*R^2/p)*(pi/p)^1.5*exp(-m*R^2/p);
    T(i,j)=0.5*(m/p)*(6-4*(m/p)*(R)^2)*S(i,j);
    Rp=a(i)*R/(p);
    V(i,j)=(-2*pi/p)*exp(-(m/p)*R^2)*(F0(p*Rp^2)+F0(p*(Rp-R)^2));
    F(i,j)=T(i,j)+V(i,j);
    end
end
%bottom-right block: 1s_b:1s_b type integrals
for i=4:6
    for j=4:6
    p=a(i)+a(j);
    m=a(i)*a(j);
    K=1;
    S(i,j)=K*(pi/p)^1.5; 
    %%T(i,j)=c(i)*c(j)*(m/p)*3*(pi/p)^1.5;
    T(i,j)=0.5*(m/p)*(6-4*(m/p)*(0)^2)*S(i,j);
    Rp=R;
    V(i,j)=(-2*pi/p)*(F0(p*Rp^2)+1);
    F(i,j)=T(i,j)+V(i,j);
    end
end
% we do not require the SCF loop since we don't have 2 electron integrals,
% hence the entire Fock matrix is independent of coefficients
% reformulating the eigenvalue equation

[U,s]=eig(S); %U is the unitary matrix that diagonalises S
for j=1:6
    X(:,j)=U(:,j)/s(j,j)^0.5; %X=Us^-0.5 is the orthogonalising transformation matrix
end
F_mod = X'*F*X; %modified Fock Matrix
[c_mod,e]=eig(F_mod);
C=X*c_mod;
eigval= diag(e); 
[e_min,index]=min(eigval);%taking the minimum eigenvalue since we are computing the ground state energy
c=C(:,index); %eigenvector corresponding to the ground state
c=c'; %updating the coefficients 
E(ii)=c*F*c'+1/R; %total energy
plot(R*0.529,E(ii),'*')
ylabel('Energy (Hartrees)')
xlabel('Nuclear separation (A)')
title('Ground state energy curve') 
hold on
ii=ii+1;
end
R=0.5:0.05:5;
[E_min,index]=min(E);
Minimum_energy=E_min
Bond_length=R(index)*0.529

r = -4:0.01:4;  
Y1 =  (c(1).*exp(-a(1).*r.^2) + c(2).*exp(-a(2).*r.^2) + (c(3)).*exp(-a(3).*r.^2));
Y2 = (c(4).*exp(-a(4).*(r-1).^2) + c(5).*exp(-a(5).*(r-1).^2) + (c(6)).*exp(-a(6).*(r-1).^2));
figure()
plot(r,Y1,'LineWidth',2)
ylabel('Wave functions')
xlabel('r (A)')
title('For the two centres separated by 1A')
hold on
plot(r,Y2,'LineWidth',2)
legend('Centred at r=0','Centred at R=1')
hold off 

figure()
Y33 = 4*pi.*(r.^2.*(Y1).^2 + (r-1).^2.*Y2.^2);
plot(r,Y33,'LineWidth',2)
ylabel('Radial distribution functions')
xlabel('r (A)')
title('For the two centres separated by 1A')
