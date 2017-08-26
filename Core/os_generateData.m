% Generates simulated data
% openSID
% Keith Soal   01-08-2017

% devnotes
% 1. generalise for higher DOFs

function [out,Y] = os_generateData(m,k,z,fv,fs,ndof,L,sil)
%
%    out:    time data
%      Y:    spectrum
%     EW:    complex eigenvalues
%     EV:    eigen vector
%     fn:    natural frequency (Hz)
%     cn:    damping ratio (%)
%      E:    excitation force
%
%      m:    mass (kg)
%      k:    stiffness (N/m)
%      z:    damping (%)
%     fv:    input force
%              1. random in amplitude and phase
%              2. burst random in amplitude no phase
%              3. multiple sine sweep
%              4. dirac delta
%     fs:    sample frequency (Hz)
%   ndof:    number of dofs
%    sil:    display if 0 or []

% Check the arguments
if (nargin < 8);sil = 1;end
if (nargin < 7);L = 1E4;end % length of vector
if (nargin < 6);ndof = 3;end
if (nargin < 5);fs = 40;end % sample frequency
if (nargin < 4);fv = 1;end


%-----------------------------------------
%         system properties
%-----------------------------------------
% m = [100 100 100]; % mass
% k = [100000 100000 100000]; % stiffness
% z = [0.1 0.12 0.14]./100; % damping ratio
[M,C,K] = sysmtx(m,z,k);

%-----------------------------------------
%        generate time data
%-----------------------------------------
% E = exp(2i*pi*rand(L,1)); % random excitation
if fv == 1
% 1.random in amplitude and phase
E = exp(2i*pi*rand(L,1)); % random excitation
elseif fv == 2
% 2.burst random in amplitude no phase
E = 2*rand(L,1); E = E - mean(E);
E(1:200) = zeros(1,200);E(end-199:end) = zeros(1,200);
elseif fv == 3
% 3.multiple sine sweep
t= 0:0.1:10;
fo=0.8;f1=20;
E = chirp(t,fo,10,f1,'logarithmic');
E = E(1:end-1)';
E = repmat(E,100,1);
elseif fv == 4
% 4.dirac delta
E = 1;
end
% sine sweep
% direc delta
[out,Y] = calcTimeData(L,fs,M,C,K,E);

%-----------------------------------------
%        analytical solution
%-----------------------------------------
[EV,EW] = ansol(M,C,K);
fn = abs(EW)/(2*pi);
cn = -real(EW)./abs(EW).*100;

if sil == 0
    figure()
    f = linspace(0,fs/2,L);
    semilogy(f,abs(Y))
    title('Frequency Response','Interpreter','latex')
    xlabel('Frequency (Hz)','Interpreter','latex')
    ylabel('Amplitude','Interpreter','latex')

    figure()
    plot(out)
    title('Time Signal','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
    ylabel('Amplitude','Interpreter','latex')
end
end

%-----------------------------------------
%         static functions
%-----------------------------------------

function [M,C,K] = sysmtx(m,z,k)
%-----------------------------------------
%         system matrices
%-----------------------------------------
    M = [m(1) 0 0;0 m(2) 0;0 0 m(3)]; % mass matrix
%     M = m*eye(ndof); % mass matrix
    K = [(k(1)+k(2)) -k(2) 0;-k(2) (k(2)+k(3)) -k(3);0 -k(3) k(3)]; % stiffness matrix
%     K = [(k+k) -k 0;-k (k+k) -k;0 -k k]; % stiffness matrix
    % damping matrix
    [~, ev] = eig(K,M);
    ev_mk = sqrt(ev);
    for i = 1:3
        c(i) = (2*z(i)*(ev_mk(i,i)))*M(i,i);
    end;
    C = [(c(1)+c(2)) -c(2) 0;-c(2) (c(2)+c(3)) -c(3);0 -c(3) c(3)];
end

function [out,Y] = calcTimeData(L,fs,M,C,K,E)
% simulate time data
    DoF = size(M,1);
    s = 2i*pi*(0:L-1)'/L*fs/2;
    H = zeros(DoF,DoF,numel(s));
    for n = 1:numel(s)
        Z = s(n).^2*M+s(n)*C+K;
        H(:,:,n) = Z\eye(DoF); % frf
    end
    
    idx = [1 2 3];
    H = reshape(H,[],size(H,3)).';
    for n = 1:3
        Y(:,n) = H(:,n).*E;
    end
    
    out = real(ifft([Y ; conj(Y(end:-1:2,:))])); % Output
end

function [EW,EV] = ansol(M,C,K)
%--------------------------------------------------------------------------
%                     Analytical State Space Solution
%--------------------------------------------------------------------------
A=[C M; M 0*M];
B=[K 0*M; 0*M -M];
[V,D]=eig(B,-A);
% Sort in descending order
[Dum,I]=sort(diag(abs(imag(D))));
p=diag(D);
p=p(I);
V=V(:,I);
% Rotate vectors to real first element (row 1)
phi=angle(V(1,:));
phi=diag(exp(-j*phi));
V=V*phi;
% Scale to unity Modal A
Ma=V.'*A*V;
for col = 1:length(V(1,:))
    V(:,col)=V(:,col)/sqrt(Ma(col,col));
end
[m,n]=size(V);
EW=p(1:2:m);
EV=V(1:m/2,1:2:n);
end

