%----------------------------------------------------------------------
%                          Data-Driven SSI 
%----------------------------------------------------------------------
%   
%   The algorithm 'os_ssid.m' identifies a reference based stochastic
%   state space model from output only data.
%
%  [sys,modal] = os_ssid(y,i,fs,Ir,W,n,AUXin,sil)
%
%              sys = [A,C,K,Ro,AUX]
%              modal = [ss,EW,fn,d,EV]
% 
%   Inputs:
%            y:   matrix of measured outputs
%            i:   block size - number of block rows in Hankel matrices
%           fs:   sample frequency
%           Ir:   projection channels index
%            W:   optional weighting function
%                      PC:    Principle component algorithm
%                      UPC:   Unweighted principal component algorithm [default]
%                      CVA:   Canonical variate analysis algorithm
%            n:   model order - optional subspace selection
%        AUXin:   optional auxilary variable to increase speed (default [])
%          sil:   when equal to 1 no text output is generated
%           
%   Outputs:
%           A,C,K,Ro: combined state space system
%           
%                  x_{k+1) = A x_k + K w_k        
%                    y_k   = C x_k + v_k
%                 cov(e_k) = Ro
%
%  
%             A:  state matrix
%             C:  output matrix
%             K:  Kalman filter gain matrix
%            Ro:  system error covariance
%           AUX:  optional auxilary variable to increase speed
%            ss:  singular values from the SVD of the projection matrix
%            EW:  complex conjugate eigenvalues (Eigenwert)
%            fn:  eigenfrequency (Hz)
%             d:  damping ratio (%)
%            EV:  mode shape (rotated but unscaled) (Eigenvektor)
%
% This algorithm was developed from the following literature:
% 1. VanOverschee, P. and De Moor, B. (1996). Subspace Identification for
% Linear Systems: Theory, Implementation, Applications. Kluwer Academic
% Publishers.
% 2. Peeters, B. and de Roeck, G. (1999). Reference-based
% stochastic-subspace identification for output-only modal analysis.
% Mechanical Systems and Signal Processing.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                         Keith Soal 2017
%                 questions to keithsoal at gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [sys,modal] = os_ssid(y,i,fs,Ir,W,n,AUXin,sil)

% Check the arguments
if (nargin < 8);sil = 0;end
if (nargin < 7);AUXin = [];end
if (nargin < 6);n = [];end
if (nargin < 5);W = 'UPC';end
if (nargin < 4);Ir = [];end

if sil == 0
disp(' ');
disp('   --------------------------------------------------------');
disp('       Data-Driven Stochastic Subspace Identification');
disp('   --------------------------------------------------------');
end

% Turn the data into row vectors
[l,ny] = size(y);if (ny < l);y = y';[l,ny] = size(y);end
if ((ny-2*i+1) < (2*l*i));error('Not enough data points');end

% determine number of projection channels r
if (~isempty(Ir)); r = length(Ir);else r = l; end

% Weighting algorithm
if strcmp(W,'PC');Wn = 1;
elseif strcmp(W,'UPC');Wn = 2;
else Wn = 3;
end
Waux = 1;
W = Wn;

% Number of columns in the Hankel matrix
j = ny-2*i+1;

% Check compatibility of AUXin
Uaux = [];

%--------------------------------------------------------------------------
%                            Hankel matrix
%--------------------------------------------------------------------------
if (isempty(AUXin));
    y = y/sqrt(j);
    ih = 2*i;
    % Make a block-row vector out of y
    if r == l
        % Block-Hankel matrix r=l  (no reference channels)
        H=zeros(l*ih,j);
        for k=1:ih
            H((k-1)*l+1:k*l,:)=y(:,k:k+j-1);
        end
    elseif r < l
        % Block-Hankel matrix r<l  (reference channels)
        Hp=zeros(r*i,j); % past part
        for k=1:i
            Hp((k-1)*r+1:k*r,:)=y(Ir,k:k+j-1);
        end
        Hf=zeros(l*i,j); % future part
        for k=1:i
            Hf((k-1)*l+1:k*l,:)=y(:,k+i:k+i+j-1);
        end
        H = [Hp;Hf];
        clearvars Hp Hf;
    else
        disp('r cannot exceed l')
        doc return
    end
    
%--------------------------------------------------------------------------
%                            QR Factorisation
%--------------------------------------------------------------------------
    if sil == 0
        disp('      Computing ... R from QR factorisation');
    end
    
    R = triu(qr(H'))'; % R factor
    Rorg = R(1:2*l*i,1:2*l*i);
    
    if r < l
        % Extract R spaces
        R21 = R(r*i+1:r*i+r,1:r*i);
        R22 = R(r*i+1:r*i+r,r*i+1:r*i+r);
        R31 = R(r*i+r+1:r*i+r+(l-r),1:r*i);
        R32 = R(r*i+r+1:r*i+r+(l-r),r*i+1:r*i+r);
        R33 = R(r*i+r+1:r*i+r+(l-r),r*i+1+r:r*i+r+(l-r));
        R43 = R(r*i+r+1+(l-r):r*i+r+(l-r)+l*(i-1),r*i+1+r:r*i+r+(l-r));
        R41 = R(r*i+r+1+(l-r):r*i+r+(l-r)+l*(i-1),1:r*i);
        R42 = R(r*i+r+1+(l-r):r*i+r+(l-r)+l*(i-1),r*i+1:r*i+r);
    else
        % Truncate R for original method
%         R = R(1:2*i*l,1:2*i*l);
    end
else
    R = AUXin(2:2*i*l+1,1:2*l*i);
    bb = 2*i*(l)+1;
end

if (isempty(AUXin)) && r == l
    Rf = R(l*i+1:2*l*i,:); % Future outputs
    Rp = R(1:l*i,:); % Past outputs
end

%--------------------------------------------------------------------------
%                            Oblique Projection 
%--------------------------------------------------------------------------
if (isempty(AUXin)) && r == l
    % no reference channels
    P = [Rf(:,1:l*i),zeros(l*i,l*i)];
elseif (isempty(AUXin)) && r < l
    % reference channels
    P = [R21;R31;R41];
else
    P = AUXin(bb+1:bb+l*i,1:2*(l)*i);
    bb = bb+l*i;
end

%--------------------------------------------------------------------------
%                           Weighting methods
%--------------------------------------------------------------------------
if (isempty(AUXin))
    % PC algorithm :
    if W == 1
        R11 = R(1:r*i,1:r*i);
        Pw = P*R11';
        % UPC algorithm
    elseif W == 2
        Pw = P;
        % CVA algorithm
    else
        R21 = R(r*i+1:2*i*r,1:r*i);
        R22 = R(r*i+1:2*i*r,r*i+1:2*i*r);
%         W1 = inv(sqrtm(R21*R21'+R22*R22'));
%         Pw = W1*P;
        wc = sqrtm(R21*R21'+R22*R22');
        Pw = wc\P;
    end
    
%--------------------------------------------------------------------------
%                    Singular Value Decomposition (SVD)
%--------------------------------------------------------------------------
    if sil == 0
        disp('      Computing ... SVD');
    end
    
    [U,S,V] = svd(Pw);
    ss = diag(S);
    
%     if W == 3;U = W1*U;end % CVA
%     clearvars V S WOW
else
    U = AUXin(bb+1:bb+l*i,1:l*i);
    ss = AUXin(bb+1:bb+l*i,l*i+1);
end


%--------------------------------------------------------------------------
%                       Subspace Selection (n)
%                           (System Order)
%--------------------------------------------------------------------------
if (isempty(n))
    figure()
    bar(1:length(ss),ss);
    title('Singular Values');
    xlabel('Order');
    n = 0;
    while (n < 1) || (n > l*i-1)
        n = input('Input system subspace:');
    end
    U1 = U(:,1:n);
elseif (isempty(n))== 0
    U1 = U(:,1:n);
end

%--------------------------------------------------------------------------
%                   Determine Oi and Ois
%--------------------------------------------------------------------------
Oi  = U1*diag(sqrt(ss(1:n))); % observability matrix
Ois = Oi(1:l*(i-1),:); % shifted observability matrix

%--------------------------------------------------------------------------
%                       Determining A and C
%--------------------------------------------------------------------------
if sil == 0
   disp(['Computing ... System matrices A,C (Order ',num2str(n),')']);
end


if r == l
Rhs = [pinv(Oi)*R(l*i+1:2*l*i,1:l*i),zeros(n,l)];
Lhs = [pinv(Ois)*R(l*i+l+1:2*l*i,1:l*i+l) ; ...
    R(l*i+1:l*i+l,1:l*i+l)];
sol = Lhs/Rhs; % Solve least square
else
Xh = pinv(Oi)*P; % state sequence
% Ps = [R41,R42,zeros(length(R41),size(R33,2))];
Ps = [R41,R42,zeros(size(R41,1),size(R33,2))];
Xhs = pinv(Ois)*Ps;  % shifted state sequence
Yii = [R21 R22 zeros(size(R22,1),size(R33,2));R31 R32 R33];
Rhs = [Xh,zeros(n,l)];
Lhs = [Xhs;Yii];
% Solve least square
sol = Lhs/Rhs;
end

% Extract the system matrices A and C
A = sol(1:n,1:n); % state matrix
C = sol(n+1:n+l,1:n); % output matrix
res = Lhs - sol*Rhs; % Residuals
sys.A = A;
sys.C = C;

%--------------------------------------------------------------------------
%                    Noise covariances (Q,S and R)
%--------------------------------------------------------------------------
if (norm(res) > 1e-10)
    % Determine QSR from the residuals
    if sil == 0
        disp(['Computing ... System matrices G,L0 (Order ',num2str(n),')']);
    end
    % Determine the residuals
    cov = res*res'; 			% Covariance
    Qs = cov(1:n,1:n);Ss = cov(1:n,n+1:n+l);Rs = cov(n+1:n+l,n+1:n+l);
    
    sig = dlyap(A,Qs);
    G = A*sig*C' + Ss;
    L0 = C*sig*C' + Rs;
    
%--------------------------------------------------------------------------
%                 Kalman gain K and system noise covariance Ro
%--------------------------------------------------------------------------
    if sil == 0
        disp('Computing ... Riccati solution')
    end
    [K,Ro] = gl2kr(A,G,C,L0);
else
    Ro = [];
    K = [];
end
sys.K = K;
sys.Ro = Ro;


%--------------------------------------------------------------------------
%                Auxilary matrix for itereative estimates    
%--------------------------------------------------------------------------
% Make AUX when needed
% if nargout > 4
if isempty(AUXin)
  AUX = zeros((4*l)*i+1,2*(l)*i);
%   Uaux = [];  % i added this...
%   if Uaux == [];Uaux = 0;end
  if isempty(Uaux);Uaux = 0;end
  ds_flag = 2;
  infon = [ds_flag,i,Uaux,y(1,1),Waux]; % in/out - i - u(1,1) - y(1,1) - W
  AUX(1,1:5) = infon;
  bb = 1;
  AUX(bb+1:bb+2*(l)*i,1:2*(l)*i) = Rorg;
  bb = bb+2*(l)*i;
  AUX(bb+1:bb+l*i,1:2*(l)*i) = P;
  bb = bb+l*i;
  AUX(bb+1:bb+l*i,1:l*i) = U;
  AUX(bb+1:bb+l*i,l*i+1) = ss;
  sys.AUX = AUX;
else
   sys.AUX = AUXin; 
end
% end

%--------------------------------------------------------------------------
%                         Modal analysis
%--------------------------------------------------------------------------
% eigenvalue problem
[Vs,Ds] = eig(A);
% Sort eigenvalues and eigenvectors in ascending order
[Dum,I] = sort(abs((log(diag(Ds)))/(1/fs)));
Dsd = diag(Ds);
Dsd = Dsd(I);
DsD = log(Dsd)/(1/fs);
Vs = Vs(:,I);
% Convert eigenvector to have physical meaning
Vs = C*Vs;
[m,n]=size(Vs);
Vs=Vs(:,1:n);
% Rotate Eigenvectors
phi=angle(Vs(1,:));
ij = 0.0000 + 1.0000i;
phi=diag(exp(-ij*phi));
% Vs=flip(Vs*phi);
% extract modal parameters from state space identification
ws = abs(log(Dsd)/(1/fs));
d = (-real(log(Dsd)/(1/fs))./ws)*100;
fn = ws/(2*pi);

modal.ss = ss;
modal.EW = DsD;
modal.fn = fn;
modal.d = d;
modal.EV = Vs;

%--------------------------------------------------------------------------
%                         Command window display
%--------------------------------------------------------------------------
if sil == 0
    no = 1:length(fn);
    no = no(:);
    disp('SSI Data')
    fprintf('No.       freq.         damp.     \n')
    for ii = 1:length(fn)
        fprintf('%4.0f %11.4f      %10.4f \n', ...
            no(ii),fn(ii),d(ii)*100)
    end
end

end




