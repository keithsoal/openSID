%--------------------------------------------------------------------------
%                        Covariance-Driven SSI 
%--------------------------------------------------------------------------
%   
%   The algorithm 'os_ssic' identifies a reference based stochastic
%   state space model from output only covariance matrices.
%
%    [sys,modal] = os_ssic(y,i,fs,Ir,n,AUXin,sil)
%  [sys,modal] = os_ssid(y,i,fs,Ir,W,n,AUXin,sil)
%
%              sys = [A,Oi,Ci,R,AUX]
%              modal = [ss,EW,fn,d,EV]
% 
%   Inputs:
%           y:    matrix of measured outputs
%           i:    block size - number of block rows in Toeplitz matrix
%          fs:    sample frequency
%          Ir:    projection channels index
%           n:    model order - optional subspace selection
%       AUXin:    optional auxilary variable to increase speed (default [])
%         sil:    when equal to 1 no text output is generated
%           
%   Outputs:
%           A,Oi,Ci: combined state space system
%           
%                  x_{k+1) = A x_k + w_k       
%                    y_k   = Ci x_k + v_k
%
%             A:  state matrix
%            Oi:  observability matrix
%            Ci:  output matrix
%           AUX: optional auxilary variable to increase speed
%           ss:  singular values from the SVD of the projection matrix
%           Dsc: complex conjugate eigenvalues
%           fn:  eigenfrequency (Hz)
%           d:   damping ratio (%)
%           Vs:  mode shape (rotated but unscaled)
%
% This algorithm was developed from the following literature:
% 1. Peeters, B. and de Roeck, G. (1999). Reference-based
% stochastic-subspace identification for output-only modal analysis.
% Mechanical Systems and Signal Processing.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                         Keith Soal 2017
%                 questions to keithsoal at gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [sys,modal] = os_ssic(y,i,fs,Ir,n,AUXin,lead,sil)

% Check the arguments
if (nargin < 6);sil = 0;end
if (nargin < 6);AUXin = [];end
if (nargin < 5);n = [];end
if (nargin < 4);Ir = [];end

if sil == 0
disp(' ');
disp('   -------------------------------------------------------------');
disp('       Covariance-Driven Stochastic Subspace Identification');
disp('   -------------------------------------------------------------');
end

% identify selected reference channels
if (~isempty(Ir));
    ir = length(Ir);
    l = size(y,2);
    Y = y(:,Ir); % Select reference outputs
else
    ir = size(y,2);
    l = ir;
    Y = y;
end

%--------------------------------------------------------------------------
%                         Correlation matrix
%--------------------------------------------------------------------------
if (isempty(AUXin));
    N = length(y);
    for nn = 1:ir
        for m = 1:l
            [R(:,m,nn),t] = xcorr(y(:,m),Y(:,nn),'unbiased');
        end
    end
    % correlation from 0 lag with time shift (lead)
    [~,lag] = find(t==0);
    Rs = zeros(round(size(R,1)/2)-lead,3,3);
    for nn = 1:ir
        for m = 1:l
            Rs(:,m,nn) = R(lag+lead:end,m,nn); % start at shifted time (Anders Brandt suggestion IOMAC)
        end
    end    
    clear R;
    R = Rs;
else
    R = AUXin;
end
AUX = R;

%--------------------------------------------------------------------------
%                         Toeplitz matrix
%--------------------------------------------------------------------------
% Build Hankel matrix
[~,l,r]=size(R);
Hank=zeros(l*i,r*i);
for c = 1:i
    for w = 1:i
        Hank(1+(w-1)*l:w*l,1+(c-1)*r:c*r)=squeeze(R(w+c-1,:,:))';
    end
end
% Flip Hankel matrix to Toeplitz matrix
Toep = flip(Hank,2);

%--------------------------------------------------------------------------
%              Singular value decomposition (SVD)
%--------------------------------------------------------------------------
if sil == 0
       disp('Computing ... SVD');
end
[U,S,V]=svd(Toep);
ss = diag(S);
% Uc = transpose(U)*U;
% Vc = transpose(V)*V;

%--------------------------------------------------------------------------
%                       Subspace Selection (n)
%                           (System Order)
%--------------------------------------------------------------------------
if (isempty(n))
    figure()
    bar(1:l*i,ss);
    title('Singular Values original method');
    xlabel('Order');
    n = 0;
        while (n < 1) || (n > l*i-1)
            n = input('Input system subspace:');
        end
    U1 = U(:,1:n);
    V1 = V(:,1:n);
else
    U1 = U(:,1:n);
    V1 = V(:,1:n);
end

if sil == 0
   disp(['Computing ... System matrices A,O,C_ref (Order ',num2str(n),')']);
end

%--------------------------------------------------------------------------
%           Observability (Oi) and controllability (Ci)
%--------------------------------------------------------------------------
Oi  = U1*diag(sqrt(ss(1:n)));
C = Oi(1:l,:);
Ci = diag(sqrt(ss(1:n)))*transpose(V1);
% The pseudo inverses
Oip = pinv(Oi);
Cip = pinv(Ci);

%--------------------------------------------------------------------------
%                     Shifted Toepltiz (Toeps)
%--------------------------------------------------------------------------
Hanks=zeros(l*i,r*i);
for c = 1:i
    for w = 1:i
        Hanks(1+(w-1)*l:w*l,1+(c-1)*r:c*r)=squeeze(R(w+c,:,:))';
    end
end
Toeps = flip(Hanks,2);

%--------------------------------------------------------------------------
%                     State matricex (A)
%--------------------------------------------------------------------------
A = Oip*Toeps*Cip;

sys.A = A;
sys.Oi = Oi;
sys.Ci = Ci;
sys.R = R;
sys.AUX = AUX;

%--------------------------------------------------------------------------
%                          Modal analysis
%--------------------------------------------------------------------------
% eigenvalue problem
[Vs,Ds] = eig(A);
% Sort eigenvalues and eigenvectors in ascending order
[Dum,I]=sort(abs((log(diag(Ds)))/(1/fs)));
Dsd=diag(Ds);
Dsd=Dsd(I);
DsD = log(Dsd)/(1/fs);
Vs=Vs(:,I);
% Convert eigenvector to have physical meaning
Vs = C*Vs;
[m,n]=size(Vs);
Vs=Vs(:,1:2:n);
% Rotate Eigenvectors
phi=angle(Vs(1,:));
ij = 0.0000 + 1.0000i;
phi=diag(exp(-ij*phi));
Vs=(Vs*phi);
% Vs=flip(Vs*phi);
% Vs2 = Vs*phi;
% extract modal parameters from state space identification
ws = abs(log(Dsd)/(1/fs));
d = (-real(log(Dsd)/(1/fs))./ws)*100;
% d2 = (-real(log(Dsd)/(1/fs))./ws)*100;
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

