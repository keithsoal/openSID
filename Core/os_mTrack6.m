% Mode tracking algorithm 6.2
% Keith Soal

function [sys,modal] = os_mTrack6(p,Psi,mcpair,sys,modal)

if (nargin < 4);sys = [];modal = [];end

%--------------------------------------------------------------------------
%                           Base family
%--------------------------------------------------------------------------
if isempty(sys)
    leeg = 1; % flag indicates sys started empty
    MACXP = ae_macxp(Psi{1,1},p{1,1},Psi{1,2},p{1,2});
    
    % find macxp values greater than mcpair
    [I,J] = find(MACXP >= mcpair); % [I - rows] [J - columns]
    Is = size(J,1); % no. of unique families
    
    % find remaining eigenvector indeces (from Psi{1,2})
    [~,J2] = ismember(1:size(Psi{1,2},2),J);
    [~,J3] = find(J2==0); % J3 index of unmatched modes (*no family match)
    
    % build modal family
    % EMACS = Psi{1,2}(:,J);
    % EMAC = [EMACS Psi{1,2}(:,J3)];
    sys.Psi_library = Psi{1,2}(:,J); % size 23 x no. families
    sys.Psi_tot = [sys.Psi_library Psi{1,2}(:,J3)]; % size 23 x unidentified families
    
    % EpS = p{1,2}(J);
    % Ep = [EpS;p{1,2}(J3)];
    sys.p_library = p{1,2}(J);
    sys.p_tot = [sys.p_library ; p{1,2}(J3)];
    
    modal.FTrack(1:length(J),1) = abs(p{1,2}(J))/(2*pi); % frequency library
    modal.DTrack(1:length(J),1) = -real(p{1,2}(J))./(abs(p{1,2}(J)))*100; % damping library
else
    leeg = 0;
end

%--------------------------------------------------------------------------
%                           Builds family
%--------------------------------------------------------------------------
if leeg == 1
    m_start = 3; % uses p1 and p2 to compile first library
    countFD = 1;
else
    m_start = 1; % library already exists
    countFD = size(modal.FTrack,2);
    Is = size(sys.p_library,1); % no. of unique families
end

Psi_tot = sys.Psi_tot;
p_tot = sys.p_tot;

for m = m_start:size(Psi,2)
    clearvars MACXP2 I J J2 J3;
    MACXP2 = ae_macxp(Psi_tot,p_tot,Psi{1,m},p{1,m});

    % find macxp values greater than mcpair
    [I,J] = find(MACXP2 >= mcpair);
%     Is = size(J,1); % no. of unique families
    
    % find remaining eigenvector indeces
    [~,J2] = ismember(1:size(Psi{1,m},2),J);
    [~,J3] = find(J2==0);

    countFD = countFD + 1; % counter to build F and Dtrack
    % build modal family
    for nn = 1:size(J,1) % no. of unique families
        
        % if current mode is outside library (i.e. new)
        % add mode one after current library size (Is + 1)
        if I(nn) > Is
            I(nn) = Is + 1;
            Is = Is + 1;
        end
        
        sys.Psi_library(:,I(nn)) = Psi{1,m}(:,J(nn)); % place new J columns in existing I columns of matrix
        sys.p_library(I(nn),1) = p{1,m}(J(nn)); % new J@nn goes to I@nn
        
        
        modal.FTrack(I(nn),countFD) = abs(p{1,m}(J(nn)))/(2*pi);
        modal.DTrack(I(nn),countFD) = -real(p{1,m}(J(nn)))./(abs(p{1,m}(J(nn))))*100;
    end
    
    modal.FTrack(modal.FTrack==0)=NaN; % convert mising frequencies to nan
    modal.DTrack(modal.DTrack==0)=NaN;

    sys.Psi_tot = [sys.Psi_library Psi{1,m}(:,J3)];
    Psi_tot = sys.Psi_tot;
    sys.p_tot = [sys.p_library ; p{1,m}(J3)];
    p_tot = sys.p_tot;
end

% sys.p_library = p_library; % Eivenvalues
% sys.Psi_library = Psi_library; % Eigenvector
% sys.EW = EpS; % Eivenvalues
% sys.EV = EMACS; % Eigenvector


% modal.FTrack = FTrack; % frequency
% modal.DTrack = DTrack; % damping

end










