% os_stabilisation build stabilisation diagram and prompts user input
% based on sdiagram.m in Abravibe from Anders Brandt
% openSID
% Keith Soal    01-08-2017
% version 8 September 2016

% devnotes
% 1. clean up

function [fn,d] = os_stabilisation(Poles,MaxModes,Mif,modal)

% Go through all requested orders
NoPoles = 1:MaxModes;
   
% Hardcode parameters
fLim=0.001;  % pole is frequency stable
zLim=0.05;   % pole is damping stable
% Matlab colors
mc2 = [0.8500 0.3250 0.0980];
mcb = [0    0.4470    0.7410];
mcy = [0.9290    0.6940    0.1250];
nSign='x'; % Marker for new poles
fSign='bd'; % Marker for stable frequency
sSign='o'; % Marker for stable pole
MarkerSize=4; % Size of markers
NoPairChk=false;

% Normalize Mode Indicator Function
MIF = mag2db(Mif);
MifMin=round(min(MIF));
MIF = MIF - MifMin;
if max(MIF(:,1)) > 1
    MifMax=max(max(MIF));
    MIF=MIF/MifMax;
elseif max(MIF(:,1)) < 1
    MifMax=max(max(MIF));
    MIF=MIF/MifMax;
end
% frequency vector
f = linspace(0,20,size(Mif,1));

if length(f) ~= length(MIF)
    fmif = f(end-length(MIF)+1:end);
else
    fmif=f;
end

% OrderStep=1;

% Clean up the poles list for non-physical poles
% Remove real poles
for on=1:length(Poles) % Order number
    pidx=1;
    for n=1:length(Poles{on})
        if ~isreal(Poles{on}(n))
            np{on}(pidx)=Poles{on}(n); % np is new pole vector
            pidx=pidx+1;
        end
    end
end

% Remove poles with positive real part
for on=1:length(np)
    idx=find(real(np{on}) < 0);
    np{on}=np{on}(idx);
end

% Now find complex conjugate pairs and save only positive frequencies
if NoPairChk == false
    clear Poles
    for on=2:length(np)
        pidx=1;
        for n=2:length(np{on})
            % Check distance between all higher poles with this model order and
            % the conjugate of current pole (current = n-1). If this distance
            % is small there is a complex conjugate, so save the positive pole
            dist=abs((np{on}(n:end)-conj(np{on}(n-1))));
            if ~isempty(find(dist < 1e-4))
%                 fprintf('found complex pole %f\n',imag(np{on}(n-1)))
                Poles{on}(pidx)=real(np{on}(n-1))+j*abs(imag(np{on}(n-1)));
                pidx=pidx+1;
            end
        end
    end
else
    Poles=np;
end
% Remove empty vectors from Poles
pidx=1;
for n=1:length(Poles)
    if ~isempty(Poles{n})
        np{pidx}=Poles{n};
        pidx=pidx+1;
    end
end
Poles=np;
clear np

% Plot MIF function and hold plot for poles to be plotted
h=figure;
plot(fmif,mag2db(MIF*max(NoPoles)),'color',mcy);
axis([min(fmif) max(fmif) 0 max(NoPoles)])
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('Number of Poles','Interpreter','latex')
title('Stabilization Diagram','Interpreter','latex')
grid on
set(gca,'fontsize',14)
hold on

% Normalise singular values
hold on
if max(modal.ss(1)) > 1
    ssMax=max(max(modal.ss(1)));
    ssScaled=modal.ss/ssMax;
elseif max(modal.ss(1)) < 1
    ssMin=max(max(modal.ss(1)));
    ssScaled=modal.ss/ssMin;
end
% plot singular values
hold on
for i = 1:length(modal.ss)
    bh = barh(ssScaled.*(f(end)*(3/4)), 'facecolor', mcy);
    hold on
end

% Check stability and plot symbols accordingly
fOffset = 0;
NPoles=NoPoles(1);
LastRow=Poles{1};               % First row stored
wp=abs(LastRow);                % Previous freqs
zp=-real(LastRow)./wp;
plot(wp/2/pi+fOffset,NPoles,nSign,'color',mc2,'MarkerSize',MarkerSize)
hold on
for n = 2:length(NoPoles)                % Each row (model order) in diagram
    NPoles=NoPoles(n);
    for m = 1:length(Poles{n})   % Each pole m from order n
        wr=sqrt(abs(Poles{n}(m))^2);
        zr=-real(Poles{n}(m))/abs(Poles{n}(m));
        % See if the pole is within limits of a pole in last row
        fDist=abs(wr-wp)./abs(wp);
        zDist=abs(zr-zp)./abs(zp);
        if min(fDist) < fLim & min(zDist) < zLim    % Stable pole
            plot(wr/2/pi+fOffset,NPoles,sSign,'color',mcb,'MarkerSize',MarkerSize)
        elseif min(abs(wr-wp)/abs(wp)) < fLim       % Stable frequency
            plot(wr/2/pi+fOffset,NPoles,fSign,'MarkerSize',MarkerSize)
        else
            plot(wr/2/pi+fOffset,NPoles,nSign,'color',mc2,'MarkerSize',MarkerSize)       % New pole
        end
    end
    LastRow=Poles{n};               % First row stored
    wp=abs(LastRow);     % Previous freqs
    zp=-real(LastRow)./abs(LastRow);
end

% Now have user select poles by clicking in vicinity of markers on plot
disp('Select poles!')
disp('Then <RETURN>')
[xx,yy]=ginput;
pidx=1;
for n = 1:length(xx)
    [dum,PoleIdx]=min(abs(yy(n)-NoPoles));
    [dum,idx]=min(abs(imag(Poles{PoleIdx})-2*pi*xx(n)));
    p(pidx)=Poles{PoleIdx}(idx);
    pidx=pidx+1;
    N{n}=NoPoles(PoleIdx);
end
% Sort and force poles to a column
% p=sort(p);              % Does not matter which order poles were selected
p=p(:);

ws = abs(p);
d = -real(p)./ws;
fn = ws/(2*pi);

disp('Frequency (fn)  Damping (z)')
disp([fn  d*100])



