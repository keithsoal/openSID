%----------------------------------------------------------------------
%                          openSID demo 
%----------------------------------------------------------------------
% This function runs through a demonstration of the openSID toolbox on
% simulated data.
% Keith Soal
% version 8 September 2016
clear;close all;clc;

% Generate simulated 3DOF data
m = [100 100 100];
k = [100000 100000 100000];
z = [0.1 0.12 0.14]./100;
[out,Y] = os_generateData(m,k,z,1);

% Visualize data
Dv = os_dataVisual;
Dv.TData = out;

% SSI-Data
i = 10; % block size
fs = 40; % sample frequency
Ir = []; % reference channels
W = 'UPC'; % weighting
n = 6; % model order
% precompute sys.AUX for increased speed in stabilisation diagram
[sys,modal] = os_ssid(out,i,fs,Ir,W,n);

% Mode indicator function
Type = 'power';
Mif = amif(Y,Type);

% Stabilisation diagram
% Go through all requested orders
MaxModes = 30;
NoModes=[];
for Order=1:MaxModes
    NoModes=[NoModes Order];
    [sys,modal] = os_ssid(out,i,fs,Ir,W,Order,sys.AUX,1);
    Poles{Order}=modal.EW;
end
[fn,d] = os_stabilisation(Poles,MaxModes,Mif,modal,fs);


% SSI-Data
i = 10; % block size
fs = 40; % sample frequency
Ir = []; % reference channels
W = 'UPC'; % weighting
n = 6; % model order
% precompute sys.AUX for increased speed in stabilisation diagram
[sys,modal] = os_ssid(out,i,fs,Ir,W,n);
%--------------------------------------------------------------------------
%                               MAC
%--------------------------------------------------------------------------
M = amac(modal.EV(:,1:2:end),modal.EV(:,1:2:end));
figure()
Mm = size(M,2); % grid resolution
x = linspace(1,Mm,Mm); % x-grid
y = linspace(1,Mm,Mm); % y-grid
os_bar3c(M,x,y,[],jet(50))
colorbar
colormap(jet(50))
view(90,90);
title('MAC Matrix','Interpreter','latex')
set(gca,'fontsize',14)

%--------------------------------------------------------------------------
%                               Mode shapes
%--------------------------------------------------------------------------
phi = modal.EV(:,1:2:end);
% Eigen vector plotting
figure(4)
set(gcf,'DefaultLineLineWidth',2)
for nn=1:size(phi,1)
    subplot(size(phi,1),1,nn)
    hold on; plot(real(phi(:,nn)));
    title(['Mode Shape' num2str(nn)],'Interpreter','latex')
    xlabel('node point','Interpreter','latex')
    ylabel('amplitude','Interpreter','latex')
    grid on
    axis tight
    set(gca,'fontsize',14)
end 

%--------------------------------------------------------------------------
%                               Complexity plot
%--------------------------------------------------------------------------
figure()
set(gcf,'DefaultLineLineWidth',2)
rp = real(phi);
ip = imag(phi);
for nn=1:size(phi,1)
    subplot(size(phi,1),1,nn)
    for i = 1:3
    hold on;plotv([rp(i,nn);ip(i,nn)]);
    end
    title(['Complexity plot mode' num2str(nn)],'Interpreter','latex')
    xlabel('real','Interpreter','latex')
    ylabel('imaginary','Interpreter','latex')
    grid on
    axis equal
    set(gca,'fontsize',13)
end 




