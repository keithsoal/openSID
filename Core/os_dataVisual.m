% os_dataVisual generates a class for data visualisation
% openSID
% Based on an example by Goran Jelicic (DLR)
% Keith Soal    01-08-2017
% version 8 September 2016

% devnotes
% 1. clean up

classdef os_dataVisual < handle
    
    properties
        TData % [Ndata,Nchannels]
        Name % object's name
        fs % sample frequency
        nfft %nfft points ds = fs/nfft
        Nw % length of smoothing window
        FData % spectra
        fvec % frequency vector
        S % stats
    end
    
    properties (SetAccess = protected, Hidden)
        h = struct('figTime',[],'hax',[],'UI',[]) % GUI handles
        hl % listeners' handles
    end
    
    events (ListenAccess = private)
        DataChanged % set.Data notifier
    end
    
    methods
        
        function obj = os_dataVisual % constructor
            fprintf('%s object created\n',mfilename)
            obj.hl = addlistener(obj,'DataChanged',@obj.updatePlot); % plot update when .Data is set
        end
        
        function set.TData(obj,val) % set method override
%             V.TData=td
            % Assign value to .Data
            if isempty(obj.TData)
                obj.TData = val;
                notify(obj,'DataChanged') % notify all listeners
            else
                notify(obj,'DataChanged')
            end
                
        end
        
        function plot(obj)
            % Plot .Data
            if isempty(obj.TData)
                warning('The property .Data is empty, use .generateData')
            else
                if isempty(obj.h.figTime) % create figure if none is available
                    obj.TData = obj.TData; % update plot (set method overridden for illustration purposes)
                end
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function updatePlot(obj,~,~) % cbk event listener
            obj.createFigure
            obj.createPlot
        end
        
        function createFigure(obj)
            % Create GUI for simple data display
            if isempty(obj.h.figTime) % create figure if none is available
                obj.h.figTime = figure('Name',[mfilename,'''s plot'],'Position',[100 100 800 400],...
                    'CloseRequestFcn',{@int_closefigTime,obj}); % assign callback when closing figure
                obj.h.hax = axes; % create axes
                obj.h.UI.plotType = uicontrol('Style','popupmenu','Position',[10 10 100 20],...
                    'String',{'Time Data','Spectrum','Statistical Moments'},'FontSize',10,...
                    'Callback',@obj.createPlot);
%                 obj.h.UI.updatePlot = uicontrol('Style','pushbutton','Position',[10 35 80 30],...
%                     'String','Update','FontSize',10,...
%                     'Callback',@obj.generateData);
            end
            
            function int_closefigTime(hObject,~,obj) % Manage closing the data figure
                delete(hObject) % delete data figure (hObject)
                obj.h.figTime = []; % clean the handles' structure
            end
            
        end
        
        function createPlot(obj,~,~)
            % Plot data or spectrum according to the popupmenu's value
            if ~isempty(obj.h.figTime)
                switch obj.h.UI.plotType.Value
                    case 1 % plot data
                        subplot(1,1,1)
                        LW = 1.5; % linewidth
                        M = rms(obj.TData);
                        [Y,I] = sort(M,'descend');
                        Cr = linspecer(size(obj.TData,2),'red'); 
                        for pp = 1:size(obj.TData,2)
                        plot(obj.TData(:,I(pp)),'color',Cr(pp,:),'LineWidth',LW)
                        hold on
                        end
                        hold off
                        title('Time Data','Interpreter','latex')
                        xlabel({'Index (no.)'},'Interpreter','latex')
                        ylabel('Acceleration ($m/s^2$)','Interpreter','latex')
                        axis tight
                        set(gca,'fontsize',12)
                        grid on
                    case 2 % plot linear spectrum
                        LW = 1.5;
                        if isempty(obj.fs) 
                        obj.fs = input('Input sample frequency:');
                        end
                        if isempty(obj.nfft)
                        obj.nfft = input('Input nfft points:');
                        end
                        w = hann(obj.nfft);
                        if isempty(obj.FData)
                            [ap3,f_vec] = obj.cpsdrefk(obj.TData,[],[],w,obj.nfft,obj.fs);
                            ap3 = reshape(ap3,[],size(ap3,3)).';
                            obj.Nw = 2; % length of smoothing window
                            a=1;
                            b=1/obj.Nw*ones(1,obj.Nw);
                            [N,Dd,R]=size(ap3);
                            Ps=zeros(N,Dd);
                            for d = 1:Dd
                                for r=1:R
                                    Ps(:,d,r)=filtfilt(b,a,ap3(:,d));
                                end
                            end
                            obj.FData = Ps;
                            obj.fvec = f_vec;
                        else
                            'continue';
                        end                       
                        
                        subplot(1,1,1)
                        C = linspecer(size(obj.FData,2),'red');
                        for pp = 1:size(obj.FData,2)
                        semilogy(obj.fvec,abs(obj.FData(:,pp)),'color',C(pp,:),'LineWidth',LW)
                        hold on
                        end
                        hold off
                        title('Power Spectral Density (PSD)','Interpreter','latex')
                        xlabel({'Frequency (Hz)'},'Interpreter','latex')
                        ylabel('Acceleration ($(m/s^2)^2 / Hz$)','Interpreter','latex')
                        grid on
                        set(gca,'fontsize',12)        
                    case 3 % statistical moments
                        LW = 1.5;
                        if isempty(obj.S)
                        ss = obj.statsMoments(obj.TData);
                        obj.S = ss;
                        else
                            'continue';
                        end
                        
%                         figure()
                        % figure('units','normalized','position',[.1 .1 0.44 0.53])
                        C = linspecer(size(obj.TData,2),'red');
                        axes('NextPlot','replacechildren', 'ColorOrder',C);
                        subplot(2,2,1)
                        hold off;
                        for ii = 1:size(C,1)
                            plot(obj.S.mm(:,ii),'color',C(ii,:),'LineWidth',LW)
                            hold on;
                        end
                        xlabel({'Index (no.)'},'Interpreter','latex')
                        ylabel('Mean','Interpreter','latex')
                        grid on
                        subplot(2,2,2)
                        hold off;
                        for ii = 1:size(C,1)
                            plot(obj.S.v(:,ii),'color',C(ii,:),'LineWidth',LW)
                            hold on;
                        end
                        xlabel({'Index (no.)'},'Interpreter','latex')
                        ylabel('Variance','Interpreter','latex')
                        grid on
                        subplot(2,2,3)
                        hold off;
                        for ii = 1:size(C,1)
                            plot(obj.S.s(:,ii),'color',C(ii,:),'LineWidth',LW)
                            hold on;
                        end
                        xlabel({'Index (no.)'},'Interpreter','latex')
                        ylabel('Skewness','Interpreter','latex')
                        grid on
                        subplot(2,2,4)
                        hold off;
                        for ii = 1:size(C,1)
                            plot(obj.S.k(:,ii),'color',C(ii,:),'LineWidth',LW)
                            hold on;
                        end
                        xlabel({'Index (no.)'},'Interpreter','latex')
                        ylabel('Kurtosis','Interpreter','latex')
                        grid on
                end
            end
        end
        
    end
    
    methods (Static)
        
        function out = calcLinearSpectrum(val)
            % Help function: calculates spectrum of data
            out = fft(val)/size(val,1);
            out = out(1:round(size(out,1)/2),:);
        end
        
        
        function [ap3,f_vec] = cpsdrefk(y,resp,ref,w,nfft,fs)
            % Efficient cross spectral density estimation - with reference channels
            % This algorithm uses the welch method with convolution and symmetry to
            % improve computational time
            % *note 50% window overlap is implemented
            % Keith Soal - PhD
            % define paramters
            NFFT=nfft/2+1;
            f_vec=linspace(0,fs/2,NFFT);
            delta_freq = fs/nfft;
            
            % Window parameters (50% overlap)
            win_size = length(w);
            overlap = length(w)/2;
            step_size = win_size - overlap;
            offset = 1:step_size:length(y)-win_size+1;
            m = length(offset);
            n = length(w);
            sy = size(y);
            
            % tic
            % preallocate space for the fft and psd results
            if isempty(resp) && isempty(ref)
                AP3 = zeros(sy(2),sy(2),nfft);
            else
                AP3 = zeros(size(ref,2),size(resp,2),nfft);
            end
            temp = zeros(nfft,sy(2));
            % loop through the data, one window at a time
            for i = 1:m
                for k = 1:sy(2)
                    % apply window
                    temp(:,k) = y(offset(i):offset(i)+win_size-1,k).*w;
                    % calculate fft
                    temp(:,k) = fft(temp(:,k));
                end
                % CPSD using all responses and references
                if isempty(resp) && isempty(ref)
                    % build lower triangular crosspower and sum blocks
                    for nn=1:sy(2)
                        for mm=1:nn
                            % convolution in frequency domain
                            % producing a periodogram of "absolute value squared"
                            ap=temp(:,nn).*conj(temp(:,mm));
                            AP3(nn,mm,:)=AP3(nn,mm,:)+(reshape(ap(1:nfft),1,1,nfft));
                            % use hermetian symmetry of matrix to populate the
                            % upper triangular with the complex conjugate
                            if mm ~= nn
                                ap_conj = conj(ap);
                                AP3(mm,nn,:)=AP3(mm,nn,:)+(reshape(ap_conj(1:nfft),1,1,nfft));
                            end
                        end
                    end
                else % CPSD using reference subset
                    % build lower triangular crosspower and sum blocks
                    for nn=1:size(ref,2)
                        for mm=1:size(resp,2)
                            % convolution in frequency domain
                            % producing a periodogram of "absolute value squared"
                            ap=temp(:,ref(nn)).*conj(temp(:,resp(mm)));
                            AP3(nn,mm,:)=AP3(nn,mm,:)+(reshape(ap(1:nfft),1,1,nfft));
                            % use hermetian symmetry of matrix to populate the
                            % upper triangular with the complex conjugate
                            %                     if mm ~= nn
                            %                         ap_conj = conj(ap);
                            %                         AP3(mm,nn,:)=AP3(mm,nn,:)+(reshape(ap_conj(1:nfft),1,1,nfft));
                            %                     end
                        end
                    end
                end
                
            end
            
            % psd scaling factor (Brandt pg. 212)
            Sp = 1 / (n*delta_freq*sum(w.^2));
            % scale and average the periodograms
            ap3 = (Sp/m)*AP3;
            % select half spectra
            ap3 = ap3(:,:,1:nfft/2 + 1);
            % scale by 2 for redundatn nyquist frequencies except ignore DC and Nyquist value
            ap3(:,:,2:end-1) = ap3(:,:,2:end-1) * 2;
        end
        
        function [S] = statsMoments(var)
            % statistics
            B = 1000;  %20000; %230400; % 1 hour
            Bovl = 0.5;%.50;
            idx = 1:B;
            Nb = fix((size(var,1)-B)/fix(B*(1-Bovl))+1);
            for n = 1:Nb
            S.mm(n,:) = nanmean(var(idx,:));
            S.v(n,:) = nanvar(var(idx,:));
            S.s(n,:) = skewness(var(idx,:));
            S.k(n,:) = kurtosis(var(idx,:));
            idx = idx + fix(B*(1-Bovl));
            end
%                         Si.k(n,:) = kurtosis(Ivar(idx,:));
        end
        
        
    end
    
end