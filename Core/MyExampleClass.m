classdef MyExampleClass < handle
 
   properties
      Data % [Ndata,Nchannels]
      Name % object's name
   end
     
   properties (SetAccess = protected, Hidden)
      h = struct('figTime',[],'hax',[],'UI',[]) % GUI handles
      hl % listeners' handles
   end
   
   events (ListenAccess = private)
      DataChanged % set.Data notifier
   end
   
   methods 
      
      function obj = MyExampleClass % constructor
         fprintf('%s object created\n',mfilename)
         obj.hl = addlistener(obj,'DataChanged',@obj.updatePlot); % plot update when .Data is set         
      end
      
      function set.Data(obj,val) % set method override
         % Assign value to .Data
         obj.Data = val;
         notify(obj,'DataChanged') % notify all listeners 
      end
      
      function plot(obj)
         % Plot .Data
         if isempty(obj.Data)
            warning('The property .Data is empty, use .generateData')
         else
            if isempty(obj.h.figTime) % create figure if none is available    
               obj.Data = obj.Data; % update plot (set method overridden for illustration purposes)
            end 
         end
      
      end
            
      function generateData(obj,~,~)
         % Generate example data         
         obj.Data = obj.calcTimeData(1E4,8);
         disp('Example data has been generated')
         
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
            obj.h.UI.plotType = uicontrol('Style','popupmenu','Position',[10 10 80 20],...
               'String',{'Data plot','Spectrum'},'FontSize',10,...
               'Callback',@obj.createPlot);
            obj.h.UI.updatePlot = uicontrol('Style','pushbutton','Position',[10 35 80 30],...
               'String','Update','FontSize',10,...
               'Callback',@obj.generateData);
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
                  plot(obj.h.hax,obj.Data,'-r')
                  obj.h.hax.YScale = 'lin';
                  title('Data plot'), xlabel('Index [num]'), ylabel('Data [unit]'), legend('Data')
               case 2 % plot linear spectrum
                  F = obj.calcLinearSpectrum(obj.Data); % calculate spectrum
                  plot(obj.h.hax,abs(F),'-b')
                  obj.h.hax.YScale = 'log';
                  title('Linear spectrum'), xlabel('Frequency bin [num]'), ylabel('|FFT| (unscaled) [unit/1]'), legend('FFT')
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
      
      function [out,Out] = calcTimeData(L,fs)
         % Generate example data         
         M = [2 0 0;0 1 0;0 0 1];
         D = [0.15 -0.10 0;-0.10 0.10 0;0 0 0]*5;
         K = [80 -50 0;-50 120 -70;0 -70 70];

         DoF = size(M,1);
         s = 2i*pi*(0:L-1)'/L*fs/2;
         H = zeros(DoF,DoF,numel(s));
         for n = 1:numel(s)
            Z = s(n).^2*M+s(n)*D+K; 
            H(:,:,n) = Z\eye(DoF);
         end
         idx = [1 2 3 5 6 9];
         H = reshape(H,[],size(H,3)).'; 
         Y = H(:,idx(randperm(6,1))).*exp(2i*pi*rand(numel(s),1));
         
         for jj = 1:3
             Yj(:,jj) = H(:,jj).*exp(2i*pi*rand(numel(s),1));
             Out(:,jj) = real(ifft([Yj(:,jj) ; conj(Yj(end:-1:2,jj))]));
         end
         
         
         out = real(ifft([Y ; conj(Y(end:-1:2,:))])); % Output
         
      end
      
   end
   
end