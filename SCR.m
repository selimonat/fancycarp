classdef SCR < handle
    
    properties (Constant = true, Hidden)
        default_run = 1;%at which run the scr data is located.
    end
    properties (Hidden)
        ledalab_defaults      = {'open', 'mat','downsample', 20, 'smooth',{'gauss' 100} 'analyze','CDA', 'optimize',10, 'overview',  1, 'export_era', [-1 7 0 1], 'export_scrlist', [0 1], 'export_eta', 1};%
        hdr
        markers
        block2phase
        p         = 0.0001;
        nfreq     = 20;
        FIR_delay = 15000;%in milliseconds
        path_scr;%path to the file
        plot_style
    end
    properties
        ledalab
        data
        fear_tuning
        ml
    end
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        y
        time
        tsample
        sampling_rate
        sampling_period
        BlockBorders_index
        BlockBorders_time
        BlockNames
        event           = [];
        event_name      = {};
        event_plotting  = {};
        tonic
        phasic
        model
    end
    
    methods
        %add a method for tonic removal a la ledalab and a la new baseline
        %correction
        %add a method to compute FIR analysis and plot the results based on
        %the corrected phasic responses.
        function scr = SCR(varargin)
            %s is either the subject number of block id.
            if length(varargin) == 1%construct
                %
                s                 = varargin{1};
                path_scr          = s.path_scr;
                path_scrmat       = regexprep(path_scr,'.smr','.mat');
                %
                if ~exist(path_scrmat)
                    a    = SONImport(fopen(path_scr),'milliseconds','scale');%will cache it.
                    if a == 0
                        fprintf('SON Import OK.\n');
                    else
                        fprintf('SON Import failed!\n')
                        keyboard
                    end
                end
                %%
                a                      = load(path_scrmat);
                data                   = a.chan3; %chan1 is HR, chan2 is breathing, chan3 is scr
                time                   = [a.head3.start:1:a.head3.stop]';
                %discard all time before the first and after the last pulses.
                pulse_times            = a.chan7;
                stim_times             = a.chan9;    %Face Onsets
                mbi_times              = a.chan2;                
                % make a robust fit for pulse times to find the intercept.
                b                      = robustfit([ [1:length(pulse_times)]'],pulse_times);%b(1) is an estimation of the time experiment starts.                
%                 plot(stim_times,'o-'),hold on;plot(pulse_times,'ro-');plot(mbi_times,'ko-');hold off;grid on
                b                      = b(1);%start of the presentation
                %
                stim_times( stim_times < b)  = [];%delete all past samples
                pulse_times( pulse_times< b) = [];%delete all past samples
                mbi_times(mbi_times < b)     = [];%delete all past samples
                %remove the last 16 stimuli 
                stim_times(end-15:end)       = [];%delete all past samples
                %if the first stimulus onset is missing, replace it with
                %the mbi onset (which is the same)
                if abs(stim_times(1)-mbi_times(1)) > 500;
                    stim_times = [mbi_times(1) ; stim_times];
                end
                %check for ghosts in the stim timings
                i = [diff(stim_times)] < 2000;
                i = [0;i];
                if any(i);
                    cprintf([1 0 0],'Removed %d ghost pulses from stim onsets...\n',sum(i));
                    stim_times(i==1) = [];
                end
                %check for ghosts in the stim timings
                i = [diff(mbi_times)] < 20000;
                i = [0;i];
                if any(i)
                    cprintf([1 0 0],'Removed %d ghost pulses from mbi...\n',sum(i));
                    mbi_times(i==1) = [];
                end                
%                 figure;
%                 plot(stim_times,'o-'),hold on;plot(pulse_times,'ro-');plot(mbi_times,'ko-');hold off;grid on
                %sanity checks
                if length(mbi_times) == 65
                    cprintf([0 1 0],'MBI count, correct.\n');
                else
                    cprintf([1 0 0],'MBI count, incorrect.\n');
                end
                if length(stim_times) == 585
                    cprintf([0 1 0],'STIM count, correct.\n');
                else
                    cprintf([1 0 0],'STIM count, incorrect.\n');
                end                
                %% now set the zero for the scr recordings
                % first cut the unwanted parts based on first and last mbi
                last_point  = max(stim_times);
                first_point = min(mbi_times);
                i           = (time >= (first_point - 15000)) & (time <= (last_point + 20000));
                time(~i)    = [];
                data(~i)    = [];
                % now shift everything so that first sample is zero.
                first_sample                            = time(1);
                time                                    = time        - first_sample;
                pulse_times                             = pulse_times - first_sample;
                mbi_times                               = mbi_times   - first_sample;
                stim_times                              = stim_times  - first_sample;
                
                %% store the stuff
                scr.hdr                = a.FileInfo;
                scr.time               = time;
                scr.y                  = data;
                scr.y                  = [scr.y    ; repmat(scr.y(end),25000,1)];
                scr.time               = [scr.time ; repmat(scr.time(end),25000,1) + [1:25000]'];                
                scr.tsample            = length(time);
                scr.sampling_period    = (a.head1.sampleinterval/1000);%in ms
                scr.sampling_rate      = 1000/scr.sampling_period;%in hertz.
                scr.BlockBorders_time  = time([1 end])';
                scr.BlockBorders_index = [1 length(scr.y)];
                scr.block2phase        = 1;
                scr.BlockNames         = {'all'};
                scr.path_scr           = path_scr;
                %
                c = 0;
                scr.event = [];
                s.paradigm{1}.presentation.dist(ismember(s.paradigm{1}.presentation.dist, [1001 1002])) = 0;
                for distance = unique(s.paradigm{1}.presentation.dist)                    
                    c                  = c + 1;
                    scr.event(:,c)     = zeros(length(scr.y),1); 
                    ii                 = find(s.paradigm{1}.presentation.dist == distance);
                    
                    %transform distance to cond id
                    if distance == -135
                        cond_id = 1;
                    elseif distance == -90
                        cond_id = 2;
                    elseif distance == -45
                        cond_id = 3;
                    elseif distance == 0
                        cond_id = 4;
                    elseif distance == 45
                        cond_id = 5;
                    elseif distance == 90
                        cond_id = 6;
                    elseif distance == 135
                        cond_id = 7;
                    elseif distance == 180
                        cond_id = 8;
                    elseif distance == 1000
                        cond_id = 9;
                    elseif distance == 1001
                        cond_id = 10;
                    elseif distance == 1002
                        cond_id = 11;
                    end
                    
                    event_times                          = stim_times(ii);
                    event_index                          = event_times;
                    scr.event(round(event_index),c)      = true;
                    scr.event_name{c}                    = sprintf('%03d',cond_id);
                    scr.event_plotting{c}.line_width     = {{'linewidth', 2}};
                    scr.event_plotting{c}.marker_size    = {{'markersize', 2}};
                    scr.event_plotting{c}.symbol         = {{'symbol','+'}};
                    scr.event_plotting{c}.line           = {{'line',{'line' '-'}}};
                    scr.event_plotting{c}.color          = {{'color',Project.colors(:,cond_id)}};
                end
                scr.event = logical(scr.event);
                
                
                
                
            elseif length(varargin) == 2 %cut
                %% will extract SCR data of a given block or blocks
                block                     = varargin{2};
                scr                       = varargin{1};
                min_i                     = min(scr.BlockBorders_index(block,1));
                max_i                     = max(scr.BlockBorders_index(block,2));
                max_i                     = min( max_i + (10*1000)*scr.sampling_period , scr.tsample);%add 10 seconds, but not if it exceeds the experiments length
                i                         = min_i:max_i;
                scr.y                     = scr.y(i);
                scr.time                  = scr.time(i);
                scr.tsample               = length(scr.y);
                scr.BlockBorders_index    = scr.BlockBorders_index(block,:);
                scr.BlockBorders_time     = scr.BlockBorders_time(block,:);
                scr.BlockNames            = scr.BlockNames(block);
                scr.event                 = scr.event(i,:);
                try
                    scr.tonic                 = scr.tonic(i);
                end
                try
                    scr.phasic                = scr.phasic(i);
                end
                try
                    scr.data                  = scr.data(i);
                end
            end
        end
    end
    methods
        function out = cut(self,block)
            out = SCR(self,block);
        end
        function out = findphase(self,filter_string);
            out = cellfun(@(x) ~isempty(regexp(x,filter_string)), self.BlockNames );
            out = find(out);
        end
        function plot_raw(self)
        end
        function plot(self)
            %% plot different event channels
            
            hhfigure;
            % plot the SCR
            plot(self.time./1000,self.y,'color','k','linewidth',1);
            hold on;
            % mark stimulus onsets
            level = mean(self.y);
            for n = 1:length(self.event_name)
                i  = self.event(:,n);
                plot(self.time(i)./1000,self.y(i), self.event_plotting{n}.symbol{1}{2} , self.event_plotting{n}.color{1}{:} , self.event_plotting{n}.marker_size{1}{:} );
                if ~isempty(self.phasic)
                    plot(self.time(i)./1000,self.phasic(i), self.event_plotting{n}.symbol{1}{2} , self.event_plotting{n}.color{1}{:} , self.event_plotting{n}.marker_size{1}{:} )
                end
            end
            
            %% phase transitions as lines
            i   = self.BlockBorders_time(:,1)./1000;%start of phases
            plot( repmat(i,1,2)', repmat(ylim,size(i,1),1)' , '-','color','k');
            %% put a # to different blocks
            cc = 0;
            for ni = 1:length(i);
                cc = cc +1;
                text( i(ni), max(ylim), mat2str(cc),'fontsize',25);
            end
            %% labels
            hold off
            xlabel('time (ms)');
            ylabel('scr');
            box off;
            legend boxoff
            axis tight;
        end
        function plot_decomposition(self)
            self.plot
            hold on;
            plot(self.time./1000,self.phasic,'k');
            plot(self.time./1000,self.tonic,'color',[.4 .4 .4]);
            plot(self.time./1000,self.tonic+self.phasic,'r');
            grid on;
            axis tight;
            hold off
        end
        function self = smooth(self,method)
            %
            if isempty(self.data);self.data=self.y;end
            
            if strcmp(method,'ma');%moving average
                self.data  = smooth(self.data,750);
            elseif strcmp(method,'sgolay');%moving average
                self.data  = sgolayfilt(self.data,5,751);
            end
        end
        function self = diff(self)
            %takes the temporal difference.
            if isempty(self.data);self.data=self.y;end
            self.data         = diff(self.data);
            self.data(end+1)  = self.data(end);
        end
        function self = tonic_lowpass(self)
            %computes the low-pass version of the signal supposed to
            %represent the tonic response.
            if isempty(self.data);self.data=self.y;end
            d            = fdesign.lowpass('Fp',0.00001);
            Hd           = design(d, 'ellip');
            self.tonic   = filtfilt(Hd.sosMatrix,Hd.ScaleValues,self.data);
            self.phasic  = self.data - self.tonic;
        end
        function xcorr(self,block)
            %%
            self.data = [];
            figure(2);clf;
            cc = 0;
            self.smooth('sgolay');
            self.diff;
            for conds = {1:11 12:16 17:27 28:29};
                cc = cc +1;
                subplot(1,4,cc)
                for n = conds{1}
                    [r lag] = xcorr(self.data,double(self.event(:,n)),round(self.sampling_rate*10));
                    plot(lag(lag>0),r(lag>0),self.event_plotting{n}.line_width{1}{:},self.event_plotting{n}.color{1}{:});
                    hold on;
                    drawnow;
                end
                xlim([0 max(lag)]);
                %add a legend
                for nn = 1:length(conds{1})
                    dummy = regexp(self.event_name{nn},'[0-9]*$','match');
                    leglab{nn} = dummy{1};
                end
                l = legend(leglab,'interpreter','none');
                set(l,'box','off','fontsize',12,'location','westoutside');
                box off;
            end
            hold off
        end
        function linespecs = LineSpecs(condition)
            %produce a string for each condition that specifies plotting
            %attributes.
            for ncond = [1:10];
                linespecs{condition_indices{ncond}} = {linetype{ncond},'color',colors(ncond,:)};
            end
            linespecs = linespecs(condition);
        end
        function downsample(self,n)
            self.sampling_period    = self.sampling_period*n;
            self.sampling_rate      = self.sampling_rate/n;
            %will downsample the object instance by N
            self.y                  = decimate(self.y,n);
            self.time               = [(0:length(self.y)-1)*self.sampling_period]'+self.time(1);
            self.tsample            = length(self.y);
            self.BlockBorders_index = self.BlockBorders_index./n;
            self.BlockBorders_time  = self.BlockBorders_time;
            %
            tevent                  = size(self.event,2);
            tsample                 = size(self.y,1);
            new_mat                 = logical(zeros(tsample,tevent));
            for event = 1:tevent
                i                = round(find(self.event(:,event))./n);
                new_mat(i,event) = 1;
            end
            self.event = new_mat;
            try
                self.data               = downsample(self.data,n);
            end
            try
                self.tonic               = downsample(self.tonic,n);
            end
            try
                self.phasic               = downsample(self.phasic,n);
            catch
                fprintf('Downsampling failed, likely cause Ledalab''s downsample shadows it\n');
            end
        end
        function [FIRmat,FIRtime] = FIR(self,events_i)
            %produces a FIR matrix using event index
            e       = self.event(:,events_i);
            tcond   = size(e,2);
            tshift  = self.FIR_delay./self.sampling_period;%number of sample shifts
            FIRtime = [0:(tshift-1)]*self.sampling_period;
            FIRmat  = logical(zeros(self.tsample,tcond*tshift));
            counter = 0;
            for nc = 1:size(e,2)
                for shift = 1:tshift
                    tshift*(nc-1)+shift;
                    counter                   = counter + 1;
                    FIRmat(shift:end,counter) = e(1:end-shift+1,nc);
                end
            end
        end
        function [P]      = FourierBasis(self,tsample)
            %
            if nargin < 2
                tsample = self.tsample;
            end
            %% Normalized Fourier basis sets
            t     = 0:tsample-1;
            nb    = self.nfreq*2+1;
            z     = zeros(nb,tsample);
            z(1,:)= ones(1,tsample);
            fs    = 0.25;
            for i=1:self.nfreq
                z(2*i,:)  = cos((2*pi*fs*t)/tsample);
                z(2*i+1,:)= sin((2*pi*fs*t)/tsample);
                if fs < 0.5
                    fs = 0.25+fs;
                elseif fs == 0.5;
                    fs = 0.5+fs;
                else
                    fs = 1+fs;
                end
                % Calculate number of data points in xab
                % calculate basis functions % calculate basis functions
            end
            for m=1:nb;
                z(m,:)=z(m,:)/norm(z(m,:));
            end
            % Orthonormal Basis sets
            [~ , ~, P]=svd((z),'econ');
        end
        function tonic_tfals(self)
            % Normalize basis functions
            % Orthogonalize basis functions
            % Generate sparse matrix of asymmetric weights
            % perform least squares fit % Estimate baseline
            % Calculate baseline corrected signal vector
            %see TFALS.m for more details.
            %%
            if ~isempty(self.data)
                %
                p      = 0.001;   %p=Asymmetry parameter (0.001>=p<=0.1) %
                Y      = self.data;
                pa     = round(self.tsample/2);%pad amount
                Y      = padarray(Y,[pa 0],'replicate','both');
                tsample= self.tsample + pa*2;
                P      = self.FourierBasis(tsample);
                %% Asymmetric least squares
                w      = ones(tsample,1);
                target = 1;
                e      = [];
                viz    = 0;
                while target
                    W          = spdiags(w,0,tsample,tsample);
                    bw         = (P'*W);
                    q          = (bw*P)\(bw*Y);
                    self.tonic = P*q;
                    w0         = w;
                    w(Y>(self.tonic))    = p;
                    w(Y<=(self.tonic))   = (1-p);
                    target               = sum(abs(w - w0)) > 0;
                    if viz
                        e = [e sum(abs(Y - self.tonic))];
                        subplot(1,2,1);plot(e);
                        subplot(1,2,2);plot(self.tonic);
                        hold on;
                        plot(Y,'r');hold off;
                        drawnow
                        pause(.1)
                    end
                end
                %clean the flankers
                self.tonic(1:pa)          = [];
                self.tonic(end-pa+1:end)  = [];
                %get the phasic
                self.phasic               = self.data-self.tonic;
                
            else
                fprintf('.data is empty honey\n')
            end
        end
        function asymmetricls_decomposition_diff(self,event_i)
            %             self.smooth('sgolay');
            self.data = self.y;
            self.tonic_tfals;%will provide both the tonic fit, the residuals are phasic.
            %model the phasic responses using FIR.
            [FIR self.model.FIRtime]= self.FIR(event_i);
            %
            self.phasic             = diff(self.phasic);
            self.phasic(end+1)      = self.phasic(end);
            %
            self.model.betas        = FIR\self.phasic;%now model the phasic response
            %             self.model.betas        = cumsum(reshape(self.model.betas,[150 11]));
            %             self.model.betas        = self.model.betas(:);
            self.model.fit          = cumsum(FIR*self.model.betas) + self.tonic;%fit = fit_tonic + fit_phasic
            self.model.fit_phasic   = cumsum(FIR*self.model.betas);
            self.model.fit_tonic    = self.tonic;
            self.model.r2           = corr2(self.model.fit, self.data);
            self.model.r2_phasic    = corr2(self.phasic   , FIR*self.model.betas);
        end
        function asymmetricls_decomposition(self,event_i)
            self.tonic_tfals;%will provide both the tonic fit, the residuals are phasic.
            %model the phasic responses using FIR.
            [FIR self.model.FIRtime]= self.FIR(event_i);
            self.model.betas        = FIR\self.phasic;%now model the phasic response
            self.model.fit          = FIR*self.model.betas + self.tonic;%fit = fit_tonic + fit_phasic
            self.model.fit_phasic   = FIR*self.model.betas;
            self.model.fit_tonic    = self.tonic;
            self.model.r2           = corr2(self.model.fit, self.data);
            self.model.r2_phasic    = corr2(self.phasic   , FIR*self.model.betas);
        end
        function plot_model(self)
            
            self.plot;
            hold on;
            plot(self.time./1000,self.model.fit_phasic,'r');
            plot(self.time./1000,self.phasic,'k');
            plot(self.time./1000,self.model.fit_tonic,'color',[.4 .4 .4]);
            plot(self.time./1000,self.model.fit_tonic+self.model.fit_phasic,'r');
            grid on;
            axis tight;
            hold off
        end
        function leastsquare_decomposition(self,event_i)
            %will model the data using FIR and FourierBasis set.
            if ~isempty(self.data)
                FIR              = self.FIR(event_i);
                P                = self.FourierBasis;
                self.model.betas = [FIR P]\self.data;
                self.model.fit   = [FIR P]*self.model.betas;
                self.model.r2    = corr2(self.model.fit,self.data);
                % estimated tonic and phasic responses.
                self.tonic       = P*self.model.betas(end-size(P,2)+1:end);
                self.phasic      = FIR*self.model.betas(1:end-size(P,2));
            else
                fprintf('data field is empty dickhead\n');
            end
        end
        function run_ledalab(self)
            %will run ledalab on all the data. Use cut method to restrict
            %the analysis to a single block phase or so.
            %             addpath('/Users/onat/Documents/Code/Matlab/ledalab/');
            %
            foldername            = regexprep(self.path_scr,'data.smr','ledalab');%storage of ledalab related files
            if exist(foldername) == 0;mkdir(foldername);end%create it if necessary.
            filename              = sprintf(['%s%sdata_' repmat('%s_',1,length(self.BlockNames))],foldername,filesep,self.BlockNames{:});%used by ledalab
            filename(end)         = [];
            filename_results      = [filename '_results.mat'];%produced by ledalab
            filename              = [filename '.mat'];
            fprintf('filename to ledalab is %s\n',filename);
            
            if exist(filename_results) == 0
                %% convert to data format that ledalab understands.
                data.conductance            = self.y;
                data.time                   = self.time/1000;%in seconds
                data.time                   = data.time - min(data.time(:));
                data.samplingrate           = self.sampling_rate;%Hz
                data.samplingperiod         = self.sampling_period/1000;%in s
                %% transform events to ledalab format

                %store
                c = 0;
                for ncond = 1:9
                    for nEvent = find(self.event(:,ncond))'
                        c                   = c+1;
                        data.event(c).time  = (nEvent*data.samplingperiod);%in ms
                        data.event(c).nid   = ncond;%add a constant to differentiate phases.
                        data.event(c).name  = self.event_name{ncond};
                    end
                end
                save(filename,'data');
                %%
                addpath('/home/onat/Documents/Code/Matlab/ledalab/');
                Ledalab({filename},self.ledalab_defaults{:});
            end
            leda         = load(filename_results);
            % before assigning tonic and phasic values upsample the data to
            % original sampling rate so that we can use normally the
            % plotting methods.
            self.phasic  = interp1(leda.data.time.data,leda.analysis.driver,(self.time-self.time(1))./1000);
            self.tonic   = interp1(leda.data.time.data,leda.analysis.tonicData,(self.time-self.time(1))./1000);
            self.ledalab = leda.analysis.split_driver;
        end
        function plot_tuning_ledalab(self,varargin)
            if nargin > 1
                conds = varargin{1};
            else
                
                conds = 1:(length(find(sum(self.event)~=0))-2);
            end
            %will return average SCR values for conditions CONDS
            %(optional).
            if ~isempty(self.ledalab)%if ledalab analysis is done
                %detect time window, based on averaged data we take [1.5 4]
                %seconds. Could be improved for single-subject variations
                i                    = (self.ledalab.x(:,1) >= 1.5)&(self.ledalab.x(:,1) <= 4);
                self.fear_tuning     = mean(self.ledalab.mean(i,:));%take out the average in that window
                self.fear_tuning = self.fear_tuning(:,conds);
            else%if the analysis not done yet,
                self.run_ledalab;%first do it
                self.plot_tuning_ledalab(conds);%and call yourself.
            end
            bar(self.fear_tuning)
        end
        function ML_fit(self)
            
            %run the ML estimation
            params = [];
            %             if length(self.BlockNames) == 1
            
            filter_string                            = self.BlockNames{1};
            event_id                                 = cellfun(@(x) ~isempty( regexp(x,filter_string)), self.event_name );%detect the relevant event column
            [onset_sample x]                         = find(self.event(:,event_id));%in samples
            ml.onsets                                   = self.time(onset_sample)./1000;%onset times in s
            ml.data                                  = self.phasic(:);%data we want to model
            ml.time                                  = self.time(:)./1000;
            ml.data(ml.time < (min(ml.onsets)-10))      = 0;
            ml.data(ml.time > (max(ml.onsets)+10))      = 0;
            
            ml.funlsq                             = @(params) (abs(ml.data - scr_model( ml.time , ml.onsets(:), params ))).^2;%squared deviations
            ml.params0                            = [rand(1,length(ml.onsets)) 2 6 2];%initial values.
            ml.LB                                 = [zeros(1,length(ml.onsets)) 0 0 0];%lower boundaries
            %                 options                            = optimset('algorithm',{'levenberg-marquardt',.01},'display','iter','MaxFunEvals',50000,'maxiter',50000,'tolx',10^-12,'tolfun',10^-12,'OutputFcn',@scr_optimizer_plot);
            ml.options                            = optimset('algorithm',{'levenberg-marquardt',.01},'display','iter','MaxFunEvals',2500,'maxiter',2500,'tolx',10^-12,'tolfun',10^-12);
            ml.params                             = lsqnonlin(ml.funlsq,ml.params0,ml.LB,[],ml.options);
            self.ml = ml;
            %             else
            %                 fprintf('Please cut me so that I have exactly one block\n');
            %             end
        end
    end
end
