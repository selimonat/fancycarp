classdef SCR < handle
    
    properties (Constant = true, Hidden)
        default_run = 4;%at which run the scr data is located.
    end
    properties (Hidden)
        ledalab_defaults      = {'open', 'mat','downsample', 5, 'analyze','CDA', 'optimize',10, 'overview',  1, 'export_era', [-1 7 0 1], 'export_scrlist', [0 1], 'export_eta', 1 };%
        hdr
        markers
        block2phase
        p         = 0.0001;
        nfreq     = 40;
        FIR_delay = 15000;%in milliseconds
        path_acqfile;%path to the file
        ledalab_all
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
                s = varargin{1};
                scr.path_acqfile  = s.path2data(s.id,SCR.default_run,'scr','acq');
                if exist(scr.path_acqfile);
                    if ~exist(sprintf('%s.mat',scr.path_acqfile))
                        dummy     = load_acq(scr.path_acqfile);
                        save(sprintf('%s.mat',scr.path_acqfile),'dummy')
                    else
                        load(sprintf('%s.mat',scr.path_acqfile));%will spawn dummy
                    end
                    scr.hdr       = dummy.hdr;
                    data          = dummy.data;
                    data          = data(1:end-10,:);%remove the ends, coz of big jumps
                    scr.y         = data(:,1);
                    scr.markers   = dummy.markers;
                    scr.tsample   = length(data);
                    %% add some goodies
                    chan_names        = {scr.hdr.per_chan_data(:).comment_text};
                    %             'SCR'    'FixOnset'    'StimOnset'    'ShockOnset'    'Oddball'    'InitPhase'    'Keypress'
                    scr.sampling_rate     = 1000/scr.hdr.graph.sample_time;%Hz
                    scr.sampling_period   = scr.hdr.graph.sample_time;%ms
                    scr.time              = [0:scr.sampling_period:(scr.tsample-1)*scr.sampling_period]';%ms
                    %% add also the start of sub-phases.
                    %find the channel where phase starts are stored
                    c                     = find(strcmp(chan_names,'InitPhase'));
                    %and detect the time index
                    i                      = find(diff(data(:,c)) > 0);
                    %most probably 2 pulses faster than 30s is a falschbeginn
                    i(diff(scr.time(i)/1000) < 40)            = [];
                    %%
                    scr.BlockBorders_index = [i [i(2:end); scr.tsample]];
                    scr.BlockBorders_time  = scr.time(scr.BlockBorders_index);
                    tblock                 = length(i);
                    % give names to blocks (end to start)
                    block_names            = s.scr_blocknames;%inherited from the Project
                    %append UCS calibration blocks as much as necessary
                    while length(block_names) ~= tblock
                        block_names      = [block_names {sprintf('calibration%02d',tblock - length(block_names))}];
                    end
                    block_names    = fliplr(block_names);%correct for the order.
                    %% label each SCR block in terms of experimental phase
                    scr.block2phase                                     = nan(1,tblock);
                    scr.block2phase([(tblock-5) (tblock-3) (tblock-1)]) = [2 3 4];
                    scr.BlockNames                                      = block_names;
                    %%
                    for nblock = find(~isnan(scr.block2phase))
                        % event times in units of samples
                        i                      = scr.BlockBorders_index(nblock,1):scr.BlockBorders_index(nblock,2);%range of the block
                        onset_sample           = find(data(:,3));%all stim onsets
                        onset_sample           = onset_sample(ismember(onset_sample,i));%take only this blocks's onsets
                        onset_sample           = onset_sample(diff([-Inf; onset_sample]) > 1);%clean repeats
                        onset_sample           = onset_sample(ismember(onset_sample,i));%take only this phase's onsets
                        % get stimulus indices for the above onsets
                        current_phase          = scr.block2phase(nblock);
                        cond_seq               = s.paradigm{current_phase}.presentation.dist;
                        if length(cond_seq) ~= length(onset_sample);keyboard;end;%if problem then stoPPP
                        cond_ids               = unique(cond_seq);
                        t_stim_id              = length(cond_ids);%total stim id
                        % run through different stimulus indices
                        zero_vec               = false(scr.tsample,t_stim_id);
                        for ii = 1:t_stim_id
                            ind                            = ismember(cond_seq,cond_ids(ii));%samples for this stim_id
                            zero_vec(onset_sample(ind),ii) = true;%
                            %replace stimall with face index
                            scr.event_name                 = [ scr.event_name         sprintf('%s_%04d',block_names{nblock},cond_ids(ii))];
                            scr.event_plotting             = [ scr.event_plotting     s.plot_style(cond_ids(ii))];
                        end
                        scr.event                          = logical([scr.event zero_vec]);
                    end
                    %% create event channels also for the UCS during the calibration
                    %these are not in the paradigm file.
                    for nblock = 1:tblock-6%run through the calibration blocks
                        zero_vec               = false(scr.tsample,1);
                        %samples where this phase recorded
                        i                      = scr.BlockBorders_index(nblock,1) : scr.BlockBorders_index(nblock,2);
                        onsets                 = data(:,4);%UCS onsets
                        onset_sample           = find(onsets);
                        onset_sample           = onset_sample(diff([-Inf; onset_sample]) > 1);%clean repeats
                        onset_sample           = onset_sample(ismember(onset_sample,i));%take only this phase's onsets
                        zero_vec(onset_sample) = true;
                        scr.event              = logical([scr.event logical(zero_vec)]);
                        scr.event_name         = [scr.event_name         sprintf('%s_%s',block_names{nblock},'1001')];
                        scr.event_plotting     = [scr.event_plotting            {s.plot_style(cond_ids(ii))}];
                    end
                    %% make a channel based on all the UCSs
                    zero_vec               = false(scr.tsample,1);
                    onsets                 = data(:,4);
                    onset_sample           = find(onsets);
                    onset_sample           = onset_sample(diff([-Inf; onset_sample]) > 1);%clean repeats
                    zero_vec(onset_sample) = true;
                    scr.event              = logical([scr.event              zero_vec]);
                    scr.event_name         = [ scr.event_name         sprintf('1001')];
                    scr.event_plotting     = [ scr.event_plotting            {s.plot_style(1001)}];
                    %% make a channel based on all the ODDs
                    zero_vec               = false(scr.tsample,1);
                    onsets                 = data(:,5);
                    onset_sample           = find(onsets);
                    onset_sample           = onset_sample(diff([-Inf; onset_sample]) > 1);%clean repeats
                    zero_vec(onset_sample) = true;
                    scr.event              = logical([scr.event              zero_vec]);
                    scr.event_name         = [scr.event_name         sprintf('1002')];
                    scr.event_plotting     = [scr.event_plotting            {s.plot_style(1002)}];
                else
                    fprintf('No acq file found, :(.\n');
                end
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
    methods %decomposition methods
        function xcorr(self,block)
            %%
%             self.data = [];
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
    end
    methods        
        function self = tonic_lowpass(self)
            %computes the low-pass version of the signal supposed to
            %represent the tonic response.
            if isempty(self.data);self.data=self.y;end
            d            = fdesign.lowpass('Fp',0.00001);
            Hd           = design(d, 'ellip');
            self.tonic   = filtfilt(Hd.sosMatrix,Hd.ScaleValues,self.data);
            self.phasic  = self.data - self.tonic;
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
        function run_ledalab(self)
            %will run ledalab on all the data. Use cut method to restrict
            %the analysis to a single block phase or so.
            %
%             addpath('/Users/onat/Documents/Code/Matlab/ledalab/');
            %
            foldername            = regexprep(self.path_acqfile,'data.acq','ledalab');%storage of ledalab related files
            if exist(foldername) == 0;mkdir(foldername);end%create it if necessary.
            filename              = sprintf(['%s%sdata_' repmat('%s_',1,length(self.BlockNames))],foldername,filesep,self.BlockNames{:});%used by ledalab
            filename(end)         = [];
            filename_results      = [filename '_results.mat'];%produced by ledalab
            filename              = [filename '.mat'];
            fprintf('filename to ledalab is %s\n',filename);
            
            if exist(filename_results) == 0
                %% convert to data format that ledalab understands.
                self.smooth('sgolay');
                data.conductance            = self.data;
                data.time                   = self.time/1000;%in seconds
                data.time                   = data.time - min(data.time(:));
                data.samplingrate           = self.sampling_rate;%Hz
                data.samplingperiod         = self.sampling_period/1000;%in s
                %% transform events to ledalab format
                %find all conditions that are in this batch
                conditions = [];
                for bnames = self.BlockNames;
                    conditions                  = [conditions find(cellfun(@(x) ~isempty(regexp(x,bnames{1})), self.event_name ))];%detect only the required conditions
                end
                %store
                c = 0;
                for ncond = conditions
                    for nEvent = find(self.event(:,ncond))'
                        c                   = c+1;
                        data.event(c).time  = (nEvent*data.samplingperiod);%in ms
                        data.event(c).nid   = ncond;%add a constant to differentiate phases.
                        data.event(c).name  = self.event_name{ncond};
                    end
                end
                save(filename,'data');
                %%
                Ledalab({filename},self.ledalab_defaults{:});
           end
            leda         = load(filename_results);
            % before assigning tonic and phasic values upsample the data to
            % original sampling rate so that we can use normally the
            % plotting methods.
            self.phasic  = interp1(leda.data.time.data,leda.analysis.driver,(self.time-self.time(1))./1000);
            self.tonic   = interp1(leda.data.time.data,leda.analysis.tonicData,(self.time-self.time(1))./1000);
            self.ledalab = leda.analysis.split_driver;
            self.ledalab_all=leda;
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
    methods %utilities
        function out              = cut(self,block)
            out = SCR(self,block);
        end       
        function out              = findphase(self,filter_string);
            out = cellfun(@(x) ~isempty(regexp(x,filter_string)), self.BlockNames );
            out = find(out);
        end
        function linespecs        = LineSpecs(condition)
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
        function self             = diff(self)
            %takes the temporal difference.
            if isempty(self.data);self.data=self.y;end
            self.data         = diff(self.data);
            self.data(end+1)  = self.data(end);
        end       
        function self             = smooth(self,method)
            %
            if isempty(self.data);self.data=self.y;end
            
            if strcmp(method,'ma');%moving average
                self.data  = smooth(self.data,750);
            elseif strcmp(method,'sgolay');%moving average
                self.data  = sgolayfilt(self.data,3,751);
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
        function [P]              = FourierBasis(self,tsample)
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
        function [out_z out_raw] = getTuning(self,timeframe)
            out_raw = nan(3*8,1);
            condcollector = { 'base_0045','base_0090','base_0135','base_0180','base_0225','base_0270','base_0315','base_0360',...
                'cond_0180','cond_0360',...
                'test_0045','test_0090','test_0135','test_0180','test_0225','test_0270','test_0315','test_0360'};
            index  = [1:8 12 16 17:24];
            for c = 1:length(condcollector)
                timey = (self.ledalab.x(:,1) >= min(timeframe))&(self.ledalab.x(:,1) <= max(timeframe));%time window
                condy = strcmp(condcollector{c},self.ledalab.condnames)';
                dummy = mean(self.ledalab.y(timey,condy));
                out_raw(index(c),1)       = mean(dummy);
                %out_rawsd(index(c),1)     = std(dummy)./sqrt(length(dummy));
            end
            out_z = nanzscore(out_raw);
        end
    end
    methods %plotters
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
                i                    = (self.ledalab.x(:,1) >= 2)&(self.ledalab.x(:,1) <= 5.0);%time window
                M = [];S=[];
                for c = conds
                    dummy = mean(self.ledalab.y(i,self.ledalab.c == c));
                    M(c)     = mean(dummy);
                    S(c)     = std(dummy)./sqrt(length(dummy));
                end
                M                    = M(conds);
                self.fear_tuning     = M;%take out the average in that window
                S                    = S(conds);
            else%if the analysis not done yet,
                self.run_ledalab;%first do it
                self.plot_tuning_ledalab(conds);%and call yourself.
            end
            h = bar(self.fear_tuning);
            try
                SetFearGenBarColors(h);
            end
            hold on;
            errorbar(M,S,'ok')
            hold off;
        end
    end
end