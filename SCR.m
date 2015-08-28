classdef SCR
    
    properties (Constant = true, Hidden)
        default_run = 4;%at which run the scr data is located.
    end
    properties (Hidden)
        hdr
        markers
        block2phase
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
        event           = [];
        event_name      = {};
        event_plotting  = {};
        y_smooth        = [];
        y_tonic         = [];
        y_phasic        = [];
        y_diff          = [];
        y_diff2         = [];
        y_tonic_spline  = [];
    end
    
    methods
        %add a method for tonic removal a la ledalab and a la new baseline
        %correction
        %add a method to compute FIR analysis and plot the results based on
        %the corrected phasic responses.
        function scr = SCR(s)
            path_acqfile  = s.path2data(SCR.default_run,'scr','acq');
            if exist(path_acqfile);
                dummy         = load_acq(path_acqfile);
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
                i(diff(scr.time(i)/1000) < 30)            = [];
                scr.BlockBorders_index = [i [i(2:end); scr.tsample]];
                scr.BlockBorders_time  = scr.time(scr.BlockBorders_index);
                tblock                 = length(i);
                % give names to blocks (end to start)
                block_names          = {'test_rating' 'test' 'cond_rating' 'cond' 'base_rating' 'base' };
                %append UCS calibration blocks as much as necessary
                while length(block_names) ~= tblock
                    block_names      = [block_names {sprintf('calibration%02d',tblock - length(block_names))}];
                end
                block_names = fliplr(block_names);%correct for the order.
                %% label each SCR block in terms of experimental phase
                scr.block2phase = nan(1,tblock);
                scr.block2phase([(tblock-5) (tblock-3) (tblock-1)]) = [2 3 4];                
                %%
                for nblock = find(~isnan(scr.block2phase))
                    % event times in units of samples                                        
                    i                      = scr.BlockBorders_index(nblock,1):scr.BlockBorders_index(nblock,2);%range of the block
                    onset_sample           = find(data(:,3));%all stim onsets
                    hold on
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
                        scr.event_name                 = [ scr.event_name         sprintf('%s_%d',block_names{nblock},cond_ids(ii))];
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
        end
    end
    methods
        
        function plot(self)
            %% plot different event channels
            
            figure(1);clf;
            set(gcf,'position',[1         444        1440         362]);
            % plot the SCR
            plot(self.time,self.y,'color','k');
            hold on;
            % mark stimulus onsets
            for n = 1:length(self.event_name)                
                i  = self.event(:,n);                 
                plot(self.time(i),self.y(i), self.event_plotting{n}.symbol{1}{2} , self.event_plotting{n}.color{1}{:} , self.event_plotting{n}.marker_size{1}{:} );
            end            
            %% phase transitions as lines
            i   = self.BlockBorders_time(:,1);%start of phases
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
        end        
        function plot2(self)    
            %plots tonic, phasic etc etc.
            self.plot
            hold on
%             plot(self.time,self.y_smooth,'r');
            plot(self.time,self.y_phasic,'k');            
            plot(self.time,self.y_diff,'color',[.4 .4 .4]);
            plot(self.time,self.y_diff2,'color','r');
            % mark stimulus onsets
            for n = 1:length(self.event_name)                
                i  = self.event(:,n);                 
                plot(self.time(i),0, self.event_plotting{n}.symbol{1}{2} , self.event_plotting{n}.color{1}{:} , self.event_plotting{n}.marker_size{1}{:} );
            end            
            grid on;    
            axis tight;            
        end
        function o = get.y_tonic_spline(self)
            keyboard    
            Y     = self.y(self.BlockBorders_index(1,1):self.BlockBorders_index(1,2));
            X     = self.time(self.BlockBorders_index(1,1):self.BlockBorders_index(1,2));
            [~,m] = DetectPeaks(Y);
            s     = std(m);
            i     = find(s > mean(std(rand(3000)+1)));
            peaks = find(diff(i) > 20);
            
        end
        function ys = get.y_smooth(self)
            %smoothing has to be done with splines.
            ys      = smooth(self.y,750,'lowess');
            self.y_smooth = ys;
        end
        function ytonic = get.y_tonic(self)
            d       = fdesign.lowpass('Fp',0.00001);
            Hd      = design(d, 'ellip');
            ytonic = filtfilt(Hd.sosMatrix,Hd.ScaleValues,self.y_smooth);
            self.y_tonic = ytonic;
        end
        function yphasic = get.y_phasic(self)
            yphasic        = self.y_smooth - self.y_tonic;
            yphasic        = diff(yphasic);
            yphasic(end+1) = yphasic(end);
            yphasic        = zscore(yphasic);
            self.y_phasic  = yphasic;
        end        
        function ydiff = get.y_diff(self)            
            ydiff        = diff(self.y_phasic);
            ydiff(end+1) = ydiff(end);
            ydiff        = zscore(ydiff);            
            self.y_diff = ydiff;
        end
        function ydiff = get.y_diff2(self)            
            ydiff              = diff(diff(self.y_phasic));
            ydiff(end+1:end+2) = ydiff(end);
            %
            m = mean(ydiff);
            s = median(abs(ydiff - m));
            ydiff = (ydiff-m)*100000;
            self.y_diff2        = ydiff;
            ylim([-10 50]);
        end
        function xcorr(self,block)                 
            %%
            figure(2);clf;
            cc = 0;
            for conds = {1:11 12:16 17:27 28:29};
                cc = cc +1;
                subplot(1,4,cc)
                for n = conds{1}                    
                    [r lag] = xcorr(self.y_diff2,double(self.event(:,n)),3000);                                        
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
        
    end
    
end