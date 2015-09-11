classdef Fixmat < Project
    properties (Hidden,Constant)
        kernel_fwhm = Fixmat.PixelPerDegree*1.75;
        dummy       = make_gaussian2D(Fixmat.kernel_fwhm*2.65,Fixmat.kernel_fwhm*2.65,Fixmat.kernel_fwhm,Fixmat.kernel_fwhm,Fixmat.kernel_fwhm*2.65/2,Fixmat.kernel_fwhm*2.65/2);
        kernel      = Fixmat.dummy./sum(Fixmat.dummy(:));       
        %         kernel = Fixmat.dummy;
    end
    
    properties (Hidden,SetAccess = public)
         baseline_correction = 0;
         unitize             = 0;
         mapsize             = 250;
         maptype             = 'bin';%conv or bin
         binfactor           = 25;
    end
    properties (Hidden,SetAccess = private)
        maps
        query
        map_titles
        selection
    end
    
    properties (SetAccess = private)
        subject = [];
        phase   = [];
        start   = [];
        stop    = [];
        x       = [];
        y       = [];
        eye     = [];
        deltacsp= [];
        file    = [];
        fixx    = [];
        fixy    = [];
        oddball = [];
        trialid = [];
        ucs     = [];
        fix     = [];
        rect    = [];
    end
    
    methods
        function obj = Fixmat(subjects,runs)%constructor
            %%
            %initialize
            for subject = subjects
                for run = runs
                    %%                    
                    dummy = load(obj.path2data(subject,run,'eye'));
                    %this is necessary to expand the PTB message to
                    %something that is understandable by the fixations
                    %method.
                    for nt = 1:length(dummy.trials)
                        %extract relevant info first from the text message
                        T = dummy.trials(nt).TRIALID.msg(1,:);
                        %add them as msg so that fixations.m can create the fields
                        pairs = regexp(T,'[\w-]*','match');
                        for np = 1:2:length(pairs)-1
                            dummy.trials(nt).(lower(pairs{np})).msg = pairs{np+1};
                        end
                    end
                    %
                    dummy.trials = rmfield(dummy.trials,{'TRIALID'});%make fixations happy
                    dummy        = obj.getfixmat(dummy.trials,dummy.info);%get the fixmat
                    dummy.subject= uint32(repmat(subject,1,length(dummy.x)));
                    %and append it to the previous fixmat
                    for fns = fieldnames(dummy)'
                        if isempty(obj.(fns{1}))
                            obj.(fns{1}) = [];
                        end
                        obj.(fns{1}) = [obj.(fns{1}) dummy.(fns{1})];
                    end
                end
                obj.x = round(obj.x);
                obj.y = round(obj.y);
            end
            %% add the rect
            obj.rect = [1 obj.screen_resolution(1) 1 obj.screen_resolution(2)];
            %% take fixations which are only coming from the required phase.            
            obj.UpdateSelection('phase',runs);
            obj.ApplySelection;
            %% remove fixations outside of the image border
            obj.selection = ~(obj.x <= obj.rect(3) | obj.x > obj.rect(4) | obj.y <= obj.rect(1) | obj.y > obj.rect(2));
            obj.ApplySelection;
        end
        function UpdateSelection(obj,varargin)
            %takes VARARGIN pairs to update the selection vektor            
            obj.selection = logical(ones(1,length(obj.subject)));
            for n = 1:2:length(varargin)
                obj.selection = obj.selection .* ismember(obj.(varargin{n}),varargin{n+1});
            end                        
            fprintf('Selection vector updated...\n');
            obj.selection = logical(obj.selection);
            obj.query     = varargin;
        end
        
        function ApplySelection(obj)
            %removes fixations based on the fixation vecktor, only used in
            %the constructor.            
            %removes fixations which are FALSE in selection            
            for p = properties(obj)'
                if ~strcmp(p{1},'rect') && ~strcmp(p{1},'maps') && ~strcmp(p{1},'selection')
                    obj.(p{1})(~obj.selection) = [];
                end
            end
            obj.query = [];
            fprintf('Selection removed from the object...\n');
        end
        function plot(obj)            
            [d u] = GetColorMapLimits(obj.maps,9);            
            if ~obj.baseline_correction
                d = 0;
            end                                
            %
            tmaps = size(obj.maps,3);
            hhfigure;
            for nc = 1:size(obj.maps,3)
                h     = subplot(2,4,nc);
                %plot the image;
                bild  = repmat(obj.cropimage(imread(obj.find_stim)),[1 1 3]);
                x     = [1:size(obj.maps,1)]';
                y     = [1:size(obj.maps,2)]';
                h     = imagesc(x,y,bild);
                hold on;
                %                 subplotChangeSize(h,0.05,0.05);
                h     = imagesc(x,y,obj.maps(:,:,nc),[d u]);
                if ~obj.baseline_correction
                    set(h,'alphaData',Scale(obj.maps(:,:,nc)));
                else
                    set(h,'alphaData',Scale(abs(obj.maps(:,:,nc))));
                end
                axis image;
                axis off;                
                t     = sprintf('%s%d/',obj.map_titles{nc}{:});
                title(t,'interpreter','none');
            end
        end
        function getmaps(obj,varargin)
            %will populate the maps property based on the filter in
            %VARARGIN, example:
            %fs.getmaps({'phase',4,'deltacsp',-135});
            
            obj.maps = [];%clear whatever is there
            c  = 0;
            for v = varargin
                c                 = c+1;

                obj.UpdateSelection(v{1}{:});                
                
                if strcmp(obj.maptype,'conv')
                    %accum and conv
                    FixMap            = accumarray([obj.current_y' obj.current_x'],1,[obj.rect(2) obj.rect(4)]);
                    FixMap            = conv2(sum(obj.kernel),sum(obj.kernel,2),FixMap,'same');
                elseif strcmp(obj.maptype,'bin')
                    %divide by a factor (i.e. binning) and accum, but no conv
                    y                 = ceil(double(obj.current_y')./obj.binfactor);
                    x                 = ceil(double(obj.current_x')./obj.binfactor);
                    FixMap            = accumarray([y x],1,[obj.rect(2)./obj.binfactor obj.rect(4)./obj.binfactor]);
                end
                %
                FixMap                = obj.cropmaps(FixMap);
                if obj.unitize
                    FixMap            = FixMap./sum(FixMap(:));
                end
                obj.maps(:,:,c)       = FixMap;
                %
                obj.map_titles{c} = obj.query;
            end
            
            if obj.baseline_correction
                obj.maps = obj.maps - repmat(mean(obj.maps,3),[1 1 size(obj.maps,3)]);
            end            
        end
      
        function out = cropimage(obj,im)            
            out      = im(obj.rect(2)/2-obj.mapsize : obj.rect(2)/2+obj.mapsize, obj.rect(4)/2-obj.mapsize: obj.rect(4)/2+obj.mapsize);
        end
        
        function out = cropmaps(obj,map)
            if strcmp(obj.maptype,'conv')
                out      = map(obj.rect(2)/2-obj.mapsize : obj.rect(2)/2+obj.mapsize, obj.rect(4)/2-obj.mapsize: obj.rect(4)/2+obj.mapsize);
            elseif strcmp(obj.maptype,'bin')
                new_size = obj.mapsize./obj.binfactor;
                out      = map(obj.rect(2)/obj.binfactor/2-new_size : obj.rect(2)/obj.binfactor/2+new_size, obj.rect(4)/obj.binfactor/2-new_size: obj.rect(4)/obj.binfactor/2+new_size);
            end
        end
        
        function maps   = vectorize_maps(obj)
            if ~isempty(obj.maps)
                maps = reshape(obj.maps,size(obj.maps,1)*size(obj.maps,2),size(obj.maps,3));
            else
                fprintf('no maps stored here\n');
            end
        end
        function cov(obj)
            set(gcf,'position',[1104         152         337         299]);
            imagesc(cov(obj.vectorize_maps));
            axis image;
            colorbar
        end
        function corr(obj)
            set(gcf,'position',[1104         152         337         299])
            imagesc(corr(obj.vectorize_maps));
            axis image;
            colorbar
        end
        function [x] = current_x(obj)
            %returns the current x y coordinates            
            x = obj.x(obj.selection);
        end
        function [y] = current_y(obj)
            %returns the current x y coordinates
            y = obj.y(obj.selection);            
        end
        function histogram(obj)
            subject   = unique(obj.subject);
            phase     = unique(obj.phase);
            condition  = unique(obj.deltacsp);
            count = nan(length(subject),length(condition));
            sub = 0;
            for ns = subject
                sub = sub + 1;
                ph = 0;
                for np = phase
                    ph = ph + 1;
                    cond = 0;
                    for nc = condition
                        cond = cond + 1;
                        dummy = obj.FilterFixmat('subject',ns,'phase',np,'deltacsp',nc);
                        count(sub,cond,ph) = length(dummy.x);
                    end
                end
            end
            figure;
            for phases = 1:ph
                subplot(1,3,phases)
                imagesc(count(:,:,phases));
            end
        end
        function fixmat = getfixmat(obj,edfdata, edfmeta)
            % fixmat = fixations(edfdata, edfmeta)
            %   where [edfdata, edfmeta] = edfread('file.edf')
            %   returns a flat fixation matrix
            %   returns user defined meta data per trial and per experiment
            %   removes pre-stimulus-onset fixations
            
            % parameter with_pupil: return pupil data
            with_pupil = 0;
            % number of trials
            trials = length(edfdata);
            
            % get field names of custom per-trial meta data
            function my_fields = trial_metadata()
                fields = fieldnames(edfdata);
                my_fields = setdiff(fields, {'left', 'right', 'button'});
            end
            
            % get field names of custom per-experiment meta data
            function my_fields = experiment_metadata()
                fields = fieldnames(edfmeta);
                my_fields = setdiff(fields, {'calib', 'header'});
            end
            
            % get empty trials for each eye
            function empty = empty_trials()
                empty = zeros(2, trials);
                for i=(1:trials)
                    empty(1, i) = isstruct(edfdata(i).left);
                    empty(2, i) = isstruct(edfdata(i).right);
                end
            end
            
            % get eye with best calibration
            function better = eye_calib()
                % init result:
                better = zeros(2, trials);
                % number of calibrations realized
                ncalib = size(edfmeta.calib, 1);
                % running from the back, in which trial we calibrated the last time
                lastc  = size(edfdata, 2) + 1;
                % run over calibrations, from the last one to the first one
                for c=(ncalib:-1:1)
                    % get current calibration
                    cal = edfmeta.calib(c);
                    % is the calibration more recent than the last one?
                    if cal.trial ~= lastc
                        % we have a fresh calibration
                        if isstruct(cal.left) && ((~isstruct(cal.right) || cal.left.err_avg <  cal.right.err_avg) || isnan(cal.right.err_avg)) && ( ~ isnan(cal.left.err_avg) )
                            % left eye has better calibration
                            better(1, max(1, cal.trial):lastc-1) = 0.5;
                            lastc = cal.trial;
                        elseif ~ isnan(cal.right.err_avg)
                            % right eye rocks
                            better(2, max(1, cal.trial):lastc-1) = 0.5;
                            lastc = cal.trial;
                        end % no suitable calibration. we will skip this
                    end
                end
                if lastc~=0
                    display(['Error: no calibration before trial ', int2str(lastc), ' - using left eye']);
                end
            end
            
            % initialize the output structure
            function fixmat = init_struct(fields)
                fixmat.start = zeros(0, 0, 'int32');
                fixmat.stop  = zeros(0, 0, 'int32');
                fixmat.x     = zeros(0, 0, 'single');
                fixmat.y     = zeros(0, 0, 'single');
                fixmat.fix   = zeros(0, 0, 'int32');
                if with_pupil
                    fixmat.pupil = zeros(0, 0, 'single');
                end
                fixmat.eye   = zeros(0, 0, 'uint8');
                for i=(1:length(fields))
                    fixmat.(fields{i}) = zeros(0, 0, 'int32');
                end
            end
            
            % determine which eye to use
            eye = eye_calib() + empty_trials();
            want_eye = 1+(eye(1,:)<eye(2,:));
            eyes = {'left'; 'right'};
            
            % read user defined meta data
            trial_fields = trial_metadata();
            trial_fields_length = length(trial_fields);
            exp_fields = experiment_metadata();
            exp_fields_length = length(exp_fields);
            
            % prepare output structure
            fixmat = init_struct(union(trial_fields, exp_fields));
            
            % main loop: iterate over single trials
            for trial=(1:trials)
                % get fixations of 'better' eye
                fix = edfdata(trial).(eyes{want_eye(trial)}).fixation;
                tstart = fix.start;
                % skip pre-stimulus fixations
                idx = find(tstart>0);
                fix_in_trial = length(idx);
                % output: eye (1 - left, 2 - right)
                fixmat.eye = [fixmat.eye, ones(1, fix_in_trial) * want_eye(trial)];
                % add user defined meta data
                for f=(1:trial_fields_length)
                    % check for empty fields, mark them with -maxint
                    if isempty(edfdata(trial).(trial_fields{f}))
                        edfdata(trial).(trial_fields{f}).msg = '-2147483648';
                    end
                    field_val = str2double(edfdata(trial).(trial_fields{f}).msg);
                    fixmat.(trial_fields{f}) = [fixmat.(trial_fields{f}), ones(1, fix_in_trial) * field_val];
                    
                end
                % add default data fields
                fixmat.start = [fixmat.start, tstart(idx)];
                fixmat.stop  = [fixmat.stop, fix.end(idx)];
                fixmat.x     = [fixmat.x, fix.x(idx)];
                fixmat.y     = [fixmat.y, fix.y(idx)];
                fixmat.fix   = [fixmat.fix,1:fix_in_trial];
                if with_pupil
                    fixmat.pupil = [fixmat.pupil, fix.pupil(idx)];
                end
            end;
            % add user defined meta data
            for f=(1:exp_fields_length)
                fixmat.(exp_fields{f}) = ones(1, length(fixmat.start), 'int32') * str2double(edfmeta.(exp_fields{f}));
            end
            
        end
        
        
    end
end