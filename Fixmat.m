classdef Fixmat < Project
    properties (Hidden,Constant)                
        window       = 250;%half size of the fixation maps        
    end
    
    properties (Hidden,SetAccess = public)
        %related to map computation
         bc                  = 0;%cocktail blank correction
         unitize             = 1;%sums to 1 or not         
         maptype             = 'conv';%conv or bin         
         kernel_fwhm         = Fixmat.PixelPerDegree*.8;
         binsize             = 25;
         linkage_method      = 'average';
         linkage_metric      = 'correlation';
         maps%current maps;
    end
   
    properties (Hidden,SetAccess = private)
        %internal
        query
        all_queries
        map_titles
        selection
        realcond%all conditions that are not ucs,odd,or nulltrial
    end
    
    properties (SetAccess = private,Dependent,Hidden)
        %rect: [y x y_size x_size]
        rect       %the aperture
        binedges   %bin edges for binned fixation maps        
        stimulus   %contains average stimulus              
    end
    
    properties (SetAccess = private)
        %experimental observations
        subject = [];
        phase   = [];
        start   = [];
        stop    = [];
        x       = [];
        y       = [];
        eye     = [];
        deltacsp= [];
        file    = [];
        oddball = [];
        trialid = [];
        ucs     = [];
        fix     = [];
        chain   = [];
        isref   = [];
    end
    
    events
      UpdateGraph 
    end
    
    methods        
        function obj = Fixmat(subjects,runs)%constructor
            %%
            %initialize
            for run = runs(:)'
                for subject = subjects(:)'                
                    if exist(regexprep(obj.path2data(subject,run,'eye'),'.mat','.edf')) ~= 0 %does the subject has an edf file?

					%% edf2mat conversion                    
					if ~exist(obj.path2data(subject,run,'eye'))%is the edf file converted?
							fprintf('-edfread not yet ran (subject:%03d, run:%03d)\n-Roger that, will do it ASAP, sir!\n',subject,run);
							[trials info] = edfread(regexprep(obj.path2data(subject,run,'eye'),'.mat','.edf'),'TRIALID');
							save(obj.path2data(subject,run,'eye'),'trials','info');
							%if the conversion fails the code will fail, it is recommended to rename the data.edf file to something else then.
					else
						fprintf('edf is converted (subject:%03d, run:%03d)\n',subject,run);
					end

                    dummy = load(obj.path2data(subject,run,'eye'));
                    %this is necessary to expand the PTB message to
                    %something that is understandable by the fixations
                    %method.
                    for nt = 1:length(dummy.trials)
                        %extract relevant info first from the text message
                        T = dummy.trials(nt).TRIALID.msg(1,:);
                        %add them as msg so that .getfixmat can create the fields
                        pairs = regexp(T,'[\w-]*','match');
                        for np = 1:2:length(pairs)-1
                            dummy.trials(nt).(lower(pairs{np})).msg = pairs{np+1};
                        end
                    end
                    %
                    dummy.trials = rmfield(dummy.trials,{'TRIALID'});%make .getfixmat happy
                    dummy        = obj.getfixmat(dummy.trials,dummy.info);%get the fixmat                    
                    dummy.subject= repmat(uint32(subject),1,length(dummy.x));
                    %and append it to the previous fixmat
                    for fns = properties(obj)'
                        %% take fixations which are only coming from the required phase.            
                        %(e.g. ratings in baseline are coded as 5)                                    
                        valid_fix = dummy.phase == run;
                        if isfield(dummy,fns{1})%if it is not a property dont even consider                            
                            obj.(fns{1}) = [obj.(fns{1}) dummy.(fns{1})(valid_fix)];%append it to the previous
                        else
                            obj.(fns{1}) = [obj.(fns{1}) zeros(1,sum(valid_fix))];%append it to the previous
                        end
                    end
				end
                end                    
            end                                                    
            %% remove fixations outside of the image border
            obj.selection = ~(obj.x < obj.rect(2) | obj.x > (obj.rect(2)+obj.rect(4)-1) | obj.y < obj.rect(1) | obj.y > (obj.rect(1)+obj.rect(3)-1) );
            obj.ApplySelection;
            %% round coordinates to pixels
            obj.x = round(obj.x);
            obj.y = round(obj.y);
            obj.realcond = unique(obj.deltacsp(~ismember(obj.deltacsp,[500 1000 3000])));
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
                if ~strcmp(p{1},'rect') && ~strcmp(p{1},'maps') && ~strcmp(p{1},'selection');%all the properties unrelated to fixation data.
                    if ~isempty(obj.(p{1}))%sometimes the property is not filled in
                        obj.(p{1})(~obj.selection) = [];
                    end
                end
            end
            obj.query = [];
            fprintf('Selection (%04d fixations) removed from the object...\n',sum(~obj.selection));
        end
        function [cormat] = condcorrmat(obj,ph,plot)
            tsub = length(unique(obj.subject));
            cormat = nan(8,8,tsub);
            subc            = 0;
            for sub = unique(obj.subject);
                subc = subc + 1;
                %create the query cell
                v = [];
                c = 0;
                for cond = -135:45:180
                    c    = c+1;
                    v{c} = {'phase', ph, 'deltacsp' cond 'subject' sub};
                end
                % plot and save fixation maps
                obj.getmaps(v{:});
                maps = obj.maps;
                obj.maps   = maps(:,:,1:8) - repmat(mean(maps(:,:,1:8),3),[1 1 8]);%take out the average
                cormat(:,:,subc) = obj.corr;
            end
            if plot ==1
                figure;
                imagesc(median(cormat,3),[-.1 .1])
                axis square;colorbar
                set(gca,'fontsize',15)
                axis off
            end
        end
        function getsubmaps(obj)
            v = [];
            c = 0;
            for sub=unique(obj.subject)
                c    = c+1;
                v{c} = {'subject' sub 'deltacsp' obj.realcond};
            end
            obj.getmaps(v{:});
        end
        function [H,T,order,tree] = dendrogram(obj)%varargin tells if reorder by optimal leafOrder
            vecmaps = obj.vectorize_maps;
            tree = linkage(vecmaps',obj.linkage_method,obj.linkage_metric);
            D = pdist(vecmaps',obj.linkage_metric);
            leafOrder = optimalleaforder(tree,D);
            figure;
            [H,T,order] = dendrogram(tree,0,'Reorder',leafOrder);
            title('optimal leaforder')
        end
        function plotband(obj,varargin)%varargin can reorder the subjmaps
            if nargin > 1
                order = varargin{1};
            else
                order = unique(obj.subject);
            end
            N = length(order);
            figure;
            c=0;
            for sub = order;
                c=c+1;
                obj.getmaps({'subject' sub 'deltacsp' obj.realcond})
                h = subplot(1,N,c);imagesc(obj.maps);
%                 title(num2str(sub))
                axis square
                axis off
                subplotChangeSize(h,.01,.01);
                colorbar off
            end
        end
        function plot(obj,varargin)    
            
            M = obj.maps;
            if nargin > 1
                cmap = varargin{1};
            else
                cmap = 'linear';
            end
            %get colormap limits.
            if strcmp(cmap,'linear')
                [d u] = GetColorMapLimits(M,7);
                if sum(obj.maps(:) < 0) == 0%if there are no negative values
                    d = 0;
                end                                      
            elseif strcmp(cmap,'log')                            
                u = max(log10(M(:)));
                d = u - 1;
                M = log10(M);            
            end
            %
            ffigure;
            clf                        
            nsp     = obj.subplot_number;
            for nc = 1:size(M,3)
                h   = subplot(nsp(1),nsp(2),nc);          
                
                %plot the image;                                       
                imagesc(obj.bincenters_x(500),obj.bincenters_y(500),obj.stimulus);
                hold on;
                h     = imagesc(obj.bincenters_x,obj.bincenters_y,M(:,:,nc),[d u]);
                set(h,'alphaData',Scale(abs((obj.maps(:,:,nc))))*.9+.1);               
                axis image;
                axis off;
                try
                t     = sprintf('%s%d/',obj.map_titles{nc}{:});                
                title(t,'interpreter','none');
                end
            end
            %             thincolorbar('vert');
            colorbar
        end
        function getmaps(obj,varargin)
            %will populate the maps property based on the filter in
            %VARARGIN, example:
            %fs.getmaps({'phase',4,'deltacsp',-135});
            
            obj.maps = [];%clear whatever is there
            c  = 0;
            obj.all_queries = varargin;
            for v = varargin
                c                 = c+1;
                obj.UpdateSelection(v{1}{:});                
                if strcmp(obj.maptype,'conv')
                    %accum and conv
                    FixMap            = accumarray([obj.current_y-obj.rect(1)+1 obj.current_x-obj.rect(2)+1],1,[obj.rect(3) obj.rect(4)]);
                    FixMap            = conv2(sum(obj.kernel),sum(obj.kernel,2),FixMap,'same');
                elseif strcmp(obj.maptype,'bin')
                    %divide by a factor (i.e. binning) and accum, but no conv                                        
                    [FixMap]          = hist3([obj.current_y obj.current_x],'edges',obj.binedges);
                    FixMap            = FixMap/obj.current_ttrial;%divide by the number of trials
                    %remove the last column and row
                    FixMap(end,:)     = [];
                    FixMap(:,end)     = [];
                end
                %
                if obj.unitize
                    FixMap            = FixMap./sum(FixMap(:));
                end
                obj.maps(:,:,c)       = FixMap;
                %
                obj.map_titles{c}     = obj.query(:);
            end
            %correct for baseline if wanted.
            if obj.bc
                obj.maps = obj.maps - repmat(mean(obj.maps,3),[1 1 size(obj.maps,3)]);
            end            
        end                      
        function bild = get.stimulus(obj)
            
            bild    = imread(obj.find_stim);            
            bild    = bild( obj.rect(1):obj.rect(1)+obj.rect(3)-1,  obj.rect(2):obj.rect(2)+obj.rect(4)-1);
            bild    = repmat(bild,[1 1 3]);
        end        
        function out = cropmaps(obj,map)
            if strcmp(obj.maptype,'conv')
                out      = map(obj.rect(2)/2-obj.mapsize : obj.rect(2)/2+obj.mapsize, obj.rect(4)/2-obj.mapsize: obj.rect(4)/2+obj.mapsize);
            elseif strcmp(obj.maptype,'bin')
                new_size = obj.mapsize./obj.binsize;
                
                out      = map(obj.rect(2)/obj.binsize/2-new_size : obj.rect(2)/obj.binsize/2+new_size, obj.rect(4)/obj.binsize/2-new_size: obj.rect(4)/obj.binsize/2+new_size);
            end
        end        
        function out = kernel(obj)            
            dummy       = make_gaussian2D(obj.kernel_fwhm*2.65,obj.kernel_fwhm*2.65,obj.kernel_fwhm,obj.kernel_fwhm,obj.kernel_fwhm*2.65/2,obj.kernel_fwhm*2.65/2);
            out         = dummy./sum(dummy(:));       
        end
        function maps   = vectorize_maps(obj)
            if ~isempty(obj.maps)
                maps = reshape(obj.maps,size(obj.maps,1)*size(obj.maps,2),size(obj.maps,3));
            else
                fprintf('no maps stored here\n');
            end
        end
        function [out] = subplot_number(self);
            % how many subplots
            [row col] = GetSubplotNumber(size(self.maps,3));
            out       = sort([row col]);
        end
        function out = cov(obj)
%             figure(3);clf
%             set(gcf,'position',[1104         152         337         299]);
            out = cov(obj.vectorize_maps);
%             imagesc(out);
%             axis image;
%             colorbar
        end
        function out = corr(obj)
            %plots the corrmat between maps             
%             figure(2);clf
%             set(gcf,'position',[1104         152         337         299])            
            out = corr(obj.vectorize_maps);
%             imagesc(out);
%             axis image;
%             colorbar
        end        
        function [x] = current_x(obj)
            %returns the current x y coordinates            
            x = obj.x(obj.selection)';
        end
        function [y] = current_y(obj)
            %returns the current x y coordinates
            y = obj.y(obj.selection)';            
        end
        function ttrial = current_ttrial(obj)
            %computes number of trials included in the current selection
            ttrial = length(unique([obj.trialid(obj.selection) ;obj.subject(obj.selection)]','rows'));
        end
        function out = get.rect(obj)
            out = [obj.screen_resolution(1)/2-obj.window obj.screen_resolution(2)/2-obj.window [obj.window obj.window]*2];
        end
        function out = get.binedges(obj)
            out = {obj.rect(1):obj.binsize:(obj.rect(1)+obj.rect(3)) ...
                               obj.rect(2):obj.binsize:(obj.rect(2)+obj.rect(4))};                           
        end
        function [out] = bincenters_y(obj,varargin)
            %
            if nargin < 2
                resolution = length(obj.binedges{1});
            else
                resolution = varargin{1};
            end
            %transforms edges to centers with a required resolution.
            out = linspace( min(obj.binedges{1}), max(obj.binedges{1}),resolution);
            out = out + mean(unique(diff(out)))/2;
            out(end) =[];
        end
        function [out] = bincenters_x(obj,varargin)
            if nargin < 2
                resolution = length(obj.binedges{2});
            else
                resolution = varargin{1};
            end
            %transforms edges to centers with a required resolution.            
            out = linspace( min(obj.binedges{2}), max(obj.binedges{2}),resolution);            
            out = out + mean(unique(diff(out)))/2;
            out(end) =[];
        end
        function [out] = mapsize(obj,varargin)
            if ~isempty(obj.maps)
                %returns size of fixmap along VARARGIN dimension
                out = size(obj.maps,varargin{1});
            else 
                out =[];
            end
        end
        function [count index] = histogram(obj)
            %will generate an histogram of fixation numbers
            subject    = unique(obj.subject);
            phase      = unique(obj.phase);
            condition  = unique(obj.deltacsp);
            count      = nan(length(subject),length(condition));
            sub = 0;
            for ns = subject
                sub = sub + 1;
                index.sub(sub) = ns;
                ph  = 0;
                for np = phase
                    ph = ph + 1;
                    cond = 0;                    
                    for nc = condition                        
                        cond = cond + 1;
                        index.cond(cond) = double(nc);%column identities
                        obj.UpdateSelection('subject',ns,'phase',np,'deltacsp',nc);
                        %number of trials
                        repet              = obj.current_ttrial;
                        %average number of fixation
                        count(sub,cond,ph) = sum(obj.selection)./repet;
                    end
                end
            end
            figure;
            set(gcf,'position',[440   393   920   405]);
            for phases = 1:ph
                subplot(1,ph,phases)
                imagesc(count(:,:,phases),[0 7]);
                thincolorbar('vert');
                ylabel('subjects')
                xlabel('condition');
                axis image;
                box off;
                title(sprintf('Phase: %03d',phases))
            end
            supertitle('Fixation Counts',1);
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
