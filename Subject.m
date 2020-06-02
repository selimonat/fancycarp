classdef Subject < Project
    %
    %   SUBJECT object defines a set of subject-related methods.
    %
    %   There are Getter methods. They get data. They are prefixed with
    %   get_ string. For example, get_total_run, would return the number of
    %   runs that are present in the subject's folder, get_hr or get_epi,
    %   would download data from the dicom server., get_motion parameters,
    %   would return estimated motion parameters, etc.
    %
    %   Preprocessing methods is the core of the pipeline and provides a
    %   standardized spm preprocessing pipeline.
    %
    %   PathTool methods returns the path to different files, for example
    %   path_epi(1), would return the path to 4D data file in run 1.
    %
    %   Analysis methods do things like first-level analysis. In order this
    %   function flawlessly with your project you will need to create and
    %   save onsets in a directory called subXXX/design/modelXX. You can
    %   save different onsets under different models and run first-level
    %   analyses on a specific model.
    %
    %   Plotter methods plot things as their name suggests.
    %
    %
    properties (Hidden)
        paradigm
        default_run       = 1;
        dicom_serie_id    = [];
        dicom_folders     = [];
        dicom_target_run  = [];
        derivatives       = [0 0];%specifies expansion degree of the cHRF when running models.
        scr_ok            = [];
        kill_unconfirmed  = 0;
        kill_nan          = 1;
    end
    properties (SetAccess = private)
        id
        path
        scr
        csp           = [];
        selectedface  = [];
        trio_session  = [];
        hr_session    = [];
        total_run     = [];
        rating_fit    = [];
        is_tuned      = [];
        thresh        = [];
        pmf = [];
    end
    %%
    methods
        function s = Subject(id,varargin)%constructor
            fprintf('Subject Constructor for id:%i is called:\n',id);
            s.id               = id;
            s.path             = s.pathfinder(s.id,[]);
            s.dicom_serie_id   = s.dicom_serie_selector{s.id};
            s.dicom_target_run = s.dicom2run{s.id};
            s.dicom_target_run = s.dicom2run{s.id};
            s.trio_session 	  = s.trio_sessions{s.id};
            s.hr_session      = s.hr_sessions{s.id};
            
            if exist(s.path)
                for nrun = 1:s.total_run
                    s.paradigm{nrun} = s.get_paradigm(nrun);
                end
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                if any(cellfun(@(i) strcmp(i,'scr'),varargin))
                    s.scr   = SCR(s);
                end
                %                 scr_ok = load(s.path
                %                 s.scr = SCR(s);
            else
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path);
                fprintf('Your path might also be wrong...\n');
            end
        end
    end
    
    methods %(Getters)
        function get_hr(self)
            %will download the latest HR for this subject to default hr
            %path (which is run000)
            %
            
            fprintf('get_hr:\nWill now dump the latest HR (%s)\n',self.current_time);
            %target location for the hr: self.dir_hr;
            %create it if necess.
            if exist(self.dir_hr) == 0
                mkdir(self.dir_hr);
            end
            self.DicomDownload(self.path_hr_dicom,self.dir_hr);
            self.ConvertDicom(self.dir_hr);
            files       = spm_select('FPListRec',self.dir_hr,'^sPRISMA');
            if ~isempty(files)
                movefile(files,regexprep(files,sprintf('%ssPRISMA.*$',filesep),sprintf('%sdata.nii',filesep)));%rename it to data.nii
            end
        end
        function p          = get_paradigm(self,nrun)
            %will load the paradigm file saved during your psychophysics
            %session.
            filename = self.path_data(nrun,'stimulation');
            p = [];
            if exist(filename)
                p = load(filename);
                p = p.p;
            end
        end
        function get_epi(self)
            %Will dump all DICOMS based on Sessions entered in the
            %Project object. trio_folders are folders in the dicom server,
            %trio2run dictates in which run these folders should be dumped
            %to.
            %
            
            %spit out some info for sanity checks
            self.dicomserver_request;
            fprintf('You told me to download the following series: ');
            fprintf('%i,',self.dicom_serie_id);
            fprintf('\nDouble check if everything is fine.\n');
            paths               = self.dicomserver_paths;
            if ~isempty(paths)
                self.dicom_folders  = paths(self.dicom_serie_id);
                fprintf('Will now dump series (%s)\n',self.current_time);
            end
            %% save the desired runs to disk
            n = 0;
            for source = self.dicom_folders(:)'
                %
                n 				 = n+1;
                dest             = self.dir_epi(n);
                if exist(dest)
                    self.DicomDownload(source{1},dest);
                    self.DicomTo4D(dest);
                else
                    fprintf('Stopped here as a sanity check\nIt seems the destination folder doesn''t exist.')
                    keyboard
                end
            end
        end
        function [o]    = get.total_run(self)
            %% returns the total number of runs in a folder (except run000)
            o      = length(dir(self.path))-4;%exclude the directories .., ., and run000 and run005 (bc no brain data)
        end
        function L      = get_log(self,run)
            % Loads the ptb log. Removes all events before/after the
            % first/last scan and defines zero as the first scan.
            %
            
            L               = self.paradigm{run}.out.log;
            %sort things according to time rather than order of being logged
            [~,i]           = sort(L(:,1),'ascend');
            L               = L(i,:);
            % delete all the events that are after the last scanning..
            scan_times      = L(find(L(:,2) == 0),1);
            %
            first_scan_time = min(scan_times);
            last_scan_time  = max(scan_times);
            %
            L(L(:,1) < first_scan_time,:) = [];
            %             L(L(:,1) > last_scan_time,:)  = [];
            L(:,1)          = L(:,1) - first_scan_time;
        end
        function o      = get_param_motion(self,run)
            %will load the realignment parameters, of course you have to
            %realign the EPIs first.
            
            filename = sprintf('%smrt%srp_data.txt',self.path_data(run),filesep);
            if exist(filename)
                o = load(filename);
            else
                fprintf('File:\n %s doesn''t exist.\n Most likely realignment is not yet done.\n',filename);
            end
        end
        function cond   = get_modelonsets(self,nrun,model_num)
            %returns stimulus onsets for NRUN defined by model specified by
            %MODEL_NUM
            cond = [];
            load(self.path_model(nrun,model_num));
        end
        function N      = get_param_nuissance(self,nrun)
            %Computes nuissances parameters from MotionParameters.
            N = self.get_param_motion(nrun);
            N = zscore([N [zeros(1,size(N,2));diff(N)] N.^2 [zeros(1,size(N,2));diff(N)].^2 ]);
        end
        function [t]    = get_total_volumes(self,run)
            % will tell you how many volumes are in a 4D image.
            bla = spm_vol_nifti(self.path_epi(run),1);%simply read the first images header
            t   = bla.private.dat.dim(4);
        end
        function [out, out_corr]   = get_totalvolumeslogged(self,run)
            %returns number of pulses logged in stimulus computer during the experiment
            
            verbalize = 1;
            
            L   = self.get_log(run);
            loggedpuls = L(:,2) == 0;
            out        = sum(loggedpuls);
            scan_times = L(loggedpuls,1);
            Ngaps   = sum(diff(scan_times) > self.TR*1.1);
            diffs      = sort(diff(scan_times),'descend');
            gap2pulse = [diffs(1:Ngaps)./self.TR]-1;
            addpuls   = sum(gap2pulse);
            out_corr = out + round(addpuls);
            if verbalize
                fprintf('Subject %g: \n',self.id)
                fprintf('Phase %g - %02d diffs bigger than TR, namely: ',run,Ngaps)
                disp(diffs(1:Ngaps))
                fprintf('Corresponds to N pulses: \n')
                disp(gap2pulse)
                fprintf('Total of %04.2f additional pulses.\n',addpuls)
            end
            if self.id == 6
                out_corr = out-7;
                if run == 4
                    out_corr = 890;
                end
            end
        end
        function [Nvols] = get_lastscan_cooldown(self,nrun)
            modelnum = 4; %get condfile from there and find out where CoolDown phase started
            a = load(self.path_model(nrun,modelnum));
            cond = a.cond;
            if self.id == 32 && nrun == 1
                Nvols = 409;
            else
                Nvols = floor(cond(end).onset);
            end
        end
        function XYZmm  = get_nativeatlas2mask(self,mask_id)
            %Will return XYZ coordinates from ROI specified by MASK_INDEX
            %thresholded by the default value. XYZ values are in world
            %space, so they will need to be brought to the voxel space of
            %the EPIs.
            mask_handle = spm_vol(self.path_native_atlas(mask_id));%read the mask
            mask_ind    = spm_read_vols(mask_handle) > self.atlas2mask_threshold;%threshold it
            [X Y Z]     = ind2sub(mask_handle.dim,find(mask_ind));%get voxel indices
            XYZ         = [X Y Z ones(sum(mask_ind(:)),1)]';%this is in hr's voxels space.
            XYZmm       = mask_handle.mat*XYZ;%this is world space.
            XYZmm       = unique(XYZmm','rows')';%remove repetititons.
        end
        function L      = get_physio2log(self)
            %returns logged events from the physio computer in the same
            %format as the log file. Events are aligned to the first valid
            %scan pulse.
            %%
            
            filename = sprintf('%s/run001/scr/data.smr',self.path);
            fh       = fopen(filename);
            % get times for trigger channels.
            chan2L   = [NaN 9 3 NaN 5 0 NaN NaN 5];%transform channels to event_types as logged by the stim computer.
            L        = [];
            for chan     = [2 3 4 5 6 9];%all event channels.
                dummy    = SONGetEventChannel(fh,chan);
                L        = [L ;[dummy(:) repmat(chan2L(chan),length(dummy),1)]];%returns time in seconds.
            end
            % include only the period starting and ending with pulses. To
            % this end find the period where the distance between the two
            % pulse is smaller than the TR times 1.1. In the physio
            % computer events might have been recorded before and after
            % scanning session.
            scan_times      = L(find(L(:,2) == 0),1);
            last_scan_ind   = find(diff(scan_times) < self.TR*1.1,1,'last');
            first_scan_ind  = find(diff(scan_times) < self.TR*1.1,1,'first');
            %
            first_scan_time = scan_times(first_scan_ind);
            last_scan_time  = scan_times(last_scan_ind);
            %
            L(L(:,1) < first_scan_time,:) = [];
            L(L(:,1) > last_scan_time,:)  = [];
            %
            L(:,1)   = L(:,1) - first_scan_time;
            %
            fclose(fh);
        end
        
        
        function selected = get.selectedface(self)
            selected = self.get_paradigm(5).out.selectedface;
        end
        
        function scr_ok = get_scr_ok(self)
            path =fullfile(self.pathfinder(4,0),'midlevel','scr_ok.mat');
            if ~exist(path)
                path2xls = fullfile(self.path_project,'midlevel','scr_ok_info.mat');
                xlsread(path2xls);
                scr_ok = [];
            else
                load(path);
            end
        end
    end
    methods %(behavioral analysis)
        function [threshold, completelist] = get_threshold(self)
            path2calib = fullfile(self.pathfinder(self.id,0),'calib2','stimulation','data.mat');
            load(path2calib);
            threshold = p.presentation.limits.threshold_ave;
            self.thresh = threshold;
            if nargout > 1
                for n = 1:2
                    path2p = fullfile(self.pathfinder(self.id,0),sprintf('calib%d',n),'stimulation','data.mat');
                    load(path2p)
                    completelist(n,:) = [p.presentation.limits.threshold_m p.presentation.limits.threshold_ave];
                end
            end
        end
        function out = deltacsp2faceID(self,deltacsplist)
            old_dim = size(deltacsplist);
            vec_list = deltacsplist(:);
            out = nan(size(vec_list));
            logf = self.get_paradigm(3);
            list = [logf.presentation.dist' logf.presentation.stim_id'];
            lookup = unique(list,'rows');
            for cond =unique(deltacsplist)
                out(deltacsplist==cond) = lookup(lookup(:,1)==cond,2);
            end
            out = reshape(out,old_dim);
        end
        function [out, raw, triallist] = get_rating(self,varargin)
            % this function is mostly to get ratings so that we can fit
            % it.. OUT is a structure feedable to the Tuning object, so
            % out.y, out.x, out.id, while RAW gives the
            % ratings for each trial (cond indicated in TRIALLIST).
            out = [];
            raw = [];
            if isempty(varargin)
                runs = 1:5;
                fprintf('No run selected, collecting all runs.\n')
            else
                runs = varargin{1};
            end
            if nargin > 2
                conds = varargin{2};
            else
                conds = self.realconds;
            end
            
            for run = runs(:)'
                if run < 5
                    try
                        a = self.get_paradigm(run);
                        dummy = a.log.ratings.relief;
                        raw{run} = dummy(:,3);
                        triallist{run} = a.presentation.dist';
                        confirmed{run} = a.log.ratings.relief(:,4);
                    catch
                        warning('Problem while loading paradigm for sub %d at run %d.',self.id, run)
                        raw{run} = NaN;
                        triallist{run} = NaN;
                        confirmed{run} = NaN;
                    end
                elseif run == 5
                    a = self.get_paradigm(3);
                    raw{run} = a.log.ratings.relief(:,3);
                    triallist{run} = a.presentation.dist(:);
                    confirmed{run} = a.log.ratings.relief(:,4);
                    try
                        a = self.get_paradigm(4);
                        raw{run} = [raw{run}; a.log.ratings.relief(:,3)];
                        triallist{run} = [triallist{run}; a.presentation.dist'];
                        confirmed{run} = [confirmed{run}; a.log.ratings.relief(:,4)];
                    end
                end
            end
            %now create the out struct
            for r = 1:length(raw)
                out{r}.x = [];
                out{r}.y = [];
                out{r}.ids = [];
                conf{r} = [];
                cc = 0;
                for cond = conds
                    cc = cc+1;
                    ind = triallist{r} == cond;
                    out{r}.x = [out{r}.x repmat(cond,1,sum(ind))];
                    out{r}.y = [out{r}.y raw{r}(ind)'];
                    out{r}.ids = self.id;
                    conf{r} = [conf{r} confirmed{r}(ind)'];
                end
            end
            
            if numel(runs) == 1
                out = out{end};
                raw = raw{end};
                triallist = triallist{end};
                confirmed = conf{end};
            end
            if self.kill_unconfirmed ==1
                warning('Kicking unconfirmed trials!')
                if numel(runs)==1
                    out.x(confirmed==0)=[];
                    out.y(confirmed==0)=[];
                else
                    for r = 1:numel(runs)
                        out{r}.x(confirmed==0)=[];
                        out{r}.y(confirmed==0)=[];
                    end
                end
            end
            if self.kill_nan
                if numel(runs)==1
                    warning('Kicking %d Nans out for sub %d, run %d.',sum(isnan(out.y)),self.id, runs)
                    out.x(isnan(out.y))=[];
                    out.y(isnan(out.y))=[];
                else
                    for r = 1:numel(runs)
                        warning('Kicking %d Nans out for sub %d, run %d.',sum(isnan(out{r}.y)),self.id, r)
                        out{r}.x(isnan(out{r}.y))=[];
                        out{r}.y(isnan(out{r}.y))=[];
                    end
                end
            end
        end
        
        function [out, trialind, temp] = get_pain(self,varargin)
            out = [];
            trialind = [];
            if isempty(varargin)
                %                 fprintf('No run selected, collecting all runs.\n')
                for run = 1:4
                    try
                        a = self.get_paradigm(run);
                        dummy = a.log.ratings.pain;
                        temp(run) = dummy(1,1); %get the stored temperature
                        dummy = dummy(~isnan(dummy(:,3)),3); %get rid of all unnecessary data, just take pain rating
                        if run < 3
                            out(:,run) = [dummy(1:3); NaN];
                        else
                            out(:,run) = dummy(1:4);
                        end
                        findtrial = find(a.presentation.ratepain(:));
                        trialind{run} = [findtrial' a.presentation.tTrial];
                    catch
                        warning('Problem while loading paradigm for sub %d at run %d.',self.id, run)
                        out(:,run) = nan(4,1);
                        temp(run) = nan;
                    end
                end
            else
                runs = varargin{1};
                rc = 0;
                for run = runs(:)'
                    rc = rc+1;
                    fprintf('Collecting ratings for run %d.\n',run)
                    try
                        a = self.get_paradigm(run);
                        dummy = a.log.ratings.pain;
                        temp(rc) = dummy(1,1); %get the stored temperature
                        dummy = dummy(~isnan(dummy(:,3)),3);
                        if run < 3
                            out = dummy(1:3);
                        else
                            out = dummy(1:4);
                        end
                        findtrial = find(a.presentation.ratepain(:));
                        trialind = [findtrial' a.presentation.tTrial];
                    catch
                        warning('Problem while loading paradigm at run %d.',run)
                        out = NaN;
                        temp = NaN;
                    end
                end
            end
            
        end
        function [PRindex,paininterp,relief] = get_pmod_PR(self,run)
            [pain,paintrials] = self.get_pain(run);
            [~,relief] = self.get_rating(run);
            
            paininterp = interp1(paintrials,pain,1:length(relief),'linear');
            PRindex = paininterp./100.*relief';
            %zscoring and nulltrial-kick will be done in PrepareCond method
        end
        function [M, S] = get_reliefmeans(self,run,varargin)
            conds = self.allconds;
            if nargin > 2
                conds = varargin{1};
            end
            if nargin > 3
                corrtype = varargin{2};
            else
                corrtype = 'raw';
            end
            out = self.get_rating(run,conds);
            
            if strcmp(corrtype,'mc')
                
                if run == 5
                    out.y = [self.get_rating(3,conds).y-nanmean(self.get_rating(3,conds).y) self.get_rating(4,conds).y-nanmean(self.get_rating(4,conds).y)];
                else
                    out.y = out.y - nanmean(out.y);
                end
                %LK check here
            elseif strcmp(corrtype,'zscore')
                if length(conds)<10
                    warning('ZSCORE might be computed on SELECTION of conds')
                end
                if run == 5
                    out.y = [nanzscore(self.get_rating(3,conds).y) nanzscore(self.get_rating(4,conds).y)];
                    out.x = [self.get_rating(3,conds).x self.get_rating(4,conds).x]; %otherwise its sorted by condition and not matching the y anymore.
                else
                    out.y = nanzscore(out.y);
                end
            elseif strcmp(corrtype,'raw')
                out.y = out.y;
            else
                warning('Unknown correction type as input, returning raw values.')
                out.y = out.y;
            end
            
            M = [];
            S = [];
            
            
            cc = 0;
            for cond = conds(:)'
                cc = cc+1;
                ind = out.x == cond;
                M(cc) = nanmean(out.y(ind));
                S(cc) = nanstd(out.y(ind));
            end
        end
        function [ratings] = get_relief_percond(self,run,varargin)
            conds = self.realconds;
            if nargin > 2
                type = varargin{1};
            else
                type = 'raw';
            end
            if strcmp(type,'raw')
                relief = self.get_rating(run,conds);
                reps = histc(relief.x,unique(relief.x));
                ratings = nan(max(reps),8);
                nc = 0;
                for cond = conds(:)'
                    nc = nc+1;
                    ratings(:,nc) = relief.y(relief.x == cond);
                end
            elseif strcmp(type,'zscore')
                if run == 5
                    run = [3 4];
                    if self.id == 15
                        run = 3;
                    end
                end
                ratings = [];
                rc = 0;
                for nr = run(:)'
                    rc = rc+1;
                    clear run_ratings
                    relief = self.get_rating(nr,self.allconds);
                    relief.y = nanzscore(relief.y);
                    reps = histc(relief.x,unique(relief.x));
                    run_ratings = nan(max(reps),length(conds));
                    nc = 0;
                    for cond = conds(:)'
                        nc = nc+1;
                        run_ratings(1:reps(nc),nc) = relief.y(relief.x == cond);
                    end
                    ratings = [ratings; run_ratings];
                end
            end
        end
        function [out, R] = fit_rating(self,run,varargin)
            %will load the rating fit (saved in runXXX/rating) if computed other
            %wise will read the raw ratingdata (saved in runXXX/stimulation)
            %and compute a fit.
            
            if nargin > 2
                datatype = varargin{1};
            else
                datatype = 'zscore';
            end
            
            force      = 1;%repeat the analysis or load from cache
            fun        = self.selected_fitfun;
            write_path = sprintf('%s/midlevel/rating_fun_%i_%s.mat',self.pathfinder(self.id,run),fun,datatype);
            
            
            if exist(write_path) && force ==0
                %load directly or
                load(write_path);
                fprintf('Rating Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif force == 1 || ~exist(write_path)
                %compute and save it.
                fprintf('Fitting Rating for run %02d..\n',run)
                if strcmp(datatype,'raw')
                    R = self.get_rating(run);
                elseif strcmp(datatype,'zscore')
                    R.y = self.get_relief_percond(run,'zscore');
                    R.y(isnan(R.y))=0;
                    R.x = repmat(self.realconds,size(R.y,1),1);
                    R.ids = self.id;
                    R.y = R.y(:)';
                    R.x = R.x(:)';
                elseif strcmp(datatype,'zscore_bc')
                    R.y_t = self.get_relief_percond(run,'zscore');
                    R.y_b = self.get_relief_percond(1,'zscore');
                    R.y   = nanmean(R.y_t,1)-nanmean(R.y_b,1);
                    R.x = self.realconds;
                    R.ids = self.id;
                end
                
                if isempty(R.y)
                    warning('No ratings, so no fit for this person and run.\n');
                    out = [];
                else
                    %create a tuning object and make a single subject fit
                    T                            = Tuning(R);
                    T.SingleSubjectFit(fun);
                    %prepare data for outputting.
                    out.params                   = T.fit_results.params(1,:);
                    out.LL                       = T.fit_results.Likelihood(1,:);
                    out.pval                     = T.fit_results.pval(1,:);
                    out.exitflag                 = T.fit_results.ExitFlag(1,:);
                    out.y_hd                     = T.fit_results.fit_HD(1,:);
                    out.x_hd                     = T.fit_results.x_HD(1,:);
                    out.y                        = T.fit_results.y_fitted(1,:);
                    out.x                        = T.fit_results.x(1,:);
                    %                     out.fitfun                   = T.fit_results.fitfun;
                    %
                    if fun == 8%if vM, then transform kappa to FWHM.
                        fprintf('Doing kappa to FWHM and Amplitude (+/-) transformation.\n');
                        out.params(:,2)            = vM2FWHM(out.params(:,2));
                        %out.params(:,3)            = abs(out.params(:,3));
                        
                        %also put the amplitude to negative if CS- > CS+
                        if out.y(4) < out.y(8)
                            out.params(1) = -out.params(1);
                        end
                    end
                end
                save(write_path,'out','R')
            end
            
            self.rating_fit{run} = out;
            self.rating_fit{run}.data = R;
            
            self.is_tuned{run} = (10.^-out.pval) < .05;
        end
        function plot_ratings(self)
            force = 0;
            corrtype = 'raw';
            
            savepath = [self.path 'run005/figures/'];
            if ~exist(savepath)
                mkdir(savepath)
            end
            savefigf = [savepath sprintf('ratingplot_sub%02d_%s.fig',self.id,corrtype)];
            savebmp = [savepath sprintf('ratingplot_sub%02d_%s.bmp',self.id,corrtype)];
            savepng = [savepath sprintf('ratingplot_sub%02d_%s.png',self.id,corrtype)];
            if exist(savefigf,'file') && force == 0
                fprintf('Found saved .fig file, opening it as figure %02d\n',self.id)
                openfig(savefigf);
            else
                
                titles = {'Base','Cond','Test1','Test2','Test1/2'};
                f=figure(self.id);
                f.Position = [73 410 1763 476];
                clf;
                
                for run = 1:5
                    [M,S] = self.get_reliefmeans(run,self.allconds,corrtype);
                    subplot(1,5,run)
                    self.plot_bar(self.plotconds,M,S)
                    hold on
                    axis square
                    set(gca,'XTickLabels',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 't0'},'XTickLabelRotation',45,'FontSize',12);
                    xlim([min(self.plotconds)-40 max(self.plotconds+40)])
                    title(titles{run},'FontSize',14)
                    if run ~=2
                        self.fit_rating(run); % for cond, this doesn't make sense
                        if self.rating_fit{run}.pval < .05
                            plot(self.rating_fit{run}.x_hd,self.rating_fit{run}.y_hd,'k-','LineWidth',3)
                        else
                            txt= text(min(xlim)+20,M(1)+S(1).*1.2,sprintf('p = %4.3f',self.rating_fit{run}.pval));
                            set(txt, 'rotation', 90)
                            l=line([-135,180],repmat(mean(self.rating_fit{run}.data.y),1,2));
                            set(l,'LineWidth',3,'Color','k')
                        end
                    end
                end
                st = supertitle(sprintf('sub %02d, csp = %d ',self.id, self.csp));
                set(st,'FontSize',16);
                EqualizeSubPlotYlim(gcf);
                
                savefig(gcf,savefigf)
                export_fig(gcf,savebmp)
                export_fig(gcf,savepng,'-transparent')
            end
        end
        function [delta,temp] = get_ramptemp(self,run)
            tTrial = self.paradigm{run}.presentation.tTrial;
            delta = nan(1,tTrial);
            temp  = nan(1,tTrial);
            delta(self.paradigm{run}.presentation.ucs)  = self.paradigm{run}.presentation.pain.tonic(1) - self.paradigm{run}.presentation.pain.low(1);
            delta(~self.paradigm{run}.presentation.ucs) = self.paradigm{run}.presentation.pain.tonic(1) - self.paradigm{run}.presentation.pain.middle(1);
            temp(self.paradigm{run}.presentation.ucs)   = self.paradigm{run}.presentation.pain.low(1);
            temp(~self.paradigm{run}.presentation.ucs)  = self.paradigm{run}.presentation.pain.middle(1);
            
        end
    end
    %%
    methods %(mri, preprocessing))
        function preprocess_pipeline(self,runs)
            %meta method to run all the required steps for hr
            %preprocessing. RUNS specifies the functional runs, make it a
            %vector if needed.
            if nargin > 1
                self.SegmentSurface_HR;%cat12 segmentation
                self.SkullStrip;%removes non-neural voxels
                %             	self.MNI2Native;%brings the atlas to native space
                self.Re_Coreg(runs);%realignment and coregistration
                self.Segment_meanEPI;%segments mean EPI with new segment
                self.SkullStrip_meanEPI;%creates a native mask
            else
                fprintf('One input argument is required!\n');
            end
        end
        function SkullStrip_meanEPI(self)
            %combines c{1,2,3}meanepi.nii images to create a binary mask in the native space.
            Vi = self.path_meanepi_segmented(1:3);
            Vo = self.path_skullstrip_meanepi;
            spm_imcalc(Vi,Vo,'(i1+i2+i3)>0');
        end
        function SkullStrip(self)
            %needs results of SegmentSurface, will produce a skullstripped
            %version of hr (filename: ss_data.nii). It will also
            %automatically create a normalized version as well
            %(w_ss_data.nii).
            %c
            
            %path to p1 and p2 images created by SegmentSurface
            c1         = strrep(self.path_hr,sprintf('mrt%sdata',filesep),sprintf('mrt%smri%sp1data',filesep,filesep));
            c2         = strrep(self.path_hr,sprintf('mrt%sdata',filesep),sprintf('mrt%smri%sp2data',filesep,filesep));
            
            if exist(c1) && exist(c2)
                matlabbatch{1}.spm.util.imcalc.input            = cellstr(strvcat(self.path_hr,c1,c2));
                matlabbatch{1}.spm.util.imcalc.output           = self.path_skullstrip;
                matlabbatch{1}.spm.util.imcalc.outdir           = {self.dir_hr};
                matlabbatch{1}.spm.util.imcalc.expression       = 'i1.*((i2+i3)>0.2)';
                matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask     = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp   = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype    = 4;
                self.RunSPMJob(matlabbatch);
                %normalize the skull stripped image as soon as it is
                %created.
                self.VolumeNormalize(self.path_skullstrip);
            else
                fprintf('Need to run segment first...\n')
            end
        end
        function Re_Coreg(self,runs)
            %will realign and coregister.
            
            %% collect all the EPIs as a cell array of cellstr
            c = 0;
            for nr = runs
                if exist(self.path_epi(nr))
                    c = c +1;
                    epi_run{c} = cellstr(spm_select('expand',self.path_epi(nr)));
                end
            end
            %% path to the mean_epi
            mean_epi    = self.path_meanepi;
            
            %double-pass realign EPIs and reslice the mean image only.
            matlabbatch{1}.spm.spatial.realign.estwrite.data = epi_run;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;%double pass
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];%reslice only the mean image.
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            
            %%coregister EPIs to skullstrip (only the affine matrix is modified)
            matlabbatch{2}.spm.spatial.coreg.estimate.ref    = cellstr(self.path_skullstrip);
            matlabbatch{2}.spm.spatial.coreg.estimate.source = cellstr(mean_epi);
            matlabbatch{2}.spm.spatial.coreg.estimate.other  = vertcat(epi_run{:});
            matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            % %
            %%write EPIs
            matlabbatch{3}.spm.spatial.realign.write.data            = vertcat(epi_run{:});
            matlabbatch{3}.spm.spatial.realign.write.roptions.which  = [2 1];%all images as well as the mean image.
            matlabbatch{3}.spm.spatial.realign.write.roptions.interp = 4;
            matlabbatch{3}.spm.spatial.realign.write.roptions.wrap   = [0 0 0];
            matlabbatch{3}.spm.spatial.realign.write.roptions.mask   = 1;
            matlabbatch{3}.spm.spatial.realign.write.roptions.prefix = 'r';
            self.RunSPMJob(matlabbatch);
            
        end
        function SegmentSurface_HR(self)
            % Runs CAT12 Segment Surface routine.
            % Will write to the disk:
            % mri and report folders
            % mri/iy_data.nii
            % mri/p1data.nii
            % mri/p2data.nii
            % mri/wmdata.nii
            % mri_y_data.nii
            matlabbatch{1}.spm.tools.cat.estwrite.data = {spm_select('expand',self.path_hr)};
            matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
            matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {sprintf('%sTPM.nii',self.tpm_dir)};
            matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
            matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1;
            matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
            matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 0.5;
            matlabbatch{1}.spm.tools.cat.estwrite.extopts.cleanupstr = 0.5;
            matlabbatch{1}.spm.tools.cat.estwrite.extopts.darteltpm = {self.dartel_templates(1)};
            matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
            matlabbatch{1}.spm.tools.cat.estwrite.output.surface = self.surface_wanted;
            matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
            matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0;
            matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
            matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
            matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
            matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
            matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
            matlabbatch{1}.spm.tools.cat.estwrite.output.jacobian.warped = 0;
            matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
            %
            self.RunSPMJob(matlabbatch);
        end
        function SegmentHR_SPMroutine(self)
            %based on mean epi segmentation script by CB.
            %prepare folder and file
            %original HR from run000
            hrfile = self.path_hr;
            % create new EPI Dartel thing folder
            hr_dartel_dir      = strrep( self.path_hr,sprintf('mrt%sdata.nii',filesep),sprintf('mrt%ssegm_dartel',filesep));
            %make and move.
            mkdir(hr_dartel_dir);
            copyfile(hrfile,[hr_dartel_dir filesep 'data.nii'])
            
            %------------------------
            % do the segmentation
            template = sprintf('%sTPM.nii',self.tpm_dir);
            
            sourcefile   = [hr_dartel_dir filesep 'data.nii'];
            
            matlabbatch{1}.spm.spatial.preproc.channel.vols     = cellstr(sourcefile);
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm    = {[template ',1']};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus  = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm    = {[template ',2']};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus  = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm    = {[template ',3']};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus  = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm    = {[template ',4']};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus  = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm    = {[template ',5']};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus  = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm    = {[template ',6']};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus  = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write       = [1 1];
            
            % DARTEL norm to template
            rc1_templ         = 'rc1data.nii';
            rc2_templ         = 'rc2data.nii';
            rc1_file      = [hr_dartel_dir filesep rc1_templ];
            rc2_file      = [hr_dartel_dir filesep rc2_templ];
            
            
            matlabbatch{2}.spm.tools.dartel.warp1.images = {cellstr(rc1_file),cellstr(rc2_file)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.rform = 0;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).K = 0;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).template = {self.dartel_templates(1)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).K = 0;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).template = {self.dartel_templates(2)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).K = 1;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).template = {self.dartel_templates(3)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).K = 2;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).template = {self.dartel_templates(4)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).K = 4;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).template = {self.dartel_templates(5)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).K = 6;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).template = {self.dartel_templates(6)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.optim.its = 3;
            
            %-------------------------------
            %Create warped HR
            
            u_rc1_templ       = 'u_rc1data.nii';
            u_rc1_file  = [hr_dartel_dir filesep u_rc1_templ];
            
            matlabbatch{3}.spm.tools.dartel.crt_warped.flowfields = cellstr(u_rc1_file);
            matlabbatch{3}.spm.tools.dartel.crt_warped.images = {cellstr(sourcefile)};
            matlabbatch{3}.spm.tools.dartel.crt_warped.jactransf = 0;
            matlabbatch{3}.spm.tools.dartel.crt_warped.K = 6;
            matlabbatch{3}.spm.tools.dartel.crt_warped.interp = 1;
            
            self.RunSPMJob(matlabbatch);
        end
        function SegmentMeanEPI_CB(self)
            %prepare folder and file
            %original mean epi from ph1 folder
            mean_epi    = strrep( self.path_epi(1),sprintf('mrt%sdata',filesep),sprintf('mrt%smeandata',filesep));
            % create new EPI Dartel thing folder
            mean_EPIdartel_dir      = strrep( self.path_epi(1),sprintf('mrt%sdata.nii',filesep),sprintf('mrt%smean_EPIdartel',filesep));
            %make and move.
            mkdir(mean_EPIdartel_dir);
            copyfile(mean_epi,[mean_EPIdartel_dir filesep 'meandata.nii'])
            
            %------------------------
            % do the segmentation
            template = sprintf('%sTPM.nii',self.tpm_dir);
            
            meanepi_file   = [mean_EPIdartel_dir filesep 'meandata.nii'];
            
            matlabbatch{1}.spm.spatial.preproc.channel.vols     = cellstr(meanepi_file);
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm    = {[template ',1']};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus  = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm    = {[template ',2']};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus  = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm    = {[template ',3']};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus  = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm    = {[template ',4']};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus  = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm    = {[template ',5']};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus  = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm    = {[template ',6']};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus  = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write       = [1 1];
            
            % DARTEL norm to template
            rc1_templ         = 'rc1meandata.nii';
            rc2_templ         = 'rc2meandata.nii';
            rc1_file      = [mean_EPIdartel_dir filesep rc1_templ];
            rc2_file      = [mean_EPIdartel_dir filesep rc2_templ];
            
            
            matlabbatch{2}.spm.tools.dartel.warp1.images = {cellstr(rc1_file),cellstr(rc2_file)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.rform = 0;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).K = 0;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(1).template = {self.dartel_templates(1)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).K = 0;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(2).template = {self.dartel_templates(2)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).K = 1;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(3).template = {self.dartel_templates(3)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).K = 2;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(4).template = {self.dartel_templates(4)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).K = 4;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(5).template = {self.dartel_templates(5)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).its = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).K = 6;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.param(6).template = {self.dartel_templates(6)};
            matlabbatch{2}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
            matlabbatch{2}.spm.tools.dartel.warp1.settings.optim.its = 3;
            
            %-------------------------------
            %Create warped mean EPI
            
            u_rc1_templ       = 'u_rc1meandata.nii';
            u_rc1_file  = [mean_EPIdartel_dir filesep u_rc1_templ];
            
            matlabbatch{3}.spm.tools.dartel.crt_warped.flowfields = cellstr(u_rc1_file);
            matlabbatch{3}.spm.tools.dartel.crt_warped.images = {cellstr(meanepi_file)};
            matlabbatch{3}.spm.tools.dartel.crt_warped.jactransf = 0;
            matlabbatch{3}.spm.tools.dartel.crt_warped.K = 6;
            matlabbatch{3}.spm.tools.dartel.crt_warped.interp = 1;
            
            self.RunSPMJob(matlabbatch);
        end
        function Segment_meanEPI(self)
            % Runs the new segment of SPM12 on the mean EPI image.
            % Will write to the disk:
            % meandata_seg8.mat
            % iy_meandata.nii
            % c{1-5}meandata.nii
            % y_meandata.nii
            
            if ~exist(self.dir_meanepi);mkdir(self.dir_meanepi);end
            if ~exist(self.path_meanepi);
                file_from_Coreg = strrep(self.path_meanepi,'meanEPI/','');
                if exist(file_from_Coreg)
                    copyfile(file_from_Coreg,self.path_meanepi)
                else
                    inp = input('Couldn''t find meandata.nii where I looked, do you want to run Re_Coreg again? Press y then. \n.','S')
                    if strcmp(inp,'y')
                        self.Re_Coreg(1:self.total_run);
                    else
                        return
                    end
                end
            end
            
            matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(self.path_meanepi);
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {self.path_tpm(1)};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {self.path_tpm(2)};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {self.path_tpm(3)};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {self.path_tpm(4)};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {self.path_tpm(5)};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {self.path_tpm(6)};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
            self.RunSPMJob(matlabbatch);
        end
        function VolumeNormalize(self,path2image)
            %SegmentSurface writes deformation fields (y_*), which are here used
            %to normalize the native hr images. Adds a prefix w- to
            %resampled images. path2image is the image to be resampled.
            %Example:
            %
            %s.VolumeNormalize(s.path_beta(1,1))
            %s.VolumeNormalize(s.path_skullstrip);
            
            %             %% Normalize with mean epi segmentation
            %             for nf = 1:size(path2image,1)
            %                 matlabbatch{nf}.spm.spatial.normalise.write.subj.def      = cellstr(regexprep(self.path_meanepi,'meandata','y_meandata'));
            %                 matlabbatch{nf}.spm.spatial.normalise.write.subj.resample = {path2image(nf,:)};
            %                 matlabbatch{nf}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70
            %                     78 76 85];
            %                 matlabbatch{nf}.spm.spatial.normalise.write.woptions.vox    = [Inf Inf Inf];
            %                 matlabbatch{nf}.spm.spatial.normalise.write.woptions.interp = 4;
            %                 matlabbatch{nf}.spm.spatial.normalise.write.woptions.prefix = 'wEPI_';
            %             end
            %             self.RunSPMJob(matlabbatch);
            %% Normalize with CAT12 segmentation
            matlabbatch =[];
            for nf = 1:size(path2image,1)
                matlabbatch{nf}.spm.spatial.normalise.write.subj.def      = cellstr(strrep(self.path_hr,'data.nii',sprintf('mri%sy_data.nii',filesep)));
                matlabbatch{nf}.spm.spatial.normalise.write.subj.resample = {path2image(nf,:)};
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70
                    78 76 85];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.vox    = [Inf Inf Inf];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.interp = 4;
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.prefix = 'wCAT_';
            end
            self.RunSPMJob(matlabbatch);
        end
        function MNI2Native(self)
            %brings the atlas to native space and saves it in run000/atlas.
            %Same as VolumeNormalize but uses the inverse deformation
            %fields but same batch. Currently functions only with the 120th
            %volume (which is right amygdala).
            if exist(self.path_atlas)
                %copy the atlas to subject's folder.
                copyfile(fileparts(self.path_atlas),[self.path_data(0) sprintf('atlas%s',filesep)]);
                %
                filename = self.path_native_atlas(120);%for test purposes and to gain speed I focus on amygdala right now.
                nf = 1;
                matlabbatch{nf}.spm.spatial.normalise.write.subj.def      = cellstr(strrep(self.path_hr,'data.nii',sprintf('mri%siy_data.nii',filesep)));%iy_ not y_!
                matlabbatch{nf}.spm.spatial.normalise.write.subj.resample = {filename};
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70
                    78 76 85];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.vox    = [Inf Inf Inf];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.interp = 4;
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.prefix = 'w';%inverse warped
                %
                self.RunSPMJob(matlabbatch);
                target_file = strrep(self.path_native_atlas,'data.nii','wdata.nii');%created by the above batch;
                movefile(target_file,self.path_native_atlas);
            else
                fprintf('For MNI2Native Analysis you need to have a atlas in %s\n',self.path_atlas);
            end
        end
        function DARTEL_warp_meanepi(self)
            clear matlabbatch
            matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {strrep(self.path_hr,'data.nii','segm_dartel/u_rc1data.nii')};
            matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{self.path_meanepi}};
            matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
            matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
            matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
            self.RunSPMJob(matlabbatch);
        end
        
        function Create_backwards_flowfield(self)
            clear matlabbatch
            cd(strrep(self.path_hr,'data.nii','segm_dartel'))
            matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {strrep(self.path_hr,'data.nii','segm_dartel/u_rc1data.nii')};
            matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
            matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
            matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {''};
            matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'backwards';
            matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.savepwd = 1;
            self.RunSPMJob(matlabbatch);
        end
        function [X,N,K]=spm_DesignMatrix(self,nrun,model_num)
            %will return the same design matrix used by spm in an efficient
            %way. see also: GetTimeSeries, spm_GetBetas
            
            %% Design matrix X
            cond                  = self.get_modelonsets(nrun,model_num);
            fMRI_T                = 16;
            fMRI_T0               = 1;
            xBF.T                 = fMRI_T;
            xBF.T0                = fMRI_T0;
            xBF.dt                = self.TR/xBF.T;
            xBF.UNITS             = 'scans';
            xBF.Volterra          = 1;
            xBF.name              = 'hrf';
            xBF                   = spm_get_bf(xBF);
            %
            for i = 1:length(cond);%one regressor for each condition
                Sess.U(i).dt        = xBF.dt;%- time bin (seconds)
                Sess.U(i).ons       = cond(i).onset;%- onsets    (in SPM.xBF.UNITS)
                Sess.U(i).name      = {sprintf('%02d',i)};%- cell of names for each input or cause
                %no parametric modulation here
                Sess.U(i).dur    =  repmat(0,length(Sess.U(i).ons),1);%- durations (in SPM.xBF.UNITS)
                Sess.U(i).P.name =  'none';
                Sess.U(i).P.P    =  'none';
                Sess.U(i).P.h    =  0;%- order of polynomial expansion
                Sess.U(i).P.i    =  1;%- sub-indices of u pertaining to P
            end
            %
            k                       = self.total_volumes(nrun);
            SPM.xBF                 = xBF;
            SPM.nscan               = k;
            SPM.Sess                = Sess;
            SPM.Sess.U              = spm_get_ons(SPM,1);
            %
            % Convolve stimulus functions with basis functions
            [X,Xn,Fc]               = spm_Volterra(SPM.Sess.U,SPM.xBF.bf,SPM.xBF.Volterra);
            % Resample regressors at acquisition times (32 bin offset)
            X                       = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
            %% Nuissance Parameters
            N                       = self.get_param_nuissance(nrun);
            %% Get high-pass filter
            %this is how it should be, but due to fearamy specificities, we
            %have to make work around. Note that this cannot be merge to
            %/mrt/xx or
            %% get the filtering strcture a la spm.
            run_borders              = [[0 910 910+895]+1;[910 910+895  self.total_volumes(nrun)]];
            K(1:size(run_borders,2)) = struct('HParam', self.HParam, 'row',    [] , 'RT',     self.TR ,'X0',[]);
            c = 0;
            for b = run_borders
                c        = c + 1;
                K(c).row = b(1);
                K(c)     = spm_filter(K(c));
                K(c).X0  = [ones(length(K(c).row),1)*std(K(c).X0(:)) K(c).X0];
            end
            
        end
        function beta = spm_GetBetas(self,nrun,model_num,mask_id)
            %will compute beta weights manually without calling SPM.
            [X N K ]  = self.spm_DesignMatrix(nrun,model_num);%returns the Design Matrix, Nuissiance Matrix, and High-pass Filtering Matrix
            Y         = self.TimeSeries(nrun,mask_id);
            Y         = zscore(Y);
            Y         = spm_filter(K,Y);%high-pass filtering.
            %
            DM        = [X N ones(size(X,1),1)];%append together Onsets, Nuissances and a constant
            DM        = spm_filter(K,DM);%filter also the design matrix
            DM        = spm_sp('Set',DM);
            DM        = spm_sp('x-',DM);% projector;
            beta      = DM*Y;
        end
        function [D,XYZvox] = TimeSeries(self,nrun,mask_id)
            %will read the time series from NRUN. MASK can be used to mask
            %the volume, otherwise all data will be returned (dangerous).
            
            vh          = spm_vol(self.path_epi(nrun,'r'));
            XYZmm       = self.get_nativeatlas2mask(mask_id);%in world space from native mask
            
            if spm_check_orientations(vh)
                XYZvox  = vh(1).mat\XYZmm;%in EPI voxel space.
                XYZvox  = unique(XYZvox','rows')';
            else
                keyboard;%sanity check;
            end
            XYZvox      = round(XYZvox);
            D           = spm_get_data(vh,XYZvox);
        end
        function check_mri_data(self)
            if self.id == 4
                
                files = {strrep(self.path_epi(1),'data.nii','meandata.nii');...
                    strrep(self.path_epi(1),'data.nii','rdata.nii');...
                    strrep(self.path_epi(2),'data.nii','rdata.nii');...
                    strrep(self.path_epi(3),'data.nii','rdata.nii')};
            else
                files = {strrep(self.path_epi(1),'data.nii','meandata.nii'),...
                    strrep(self.path_epi(1),'data.nii','rdata.nii'),...
                    strrep(self.path_epi(2),'data.nii','rdata.nii'),...
                    strrep(self.path_epi(3),'data.nii','rdata.nii')
                    strrep(self.path_epi(4),'data.nii','rdata.nii')};
            end
            self.CheckReg(files)
        end
    end
    methods %path_tools which are related to the subject
        function out = path_skullstrip_meanepi(self)
            out = regexprep(self.path_meanepi,'meandata','ss_meandata');
        end
        function out        = path_skullstrip(self,varargin)
            %returns filename for the skull stripped hr. Use VARARGIN to
            %add a prefix to the output, such as 'w' for example.
            if nargin == 1
                out = sprintf('%s%s',self.dir_hr,'ss_data.nii');
            elseif nargin == 2
                out = sprintf('%s%s_%s',self.dir_hr,varargin{1},'ss_data.nii');
            end
        end
        
        function out  = dir_meanepi(self)
            %returns the path to the meanepi (result of realignment).
            %returns empty if non-existent. Assumes that the first run contains the mean epi.
            first_run = self.dicom_target_run(1);
            out       = fullfile(self.pathfinder(self.id,first_run),'mrt','meanEPI');
        end
        function out  = path_meanepi(self)
            %returns the path to the meanepi (result of realignment).
            %returns empty if non-existent. Assumes that the first run contains the mean epi.
            out       = fullfile(self.dir_meanepi,'meandata.nii');
        end
        function out = path_tpm(self,n)
            %return the path to the Nth TPM image from the spm
            out = sprintf('%s/TPM.nii,%i',self.tpm_dir,n);
        end
        function out = path_meanepi_segmented(self,num)
            %Returns the path to the output of Segment_meanEPI, N can be a vector.
            mean_epi = self.path_meanepi;
            out      = '';
            for n = num(:)'
                out      = strvcat(out,regexprep(mean_epi,'meandata',sprintf('c%dmeandata',n)));
            end
        end
        
        function out        = dir_spmmat(self,nrun,model_num)
            %Returns the path to SPM folder in a given NRUN responsible for
            %the model MODEL_NUM. VARARGIN is used for the derivatives.
            out = sprintf('%sspm%smodel_%02d_chrf_%d%d%s',self.path_data(nrun),filesep,model_num,self.derivatives(1),self.derivatives(2),filesep);
            
        end
        function out        = path_model(self,run,model_num)
            %returns the path to a model specified by MODEL_NUM in run RUN.
            out = sprintf('%sdesign%smodel%02d%sdata.mat',self.path_data(run),filesep,model_num,filesep);
        end
        function out        = path_spmmat(self,nrun,model_num)
            %returns the path to spm folder for run RUN.
            dummy = self.dir_spmmat(nrun(1),model_num);
            out   = sprintf('%s%sSPM.mat',dummy,filesep);
        end
        function out        = path_native_atlas(self,varargin)
            %path to subjects native atlas, use VARARGIN to slice out a
            %given 3D volume.
            out = sprintf('%satlas%sdata.nii',self.path_data(0),filesep);
            if nargin > 1
                out = sprintf('%s,%d',out,varargin{1});
            end
        end
        function out        = dir_hr(self)
            %the directory where hr is located
            out = sprintf('%smrt%s',self.pathfinder(self.id,0),filesep);
        end
        function out        = path_hr(self)
            %the directory where hr is located
            out = sprintf('%smrt%sdata.nii',self.pathfinder(self.id,0),filesep);
        end
        function out        = path_epi(self,nrun,varargin)
            % simply returns the path to the mrt data. use VARARGIN to add
            % prefixes.
            if nargin == 2
                out = sprintf('%smrt%sdata.nii',self.pathfinder(self.id,nrun),filesep);
            elseif nargin == 3
                out = sprintf('%smrt%s%sdata.nii',self.pathfinder(self.id,nrun),filesep,varargin{1});
            else
                fprintf('Need to give an input...\n')
                return
            end
        end
        function out        = dir_epi(self,nrun)
            % simply returns the path to the mrt data.
            
            if nargin == 2
                out = sprintf('%smrt%s',self.pathfinder(self.id,self.dicom_target_run(nrun)),filesep);
            else
                fprintf('Need to give an input...\n')
                return
            end
        end
        function out        = path_beta(self,nrun,model_num,prefix,varargin)
            %returns the path for beta images computed in NRUN for
            %MODEL_NUM. Use VARARGIN to select a subset by indexing.
            %Actually spm_select is not even necessary here.
            
            out = self.dir_spmmat(nrun,model_num);
            out = spm_select('FPList',out,sprintf('^%sbeta_*',prefix'));
            if isempty(out)
                fprintf('No beta images found, probably wrong prefix is entered...\n');
                keyboard%sanity check
            end
            %select if VARARGIN provided
            if nargin > 4
                selector        = varargin{1};
                out             = out(selector,:);
            end
        end
        function out        = path_con(self,nrun,model_num,prefix,varargin)
            %returns the path for beta images computed in NRUN for
            %MODEL_NUM. Use VARARGIN to select a subset by indexing.
            %Actually spm_select is not even necessary here.
            
            out = self.dir_spmmat(nrun,model_num);
            out = spm_select('FPList',out,sprintf('^%scon_*',prefix'));
            if isempty(out)
                fprintf('No con images found, probably wrong prefix is entered...\n');
                %sanity check
            end
            %select if VARARGIN provided
            if nargin > 4
                selector        = varargin{1};
                out             = out(selector,:);
            end
        end
        function out        = path_FIR(self,nrun,model_num,order,nametag,varargin)
            %returns the path for FIR model
            out = self.dir_spmmat(nrun,model_num);
            out = strrep(out,'chrf',sprintf('FIR_%02d_%s',order,nametag));
            
            %select certain nifti if VARARGIN provided
            if nargin > 5
                out = spm_select('FPList',out,sprintf('%s',varargin{1}));
            end
        end
        function out        = path_contrast(self,nrun,model_num,prefix,type,varargin)
            %returns path to spm{T,F}_XXXX.nii contrast volumes in NRUN for
            %MODEL_NUM. Use PREFIX to select a subset, such s_ or s_w_.
            %TYPE selects for T or F.
            
            out = self.dir_spmmat(nrun,model_num);
            fprintf('Searching for beta images in:\n%s\n',out)
            out = spm_select('FPList',out,sprintf('^%sspm%s_*',prefix',type));
            if isempty(out)
                cprintf([1 0 0],'No SPM{F,T} images found, probably wrong prefix/run/etc is entered...\n');
                fprintf('%s\n',out)
                keyboard%sanity check
            end
            %select if VARARGIN provided
            if nargin > 5
                selector        = varargin{1};
                out             = out(selector,:);
            end
        end
        function [HRPath]   = path_hr_dicom(self)
            % finds the dicom path to the latest HR measurement for this
            % subject.
            
            HRPath = [];
            if ~ismac & ~ispc
                
                [status2 DicqOutputFull] = system(['env LD_LIBRARY_PATH= ' sprintf('/common/apps/bin/dicq --verbose  --series --exam=%s --folders',self.hr_session)]);
                %take the latest anatomical scan.
                %                     [status2 HRLine] = system(sprintf('env LD_LIBRARY_PATH= /common/apps/bin/dicq --series --exam=%s --folders -S "mprage, HR64*" ',self.hr_session));
                
                [status2 HRLine] = system(sprintf('env LD_LIBRARY_PATH= /common/apps/bin/dicq --series --exam=%s --folders | grep mprage | grep HR64 | grep prisma | tail -n 1',self.hr_session));
                %                 %
                if ~isempty(HRLine);
                    HRPath = regexp(HRLine,'/common/mrt.*/\S*','match');
                    HRPath = HRPath{1};
                    %HRPath = GetDicomPath(HRLine);
                    fprintf('Dicom Server returns:\n=====\n')
                    fprintf(DicqOutputFull);
                    %MP for errors spitted out
                    list_of_series = strsplit(DicqOutputFull,'\n');
                    for series = [list_of_series(~cellfun(@isempty, regexp(list_of_series, '^/')))]         % process only lines starting with '/'
                        warning('Error spitted out:')
                        disp(series);                   % ...or do whatever you want with 'series'
                    end
                    fprintf('=====\n');
                    fprintf('The latest recorded HR data:\n')
                    fprintf(HRLine);
                else
                    fprintf('There is no HR data found for this subject.\n Here is the output of the Dicq:\n');
                    fprintf(DicqOutputFull);
                end
            else
                fprintf('path_hr_dicom: To use dicom query you have to use one of the institute''s linux boxes\n');
            end
        end
        
        function path2data  = path_data(self,run,varargin)
            % s.path_data(4) will return the path to the subject's phase 4
            % s.path_data(4,'eye') return the path to the eye data file at the
            % 4th phase. VARARGIN{1} is a subfolder in the run folder e.g.
            % eye, mrt etc. VARARGIN{2} is file extension changer.
            
            if nargin < 2
                fprintf('you have to have at least one input for me...\n');
                return
            end
            %will return the path to phase/data_type/
            path2data = self.pathfinder(self.id , run);
            if length(varargin) >= 1
                path2data = sprintf('%s%s%sdata.mat',path2data,varargin{1},filesep);
            end
            if length(varargin) == 2
                path2data = strrep(path2data,'mat',varargin{2});
            end
        end
    end
    methods %analysis
        function [out] = get_beta_index(self,nrun,model_num,conds) % CAVE: for concatenated runs 3 and 4, gives only first index. so you need to compute the second yourself.
            load(self.path_model(nrun,model_num));
            
            out = [];
            if isnumeric(conds)
                cc = 0;
                for c = conds(:)'
                    cc = cc+1;
                    try
                        out(cc) = find(strcmp({cond.name}, num2str(c))==1);
                    catch
                        out(cc) = nan;
                    end
                end
            elseif iscell(conds)
                for cc = 1:length(conds)
                    out(cc) = find(strcmp({cond.name},conds(cc)));
                end
            elseif ischar(conds)
                out = find(strcmp({cond.name},conds));
            end
        end
        function [out] = get_Nbetas(self,nrun,model_num)
            out = size( self.path_beta(nrun,model_num,''),1);
        end
        function [out] = fit_pmf(self,varargin)
            %will load the pmf fit (saved in runXXX/pmf/pmf_fit_xxx) if computed other
            %wise will read the raw pmf data (saved in runXXX/pmf/stimulation)
            %and compute a fit.
            
            if exist(self.path_data(0,'pmf/pmf_fit')) && isempty(varargin)
                %load directly or
                load(self.path_data(0,'pmf/pmf_fit'));
                %                 fprintf('PMF Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif ~isempty(varargin) || ~exist(self.path_data(0,'pmf/pmf_fit'))
                addpath(self.path_palamedes);
                %load pmf file
                load(self.path_data(0,'pmf/stimulation'));
                
                %compute and save it.
                fprintf('Fitting PMF...\n')
                % define a search grid
                searchGrid.alpha  = linspace(0,100,10);    %structure defining grid to
                searchGrid.beta   = 10.^[-1:0.1:1];         %search for initial values
                searchGrid.gamma  = linspace(0,0.5,10);
                searchGrid.lambda = linspace(0,0.1,10);
                paramsFree        = [1 1 1 1];
                PF                = @PAL_CumulativeNormal;
                %combine 2 chains to 1 pooled chain
                p.psi.log.xrounded = cat(3,[p.psi.log.xrounded(:,:,1) nan(15,10)],[p.psi.log.xrounded(:,:,2) nan(15,10)],[p.psi.log.xrounded(:,:,1) p.psi.log.xrounded(:,:,2)]);
                
                %prepare some variables
                tchain            = size(p.psi.log.xrounded,3);
                xlevels           = unique(abs(p.psi.presentation.uniquex));
                xlevels_HD        = linspace(min(xlevels),max(xlevels),100);
                NumPos            = NaN(length(xlevels),tchain);
                OutOfNum          = NaN(length(xlevels),tchain);
                sd                = NaN(length(xlevels),tchain);
                %first collapse the two directions (pos/neg differences from
                %csp)
                
                for chain = 1:tchain
                    fprintf('Starting to fit chain %g...\n',chain)
                    %get responses, and resulting PMF from PAL algorithm
                    data = p.psi.log.xrounded(:,:,chain);
                    rep  = p.psi.presentation.rep;
                    if chain == 3
                        rep =  sum(~isnan(data),2);
                        data = sort(data,2);
                    end
                    cl   = 0;
                    for l = xlevels(:)'
                        cl                 = cl+1;
                        ind                = find(abs(p.psi.presentation.uniquex) == l);
                        collecttrials      = data(ind,1:rep(ind(1)));
                        collecttrials      = collecttrials(:);
                        NumPos(cl,chain)   = sum(collecttrials);% number of "different" responses
                        OutOfNum(cl,chain) = length(collecttrials);%number of presentations at that level
                        y_mean(cl,chain)   = NumPos(cl,chain)./OutOfNum(cl,chain);
                        y_var(cl,chain)    = (y_mean(cl,chain)*(1-y_mean(cl,chain)))./OutOfNum(cl,chain);%var of binomial distr. (np(1-p))
                        
                        sd(cl,chain)       = (OutOfNum(cl,chain)*NumPos(cl,chain)/OutOfNum(cl,chain)...
                            *(1-NumPos(cl,chain)/OutOfNum(cl,chain)))./OutOfNum(cl,chain);%var of binomial distr. (np(1-p))
                    end
                    %fit the function using PAL
                    %%
                    options             = PAL_minimize('options');
                    options.MaxIter     = 10.^3;
                    options.MaxFunEvals = 10.^3;
                    options.Display     = 'On';
                    options.ToX         = -10.^3;
                    options.TolFun      = -10.^3;
                    
                    [paramsValues LL exitflag output] = PAL_PFML_Fit(xlevels, ...
                        NumPos(:,chain), OutOfNum(:,chain), searchGrid, paramsFree, PF,'lapseLimits',[0 .5],'guessLimits',[0 .5]);
                    fprintf('%s.\n',output.message );
                    out.NumPos(chain,:)          = NumPos(:,chain);
                    out.OutOfNum(chain,:)        = OutOfNum(:,chain);
                    out.PropCorrectData(chain,:) = NumPos(:,chain)./OutOfNum(:,chain);
                    out.sd(chain,:)              = sd(:,chain);
                    out.params(chain,:)          = paramsValues;
                    out.LL(chain,:)              = LL;
                    out.exitflag(chain,:)        = exitflag;
                    out.PF                       = PF;
                    out.xlevels                  = xlevels;
                    out.y(chain,:)               = PF(paramsValues,xlevels_HD);
                    out.x(chain,:)               = xlevels_HD;
                end
                out_dir = strrep(self.path_data(0,'pmf/pmf_fit'),'data.mat','');
                if ~exist(out_dir)
                    mkdir(out_dir)
                end
                save(fullfile(out_dir,'data.mat'),'out');
            end
            self.pmf = out;
        end
        function plot_pmf(self,chains)
            load(self.path_data(0,'pmf/pmf_fit'));
            % plot the fits
            
            colors = {'r' 'c' 'k'};
            for chain = chains(:)'
                plot(out.xlevels,out.PropCorrectData(chain,:),'color',colors{chain},'linewidth',3);
                hold on;
                errorbar(out.xlevels,out.PropCorrectData(chain,:),out.sd(chain,:),'o','markersize',8,'color',colors{chain});
                plot([out.params(chain,1) out.params(chain,1)],[0 1],'color',colors{chain});
            end
            axis tight;box off;axis square;ylim([-0.1 1.2]);xlim([0 135]);drawnow;
            title(sprintf('sub:%02d (cs+:%d)',self.id,self.csp),'fontsize',12);%subject and face id
            hold on;plot(xlim,[0 0 ],'k-');plot(xlim,[0.5 0.5 ],'k:');plot(xlim,[1 1 ],'k-');hold off;%plot grid lines
        end
        function plot_pmf_fit(self,varargin)
            if nargin<2
                chains = 1:3;
            else
                chains = varargin{1};
            end
            clf;
            tchains = length(chains);
            nsp = [2 tchains+1];
            load(self.path_data(0,'pmf/pmf_fit'));
            
            colors = {'r' 'c' 'k'};
            chainname = {'CS+','CS-','merged'};
            % plot the raw data
            for chain = chains(:)'
                subplot(nsp(1),nsp(2),chain)
                plot(out.xlevels,out.PropCorrectData(chain,:),'color',colors{chain},'linewidth',3);
                hold on;
                errorbar(out.xlevels,out.PropCorrectData(chain,:),out.sd(chain,:),'o','markersize',8,'color',colors{chain});
                %                 plot([out.params(chain,1) out.params(chain,1)],[0 1],'color',colors{chain});
                axis tight;box off;axis square;ylim([-0.1 1.2]);xlim([0 135]);drawnow;
                title(chainname{chain})
                hold on;plot(xlim,[0 0 ],'k-');plot(xlim,[0.5 0.5 ],'k:');plot(xlim,[1 1 ],'k-');%plot grid lines
            end
            st = supertitle(sprintf('sub:%02d (cs+:%d)',self.id,self.csp));%subject and face id
            set(st,'FontSize',14,'Position',[0 .3 0]);
            %plot the fits
            for chain = chains(:)'
                subplot(nsp(1),nsp(2),chain+tchains+1)
                hold on;
                errorbar(out.xlevels,out.PropCorrectData(chain,:),out.sd(chain,:),'o','markersize',10,'MarkerFaceColor',colors{chain},'color',colors{chain},'LineWidth',2);
                plot([out.params(chain,1) out.params(chain,1)],[0 1],'color',colors{chain});
                plot(out.x(chain,:),out.y(chain,:),'color',colors{chain},'linewidth',3);
                axis tight;box off;axis square;ylim([-0.1 1.2]);xlim([0 135]);drawnow;
                hold on;plot(xlim,[0 0 ],'k-');plot(xlim,[0.5 0.5 ],'k:');plot(xlim,[1 1 ],'k-');hold off;%plot grid lines
            end
            subplot(nsp(1),nsp(2),tchains+1)
            for chain = chains(:)'
                plot(out.xlevels,out.PropCorrectData(chain,:),'color',colors{chain},'linewidth',3);
                hold on;
                errorbar(out.xlevels,out.PropCorrectData(chain,:),out.sd(chain,:),'o','markersize',8,'color',colors{chain});
                %                 plot([out.params(chain,1) out.params(chain,1)],[0 1],'color',colors{chain});
                axis tight;box off;axis square;ylim([-0.1 1.2]);xlim([0 135]);drawnow;
                title('all chains')
                hold on;plot(xlim,[0 0 ],'k-');plot(xlim,[0.5 0.5 ],'k:');plot(xlim,[1 1 ],'k-');%plot grid lines
            end
            subplot(nsp(1),nsp(2),nsp(1)*nsp(2))
            for chain = chains(:)'
                errorbar(out.xlevels,out.PropCorrectData(chain,:),out.sd(chain,:),'o','markersize',10,'MarkerFaceColor',colors{chain},'color',colors{chain},'LineWidth',2);
                plot([out.params(chain,1) out.params(chain,1)],[0 1],'color',colors{chain});
                plot(out.x(chain,:),out.y(chain,:),'color',colors{chain},'linewidth',3);
                axis tight;box off;axis square;ylim([-0.1 1.2]);xlim([0 135]);drawnow;
                hold on;plot(xlim,[0 0 ],'k-');plot(xlim,[0.5 0.5 ],'k:');plot(xlim,[1 1 ],'k-');%plot grid lines
            end
            set(gcf,'Color','w')
        end
        function out        = path_scr(self,varargin)
            %the directory where SCR is located
            if nargin == 1
                out = sprintf('%sscr%sdata_uncut.mat',self.pathfinder(self.id,0),filesep);
            else
                nrun = varargin{1};
                out = sprintf('%sscr%sleda_scr%sdata.mat',self.pathfinder(self.id,nrun),filesep,filesep);
            end
        end
    end
    methods %scr
        function [scr_trials_phasic,scr_trials_phasic_mean, scr_time] = get_scr_trials_phasic(self,nrun)
            force =  1;
            outputfile = strrep(strrep(self.path_scr,'run000',sprintf('run%03d',nrun)),'data_uncut.mat','scr_trials_phasic.mat');
            
            if ~exist(outputfile) || force == 1
                fprintf('Requested SCR not yet stored, computing phasic scr response for singletrials in sub %02d, run %d.\n',self.id,nrun)
                ledafile = strrep(strrep(self.path_scr,'run000/scr/',sprintf('run%03d%sscr%sleda_scr%s',nrun,filesep,filesep,filesep)),'data_uncut.mat','data_results.mat');
                if self.id == 4 && nrun == 4
                    scr_trials_phasic_mean = nan(161,10);
                    scr_time = nan(161,1);
                    scr_trials_phasic = nan(161,6,10);
                else
                    load(ledafile)
                    scr_trials_phasic = nan(size(analysis.split_driver.y,1),max(analysis.split_driver.n),10);
                    
                    scr_time = analysis.split_driver.x(:,1);
                    cc = 0;
                    for c = self.condsposition{nrun}
                        cc = cc+1;
                        ind = find(strcmp(self.condstring{c},analysis.split_driver.condnames));
                        singletrials = analysis.split_driver.y(:,ind); %single trials timecourse
                        scr_trials_phasic(:,1:length(ind),c) = singletrials;
                    end
                    scr_trials_phasic_mean = squeeze(nanmean(scr_trials_phasic,2));
                    
                end
                save(outputfile,'scr_trials_phasic','scr_time','scr_trials_phasic_mean')
                
            else
                fprintf('Requested SCR mat was stored, loading from %s.\n',outputfile);
                load(outputfile)
                fprintf('Your mat file is of size %d (timepoints) x %d (trials) x %d (conditions).\n',size(scr_trials_phasic,1),size(scr_trials_phasic,2),size(scr_trials_phasic,3));
            end
        end
        function [scr_score_trials] = get_scr_score_trials(self,nrun,varargin)
            force   = 1;
            lgt     = 1; %log transform
            
            %varargin expects timewindow;
            
            if nargin == 2
                timewindow = self.scr_timewin;
            else
                timewindow = varargin{1};
            end
            outputfile = strrep(strrep(self.path_scr,'run000',sprintf('run%03d',nrun)),'data_uncut.mat',sprintf('scr_trials_score_win%02dto%02d_log%d.mat',timewindow(1)*10,timewindow(end)*10,lgt));
            
            if ~exist(outputfile) || force == 1
                
                fprintf('Requested SCR not yet stored, computing scr responses in timewindow %g to %g in sub %02d, run %d.\n',timewindow(1), timewindow(end),self.id,nrun)
                
                dir_phasic_trials = strrep(strrep(self.path_scr,'run000',sprintf('run%03d',nrun)),'data_uncut.mat','scr_trials_phasic.mat');
                
                if ~exist(dir_phasic_trials)
                    [scr_trials_phasic, ~, scr_time] = self.get_scr_trials_phasic(nrun);
                else
                    load(dir_phasic_trials);
                end
                timey = logical((scr_time >= timewindow(1)).*(scr_time <= timewindow(end))); %index of samples in scr_time we want to include
                scr_score_trials = squeeze(nanmean(scr_trials_phasic(timey,:,:))); %single trials, time dimension is gone, mean over timewindow
                if lgt == 1
                    scr_score_trials = log10(1+abs(scr_score_trials)).*sign(scr_score_trials); %from M Benedek, C Kaernbach - Psychophysiology, 2010
                end
                save(outputfile,'scr_score_trials');
            else
                fprintf('Requested SCR mat was stored, loading from %s.\n',outputfile);
                load(outputfile);
            end
        end
        function [out_M, out_STD, out_SEM, alltrials, alltrials_z] = get_scr_score_final(self,varargin)
            force   = 1;
            lgt     = 1;
            zsc     = 1; %zscore transform
            alltrials = nan(12,30); %max 12 trials (2 runs with 6 reps, 10 conditions)
            
            
            %varargin expects timewindow;
            
            if nargin == 1
                timewindow = self.scr_timewin;
            else
                timewindow = varargin{1};
            end
            
            outputfile = strrep(self.path_data(0,'midlevel'),'data.mat',sprintf('scr_score_win%02dto%02d_log%d_z%d.mat',timewindow(1)*10,timewindow(end)*10,lgt,zsc));
            
            if ~exist(outputfile) || force == 1
                
                fprintf('Requested SCR not yet stored, computing scr responses in timewindow %g to %g in sub %02d.\n',timewindow(1), timewindow(end),self.id);
                alltrials(1:3,1:10)   = self.get_scr_score_trials(1,timewindow); % single trials in one run, log but not z yet
                alltrials(1:10,11:20) = self.get_scr_score_trials(2,timewindow); % single trials in one run, log but not z yet
                alltrials(1:6,21:30)  = self.get_scr_score_trials(3,timewindow); % single trials in one run, log but not z yet
                alltrials(7:12,21:30) = self.get_scr_score_trials(4,timewindow); % single trials in one run, log but not z yet
                alltrials_vec = alltrials(:);
                alltrials_z = reshape(nanzscore(alltrials_vec),[12 30]);
                out_M   = reshape(nanmean(alltrials_z),10,3);
                out_STD = reshape(nanstd(alltrials_z),10,3);
                out_SEM = reshape(nanstd(alltrials_z)./sqrt(sum(~isnan(alltrials_z))),10,3);
                save(strrep(outputfile,'.mat','_alltrials.mat'),'alltrials','alltrials_z','out_M','out_STD','out_SEM')
            else
                fprintf('Requested SCR mat was stored, loading from %s.\n',outputfile);
                load(outputfile);
            end
            
        end
    end
    methods %(plotters)
        function plot_log(self,nrun)
            %will plot the events that are logged during the experiment.
            L       = self.get_log(nrun);
            tevents = size(L,1);
            plot(L(1:tevents,1),L(1:tevents,2),'o','markersize',10);%plot events as dots.
            % plot lines for pulses
            hold on;
            scan_events = find(L(:,2) == 0);
            scan_times  = L(scan_events,1);
            plot([scan_times(5:5:end) scan_times(5:5:end)],ylim,'k','linewidth',.1);%plot every 5 th pulse event as a line
            % text pulse indices for each line as well.
            t_scan = length(scan_times);
            text(scan_times(5:5:t_scan),repmat(0,length(5:5:t_scan),1),num2str([5:5:t_scan]'),'color','r');
            % mark with a star missing pulses (if any)
            miss        = find(diff(scan_times) > self.TR*1.1);
            if ~isempty(miss)
                plot(scan_times(miss)+self.TR,0,'mp','markersize',30);
            end
            % text condition ids on dots.
            facelog     = 13;
            stim_events = find(L(:,2) == facelog);
            stim_types  = L(stim_events,3);
            stim_times  = L(stim_events,1);
            text(stim_times,repmat(facelog,length(stim_times),1),num2str(stim_types),'color','k');
            for dot = 1:numel(stim_events)
                colorwheel= self.GetFearGenColors;
                if stim_types(dot)<500
                    col = colorwheel(self.compute_deltacsp2ind(stim_types(dot)),:);
                elseif stim_types(dot)==500
                    col = colorwheel(end-1,:);
                elseif stim_types(dot)==3000
                    col = colorwheel(end,:);
                end
                plot(stim_times(dot),facelog,'ko','markersize',10,'MarkerFaceColor',col)
            end
            %
            hold off;
            set(gca,'ytick',[-1:18 30 31 99],'yticklabel',{'Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PainReached','FX Jump','CooledDown','TrialStart','TrialEnd','CooledDownEnd'});
            grid off;box off;
            ylim([-2 32]);
            drawnow;
        end
        function plot_keylog(self,nrun)
            p = self.get_paradigm(nrun);
            ntrials = max(p.log.ratings.relief(:,2));
            pain_confirmed = p.log.ratings.pain(1:ntrials,4);
            rel_confirmed  = p.log.ratings.relief(1:ntrials,4);
            first_face = p.out.log(find(p.out.log(:,2)==13,1,'first'),1);%marks onset of real experiment
            
            rate_pain_on  = p.out.log(p.out.log(:,2)==9,1);
            rate_pain_off = p.out.log(p.out.log(:,2)==10,1);
            rate_treat_on  = p.out.log(p.out.log(:,2)==11,1);
            rate_treat_off =  p.out.log(p.out.log(:,2)==12,1);
            thermode_off  = p.out.log(p.out.log(:,2)==18,1);
            ind_treat  = find((logical(rate_treat_on>first_face) & logical(rate_treat_on<thermode_off))); %pain ratings between first trial and before thermode off
            treat_on_oi = rate_treat_on(ind_treat);
            treat_off_oi = rate_treat_off(ind_treat);
            ind_pain = find((logical(rate_pain_on>first_face) & logical(rate_pain_on<thermode_off))); %pain ratings between first trial and before thermode off
            ind_pain = [ind_pain(1)-1; ind_pain]; %include tonic pain rating before first face
            pain_on_oi = rate_pain_on(ind_pain);
            pain_off_oi = rate_pain_off(ind_pain);
            n_ind= length(ind_pain);
            
            t_button = p.out.log(p.out.log(:,2)==7,1); %t_button = t_button(t_button>pain_on_oi(1));t_button = t_button(t_button<pain_off_oi(end)); %this also contains relief ratings, though
            id_button = p.out.log(p.out.log(:,2)==7,3); % id_button = id_button(t_button>pain_on_oi(1));id_button = id_button(t_button<pain_on_oi(end)); %this still contains relief ratings, though
            time_select = logical(t_button>pain_on_oi(1)) & logical(t_button<thermode_off);
            % bpoi = ismember(id_button,[confirm increase decrease]);
            id_button(id_button==0)=45; % for plotting
            figure
            plot(t_button(time_select),id_button(time_select),'bo')
            hold on;for n=1:length(pain_on_oi);l=line(repmat(pain_on_oi(n),2,1),ylim);set(l,'Color','g','LineStyle',':');end
            hold on;for n=1:length(pain_off_oi);l=line(repmat(pain_off_oi(n),2,1),ylim);set(l,'Color','r','LineStyle',':');end
            hold on;for n=1:length(treat_on_oi);l=line(repmat(treat_on_oi(n),2,1),ylim);set(l,'Color','g','LineStyle','-');end
            hold on;for n=1:length(treat_off_oi);l=line(repmat(treat_off_oi(n),2,1),ylim);set(l,'Color','r','LineStyle','-');end
            
            ylim([44 53])
            set(gca,'YTick',[44 49 51 52],'YTickLabel',{'code=0','increase','decrease','confirm'})
            title(sprintf('sub %02d, run %d',self.id,nrun));
        end
        function plot_motionparams(self,nrun)
            dummy = self.get_param_motion(nrun);
            subplot(2,1,1);
            plot(dummy(:,1:3));
            legend({'x','y' 'z'})
            legend boxoff;ylabel('mm');box off
            axis tight
            ylim([-5 5])
            set(gca,'ygrid','on');
            title(sprintf('sub %02d phase %02d.\n',self.id,nrun))
            %
            subplot(2,1,2)
            plot(dummy(:,4:6));
            legend({'pitch' 'roll' 'yaw'});
            legend boxoff;ylabel('degrees');box off;
            xlabel('volumes')
            axis tight
            ylim([-5 5].*10.^-2)
            set(gca,'ygrid','on')
            %             export_fig(strrep(self.path_data(nrun,'midlevel'),'data.mat',sprintf('motionparams_sub%03d_run%03d.png',self.id,nrun)))
        end
        function plot_logcomparison(self,nrun)
            %plots the data logged by the stim pc together with data logged
            %in the physio-computer. Will mark with a star the missing
            %pulses.
            clf;
            self.plot_log(nrun);
            hold on;
            L = self.get_physio2log;
            plot(L(:,1),L(:,2),'r+');
            hold off;
        end
        function plot_scr_raw(self)
            ff=figure(self.id);
            set(ff,'Position', [70 120 1600 900]);
            clf;
            load(self.path_scr);
            data = data-nanmean(data);
            plot([1:numel(data)]./1000,data)
            xlabel('time [secs]')
            hold on;
            title(sprintf('sub %02d',self.id));
            box off
            %plot scanner pulses
            firstscan = [1; find(diff(trigger(1).times)>10)+1]';
            lastscan = [firstscan(2:end)-1 trigger(1).length];
            for ss = 1:length(firstscan)
                rectangle('Position',[trigger(1).times(firstscan(ss)),min(data),trigger(1).times(lastscan(ss))-trigger(1).times(firstscan(ss)),range(data)],'EdgeColor','c');
            end
            %plot faces on curve
            plot(trigger(3).times,data(round(trigger(3).times*1000)),'r.')
            legend('SCR','FaceOnset','Scanner')
            cols = self.GetFearGenColors;
            for n = 1:7
                plot(trigger(n).times,-.5*n,'ko','MarkerFaceColor',cols(n,:))
                ytag{n} = trigger(n).name;
            end
            set(gca,'YTick',-4:.5:0,'YTickLabel',ytag)
        end
        function plot_scr_phases(self)
            ff=figure(1000);
            set(ff,'Position', [32 473 1302 535]);
            clf;
            load(self.path_scr);
            subplot(2,4,1:4);
            plot([1:numel(data)]./1000,data)
            xlabel('time [secs]')
            hold on;
            title(sprintf('sub %02d',self.id));
            box off
            %plot scanner pulses
            firstscan = [1; find(diff(trigger(1).times)>10)+1]';
            lastscan = [firstscan(2:end)-1 trigger(1).length];
            for ss = 1:length(firstscan)
                rectangle('Position',[trigger(1).times(firstscan(ss)),min(data),trigger(1).times(lastscan(ss))-trigger(1).times(firstscan(ss)),range(data)],'EdgeColor','c');
            end
            plot(trigger(3).times,data(round(trigger(3).times*1000)),'r.')
            %              plot(trigger(1).times,max(data)+.3,'c*');
            legend('SCR','FaceOnset','Scanner')
            for n = 1:4
                load(self.path_scr(n));
                subplot(2,4,4+n);
                plot(data.time,data.conductance);
                hold on;
                col = self.GetFearGenColors;
                for c = 1:length(data.event)
                    plot(data.event(c).time,data.conductance(round(data.event(c).time*1000)),'o','Color',col(data.event(c).nid,:),'MarkerFaceColor',col(data.event(c).nid,:))
                    %                     fprintf('Plot at %04.2f , %04.2f, col no %d\n',data.event(c).time,data.conductance(round(data.event(c).time*1000)),data.event(c).nid); %just a check
                end
                box off
                clear data;
            end
            EqualizeSubPlotYlim(gcf);
            savepath = strrep(self.path_scr,'data_uncut.mat','scr_plot.png');
            if ~exist(savepath)
                SaveFigure(savepath);
            end
        end
    end
    methods %(fmri analysis)
        function [Ndummy] = KickDummies(self)
            kick  = 0;
            runs = 1:self.total_run-1;
            rc = 0;
            for run = runs;
                rc = rc+1;
                L = self.get_log(run);
                scantimes = L(L(:,2)==0,1);
                FixOnsets  = L(L(:,2)== 2,1);
                FirstOnset = FixOnsets(1);
                Ndummy(rc) = find(scantimes > FirstOnset,1,'first')-2; %when pulse is after FixCross, it's the second pulse already (of those we want to keep).
                %                 fprintf('%02d dummies before fix onset.',Ndummy)
                if kick == 1
                    dummyfolder = [self.path_data(run) filesep 'mri' filesep 'dummyscans' filesep];
                    if ~exist(dummyfolder)
                        mkdir(dummyfolder)
                        Nscans = length(dir([self.path_data(run) filesep 'mri']))<300;
                        if length(dir([self.path_data(run) filesep 'mri']))<300 % then it's already converted or sth
                            
                            fprintf('Not many files (%03d) found, please check procedure.')
                            keyboard;
                        end
                    else
                        fprintf('Dummyfolder already found! Check how to proceed!')
                        keyboard;
                    end
                end
            end
            
            if length(Ndummy) < 4
                Ndummy(end+1) = NaN;
            end
        end
        function CompareTriallist2Log(self,run,modelnum)
            
            
        end
        function PrepareCondFile(self,run,modelnum)
            
            force = 1;
            verbalize = 0;
            modelpath = self.path_model(run,modelnum);
            
            switch modelnum
                case 1
                    modelname = 'FaceOnsets'; %tonic pain is baseline, everything else is modelled
                    if ~exist(modelpath) || force == 1
                        if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                        if verbalize == 1
                            fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                        end
                        L               = self.get_log(run);
                        %sort things according to time rather than order of being logged
                        [~,i]           = sort(L(:,1),'ascend');
                        L               = L(i,:);
                        % delete all the events that are after the last scanning..
                        scan_times      = L(find(L(:,2) == 0),1);
                        first_scan_time = min(scan_times);
                        last_scan_time  = max(scan_times);
                        L(L(:,1) < first_scan_time,:) = [];
                        L(L(:,1) > last_scan_time,:)  = [];
                        L(:,1)                        = L(:,1) - first_scan_time; % correct time logging so that it's referring to the first scan as 0. (later first col will be from > dummy trials)
                        
                        % correct timings for dummy scans
                        L(1:6,:) = []; %let dummies go
                        first_real_scan = L(1,1);
                        L(:,end+1) = L(:,1); %store it so we still know the original timing
                        L(:,1) = L(:,1)-first_real_scan;
                        
                        %event types are as follows:
                        %         %Pulse Detection      :     0    info: NaN;
                        %         %Tracker Onset        :     1
                        %         %Cross (tonic) Onset  :     2    info: position
                        %         %Cross (pain) Onset   :     3    info: position
                        %         %Ramp down Onset      :     4    info: ror
                        %         %Treatment Plateau    :     5    info: temp
                        %         %Ramp back onset      :     6    info: ror;
                        %         %Key Presses          :     7    info: keycode;
                        %         %Tracker Offset       :     8    info: NaN;
                        %         %Rate pain Onset		:     9    info: nTrial;
                        %         %Rate pain Offset     :     10   info: nTrial;
                        %         %Rate treat Onset     :     11   info: nTrial;
                        %         %Rate treat Offset    :     12   info: nTrial;
                        %         %Face Onset           :     13   info: dist;
                        %         %Face Offset          :     14   info: dist;
                        %         %FaceStim Fixcross    :     15   info: position(1)
                        %         %Tonic Pain reached   :     16   info: nTrial
                        %         %Fixcross Jump        :     17   info:
                        %         %dummy fixflip        :     22   info: NaN;
                        %         planned trialstart    :     30   info: NaN
                        %         planned trialend      :     31   info: NaN
                        
                        
                        names = {'VAS','Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PlateauReached','FixJump','CoolDown','TrialStart','TrialEnd'};
                        condnum = [-2 -1 0 1:18 30 31];
                        
                        TR = Project.TR;
                        scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
                        scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
                        
                        %% model 1 - relief
                        
                        conds = [13 11 9 18];% 13 = face, 11 = RateTreat, 9=RatePain
                        if verbalize ==1
                            fprintf('Considering the following conditions:\n')
                            for c = conds
                                fprintf('%s\n',names{condnum==c})
                            end
                        end
                        
                        
                        all_events          = find(ismember(L(:,2),conds));
                        LL = L(all_events,:);
                        
                        if verbalize ==1
                            fprintf('Overall we have %d events.\n',length(all_events))
                        end
                        
                        %Conds
                        
                        %1 -135
                        %2 -90
                        %3 -45
                        %4  CSP
                        %5 45
                        %6 90
                        %7 135
                        %8 180
                        %9 UCS
                        %10 Null
                        %11 RatePain
                        %12 RateRelief
                        %13 BaselineEnd
                        
                        FaceOnsets       = LL(LL(:,2)==13,1);
                        RateReliefOnsets = LL(LL(:,2)==11,1);
                        RatePainOnsets   = LL(LL(:,2)==9,1);
                        
                        
                        conds_presented = unique(LL(LL(:,2)==13,3));
                        cond = [];
                        cc = 0;
                        for c = conds_presented(:)';
                            cc = cc+1;
                            trial_ind = find(LL(LL(:,2)==13,3)==c);
                            
                            cond(cc).name     = mat2str(c);
                            cond(cc).onset    = FaceOnsets(trial_ind);
                            cond(cc).duration = RateReliefOnsets(trial_ind)-FaceOnsets(trial_ind);
                            if verbalize ==1
                                fprintf('Number of onsets for cond %04g: %d.\n',c,length(cond(cc).onset))
                            end
                        end
                        
                        % Rate Pain
                        cc = cc+1;
                        cond(cc).name   = 'RatePain';
                        cond(cc).onset  = RatePainOnsets;
                        cond(cc).duration = L(L(:,2)==10,1)-RatePainOnsets;
                        if verbalize ==1
                            fprintf('Number of onsets for RatePain: %d.\n',length(cond(cc).onset))
                        end
                        if self.id == 6
                            corrdur = [1.5 2 2 1]; %scanner stopped too early, during last rating.
                            cond(cc).duration(end) = corrdur(run);
                        elseif self.id == 32 && run == 1
                            cond(cc).onset(end-1:end) = [];
                            cond(cc).duration(end-1:end) = [];
                        end
                        
                        % Rate Relief
                        cc = cc+1;
                        cond(cc).name     = 'RateRelief';
                        cond(cc).onset    = RateReliefOnsets;
                        cond(cc).duration = ones(1,length(cond(cc).onset)).*5;
                        if self.id == 32 && run ==1
                            cond(cc).duration(end) = 4; % scanner stopped here
                        end
                        if verbalize ==1
                            fprintf('Number of onsets for RateRelief: %d.\n',length(cond(cc).onset))
                        end
                        
                        % CoolDown phase
                        cc = cc+1;
                        RampDown        = L(L(:,2)==4,1);
                        RampDownEnd     = RampDown(end);
                        cond(cc).name   = '999';
                        cond(cc).onset  = RampDownEnd;
                        cond(cc).duration = RatePainOnsets(end)-L(find(L(:,2)==4,1,'last'),1);
                        if self.id == 4
                            cond(cc).duration = 32.6; %logging problem.
                        elseif self.id == 32 && run == 1 %wrong seq, scanner stopped already
                            cond(cc) = [];
                        end
                        
                        if verbalize ==1
                            fprintf('Number of onsets for CoolDown %d.\n',length(cond(end).onset))
                        end
                        
                        
                        %% housekeeping
                        for counter = 1:length(cond)
                            cond(counter).onset    = cond(counter).onset./TR;
                            cond(counter).duration = cond(counter).duration./TR;
                            cond(counter).tmod      = 0;
                            cond(counter).pmod      = struct('name',{},'param',{},'poly',{});
                        end
                        
                        save(modelpath,'cond');
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 2
                    modelname = 'FaceAndRampOnset'; %tonic pain is baseline, everything else is modelled
                    
                    if ~exist(modelpath) || force == 1
                        if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                        if verbalize == 1
                            fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                        end
                        L               = self.get_log(run);
                        %sort things according to time rather than order of being logged
                        [~,i]           = sort(L(:,1),'ascend');
                        L               = L(i,:);
                        % delete all the events that are after the last scanning..
                        scan_times      = L(find(L(:,2) == 0),1);
                        first_scan_time = min(scan_times);
                        last_scan_time  = max(scan_times);
                        L(L(:,1) < first_scan_time,:) = [];
                        L(L(:,1) > last_scan_time,:)  = [];
                        L(:,1)                        = L(:,1) - first_scan_time; % correct time logging so that it's referring to the first scan as 0. (later first col will be from > dummy trials)
                        
                        % correct timings for dummy scans
                        L(1:6,:) = []; %let dummies go
                        first_real_scan = L(1,1);
                        L(:,end+1) = L(:,1); %store it so we still know the original timing
                        L(:,1) = L(:,1)-first_real_scan;
                        
                        %event types are as follows:
                        %         %Pulse Detection      :     0    info: NaN;
                        %         %Tracker Onset        :     1
                        %         %Cross (tonic) Onset  :     2    info: position
                        %         %Cross (pain) Onset   :     3    info: position
                        %         %Ramp down Onset      :     4    info: ror
                        %         %Treatment Plateau    :     5    info: temp
                        %         %Ramp back onset      :     6    info: ror;
                        %         %Key Presses          :     7    info: keycode;
                        %         %Tracker Offset       :     8    info: NaN;
                        %         %Rate pain Onset		:     9    info: nTrial;
                        %         %Rate pain Offset     :     10   info: nTrial;
                        %         %Rate treat Onset     :     11   info: nTrial;
                        %         %Rate treat Offset    :     12   info: nTrial;
                        %         %Face Onset           :     13   info: dist;
                        %         %Face Offset          :     14   info: dist;
                        %         %FaceStim Fixcross    :     15   info: position(1)
                        %         %Tonic Pain reached   :     16   info: nTrial
                        %         %Fixcross Jump        :     17   info:
                        %         %dummy fixflip        :     22   info: NaN;
                        %         planned trialstart    :     30   info: NaN
                        %         planned trialend      :     31   info: NaN
                        
                        
                        names = {'VAS','Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PlateauReached','FixJump','CoolDown','TrialStart','TrialEnd'};
                        condnum = [-2 -1 0 1:18 30 31];
                        
                        TR = Project.TR;
                        scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
                        scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
                        
                        %% model 2 - FaceAndRampOnset
                        
                        conds = [13 4 11 9 18];% 4 = RampOnset 11 = RateTreat, 9=RatePain 18= CoolDown.
                        if verbalize ==1
                            fprintf('Considering the following conditions:\n')
                            for c = conds
                                fprintf('%s\n',names{condnum==c})
                            end
                        end
                        
                        
                        all_events          = find(ismember(L(:,2),conds));
                        LL = L(all_events,:);
                        
                        if verbalize ==1
                            fprintf('Overall we have %d events.\n',length(all_events))
                        end
                        
                        %Conds
                        % all face conds with Face Onsets
                        %1 -135
                        %2 -90
                        %3 -45
                        %4  CSP
                        %5 45
                        %6 90
                        %7 135
                        %8 180
                        %9 UCS
                        %10 Null
                        % repeat the same conditions for Ramp Onset
                        %11 RatePain
                        %12 RateRelief
                        %13 BaselineEnd
                        
                        FaceOnsets          = LL(LL(:,2)==13,1);
                        RampDownOnsets      = LL(LL(:,2)==4,1);
                        RampDownOnsets(end) = []; % this is the CoolDown at the very end, no normal ramp down.
                        RateReliefOnsets    = LL(LL(:,2)==11,1);
                        RatePainOnsets      = LL(LL(:,2)==9,1);
                        
                        
                        conds_presented = unique(LL(LL(:,2)==13,3));
                        Ncols = length(conds_presented);
                        if verbalize == 1
                            fprintf('In phase: %s, we have %02d different face conds.\n',Project.plottitles{run},Ncols)
                        end
                        cond = [];
                        cc = 0;
                        for c = conds_presented(:)';
                            cc = cc+1;
                            trial_ind = find(LL(LL(:,2)==13,3)==c);
                            % RampOnset for this Cond
                            cond(cc).name     = mat2str(c);
                            cond(cc).onset    = RampDownOnsets(trial_ind);
                            cond(cc).duration = RateReliefOnsets(trial_ind)-RampDownOnsets(trial_ind);
                            % face for this Cond
                            cond(cc+Ncols).name     = [mat2str(c) 'Face'];
                            cond(cc+Ncols).onset    = FaceOnsets(trial_ind);
                            cond(cc+Ncols).duration = repmat(1.5,length(FaceOnsets(trial_ind)),1); % including the little jitter:  RampDownOnsets(trial_ind)-FaceOnsets(trial_ind);
                            
                            if verbalize ==1
                                fprintf('Number of onsets for cond %04g: %d.\n',c,length(cond(cc).onset));
                                fprintf('Mean duration: %04.2f seconds.\n',mean(cond(cc).duration));
                            end
                        end
                        cc = length(cond); %need to get to the end of the face/ramp columns.
                        
                        % Rate Pain
                        cc = cc+1;
                        cond(cc).name   = 'RatePain';
                        cond(cc).onset  = RatePainOnsets;
                        cond(cc).duration = L(L(:,2)==10,1)-RatePainOnsets;
                        if verbalize ==1
                            fprintf('Number of onsets for RatePain: %d.\n',length(cond(cc).onset))
                            fprintf('Mean duration: %04.2f seconds.\n',mean(cond(cc).duration));
                        end
                        if self.id == 6
                            corrdur = [1.5 2 2 1]; %scanner stopped too early, during last Pain rating.
                            cond(cc).duration(end) = corrdur(run);
                        elseif self.id == 32 && run == 1
                            cond(cc).onset(end-1:end) = [];
                            cond(cc).duration(end-1:end) = [];
                        end
                        
                        % Rate Relief
                        cc = cc+1;
                        cond(cc).name     = 'RateRelief';
                        cond(cc).onset    = RateReliefOnsets;
                        cond(cc).duration = ones(1,length(cond(cc).onset)).*5;
                        if self.id == 32 && run ==1
                            cond(cc).duration(end) = 4; % scanner stopped here
                        end
                        if verbalize ==1
                            fprintf('Number of onsets for RateRelief: %d.\n',length(cond(cc).onset))
                            fprintf('Mean duration: %04.2f seconds.\n',mean(cond(cc).duration));
                        end
                        
                        % CoolDown phase
                        cc = cc+1;
                        RampDown        =L(L(:,2)==4,1);
                        RampDownEnd     = RampDown(end);
                        cond(cc).name   = '999';
                        cond(cc).onset  = RampDownEnd;
                        cond(cc).duration = RatePainOnsets(end)-L(find(L(:,2)==4,1,'last'),1);
                        if self.id == 4
                            cond(cc).duration = 32.6; %logging problem.
                        elseif self.id == 32 && run == 1 %wrong seq, scanner stopped already
                            cond(cc) = [];
                        end
                        
                        if verbalize ==1
                            fprintf('Number of onsets for CoolDown %d.\n',length(cond(end).onset))
                            fprintf('Mean duration: %04.2f seconds.\n',mean(cond(cc).duration));
                        end
                        
                        %% housekeeping
                        for counter = 1:length(cond)
                            cond(counter).onset    = cond(counter).onset./TR;
                            cond(counter).duration = cond(counter).duration./TR;
                            cond(counter).tmod      = 0;
                            cond(counter).pmod      = struct('name',{},'param',{},'poly',{});
                        end
                        save(modelpath,'cond');
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 3
                    modelname = 'RampOn'; %tonic pain is baseline, everything else is modelled, except faces. they go into ISI.
                    
                    if ~exist(modelpath) || force == 1
                        if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                        if verbalize == 1
                            fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                        end
                        L               = self.get_log(run);
                        %sort things according to time rather than order of being logged
                        [~,i]           = sort(L(:,1),'ascend');
                        L               = L(i,:);
                        % delete all the events that are after the last scanning..
                        scan_times      = L(find(L(:,2) == 0),1);
                        first_scan_time = min(scan_times);
                        last_scan_time  = max(scan_times);
                        L(L(:,1) < first_scan_time,:) = [];
                        L(L(:,1) > last_scan_time,:)  = [];
                        L(:,1)                        = L(:,1) - first_scan_time; % correct time logging so that it's referring to the first scan as 0. (later first col will be from > dummy trials)
                        
                        % correct timings for dummy scans
                        L(1:6,:) = []; %let dummies go
                        first_real_scan = L(1,1);
                        L(:,end+1) = L(:,1); %store it so we still know the original timing
                        L(:,1) = L(:,1)-first_real_scan;
                        
                        %event types are as follows:
                        %         %Pulse Detection      :     0    info: NaN;
                        %         %Tracker Onset        :     1
                        %         %Cross (tonic) Onset  :     2    info: position
                        %         %Cross (pain) Onset   :     3    info: position
                        %         %Ramp down Onset      :     4    info: ror
                        %         %Treatment Plateau    :     5    info: temp
                        %         %Ramp back onset      :     6    info: ror;
                        %         %Key Presses          :     7    info: keycode;
                        %         %Tracker Offset       :     8    info: NaN;
                        %         %Rate pain Onset		:     9    info: nTrial;
                        %         %Rate pain Offset     :     10   info: nTrial;
                        %         %Rate treat Onset     :     11   info: nTrial;
                        %         %Rate treat Offset    :     12   info: nTrial;
                        %         %Face Onset           :     13   info: dist;
                        %         %Face Offset          :     14   info: dist;
                        %         %FaceStim Fixcross    :     15   info: position(1)
                        %         %Tonic Pain reached   :     16   info: nTrial
                        %         %Fixcross Jump        :     17   info:
                        %         %Cooldown end         :     18   info:
                        %         %dummy fixflip        :     22   info: NaN;
                        %         planned trialstart    :     30   info: NaN
                        %         planned trialend      :     31   info: NaN
                        
                        
                        names = {'VAS','Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PlateauReached','FixJump','CoolDown','TrialStart','TrialEnd'};
                        condnum = [-2 -1 0 1:18 30 31];
                        
                        TR = Project.TR;
                        scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
                        scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
                        
                        %% get onsets
                        
                        conds = [4 11 9 18];% (13 = face is gone) 4 = RampDown, 11 = RateTreat, 9=RatePain 18 = CoolDown at very end.
                        if verbalize ==1
                            fprintf('Considering the following conditions:\n')
                            for c = conds
                                fprintf('%s\n',names{condnum==c})
                            end
                        end
                        
                        
                        all_events          = find(ismember(L(:,2),conds));
                        LL = L(all_events,:);
                        
                        if verbalize ==1
                            fprintf('Overall we have %d events.\n',length(all_events))
                        end
                        
                        %Conds
                        
                        %1 -135
                        %2 -90
                        %3 -45
                        %4  CSP
                        %5 45
                        %6 90
                        %7 135
                        %8 180
                        %9 UCS
                        %10 Null
                        %11 RatePain
                        %12 RateRelief
                        %13 BaselineEnd
                        
                        RampdownOnsets   = LL(LL(:,2)==4,1);
                        RateReliefOnsets = LL(LL(:,2)==11,1);
                        RatePainOnsets   = LL(LL(:,2)==9,1);
                        
                        
                        conds_presented = unique(L(L(:,2)==13,3));
                        cond = [];
                        cc = 0;
                        for c = conds_presented(:)';
                            cc = cc+1;
                            trial_ind = find(L(L(:,2)==13,3)==c);
                            
                            cond(cc).name     = mat2str(c);
                            cond(cc).onset    = RampdownOnsets(trial_ind);
                            cond(cc).duration = RateReliefOnsets(trial_ind)-RampdownOnsets(trial_ind);
                            if verbalize ==1
                                fprintf('Number of onsets for cond %04g: %d.\n',c,length(cond(cc).onset))
                            end
                        end
                        
                        % Rate Pain
                        cc = cc+1;
                        cond(cc).name   = 'RatePain';
                        cond(cc).onset  = RatePainOnsets;
                        cond(cc).duration = L(L(:,2)==10,1)-RatePainOnsets;
                        if verbalize ==1trial_ind
                            fprintf('Number of onsets for RatePain: %d.\n',length(cond(cc).onset))
                        end
                        if self.id == 6
                            corrdur = [1.5 2 2 1]; %scanner stopped too early, during last rating.
                            cond(cc).duration(end) = corrdur(run);
                        elseif self.id == 32 && run == 1
                            cond(cc).onset(end-1:end) = [];
                            cond(cc).duration(end-1:end) = [];
                        end
                        
                        % Rate Relief
                        cc = cc+1;
                        cond(cc).name     = 'RateRelief';
                        cond(cc).onset    = RateReliefOnsets;
                        cond(cc).duration = ones(1,length(cond(cc).onset)).*5;
                        if self.id == 32 && run ==1
                            cond(cc).duration(end) = 4; % scanner stopped here
                        end
                        if verbalize ==1
                            fprintf('Number of onsets for RateRelief: %d.\n',length(cond(cc).onset))
                        end
                        
                        % CoolDown phase
                        cc = cc+1;
                        RampDown        = L(L(:,2)==4,1);
                        RampDownEnd     = RampDown(end);Prepare
                        cond(cc).name   = '999';
                        cond(cc).onset  = RampDownEnd;
                        cond(cc).duration = RatePainOnsets(end)-L(find(L(:,2)==4,1,'last'),1);
                        if self.id == 4
                            cond(cc).duration = 32.6; %logging problem.
                        elseif self.id == 32 && run == 1 %wrong seq, scanner stopped already
                            cond(cc) = [];
                        end
                        
                        if verbalize ==1
                            fprintf('Number of onsets for CoolDown %d.\n',length(cond(end).onset))
                        end
                        
                        %% housekeeping
                        for counter = 1:length(cond)
                            cond(counter).onset    = cond(counter).onset./TR;
                            cond(counter).duration = cond(counter).duration./TR;
                            cond(counter).tmod      = 0;
                            cond(counter).pmod      = struct('name',{},'param',{},'poly',{});
                        end
                        
                        save(modelpath,'cond');
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                    
                case 4
                    modelname = 'RampOnsetStick'; %based on model_03, onsets are RampOnsets, but we model RampOnsets as stick (but everything else as box)
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,3))
                        set2zero = self.nreliefconds(run);
                        %% set durations to zero
                        for counter = 1:set2zero
                            cond(counter).duration = zeros(size(cond(counter).onset));
                        end
                        save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 5
                    modelname = 'RampOnsetBox_NoCool'; %based on model_03, onsets are RampOnsets, but we model RampOnsets as box
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,3))
                        
                        if self.id ==32 && run ==1
                            fprintf('Stays like model 4.\n')
                            save(modelpath,'cond')
                        else
                            fprintf('Cutting conds.\n')
                            %get rid of CoolDown
                            CoolCond = self.nreliefconds(run)+3;
                            cond(CoolCond) = [];
                            %get rid of PainRating at the end
                            RatePainCond = self.nreliefconds(run)+1;
                            cond(RatePainCond).onset(end) = [];
                            cond(RatePainCond).duration(end) = [];
                            save(modelpath,'cond')
                        end
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 6
                    modelname = 'RampOnsetStick_NoCool'; %based on model_03, onsets are RampOnsets, but we model RampOnsets as stick (but everything else as box)
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,4))
                        set2zero = self.nreliefconds(run);
                        %% set durations to zero
                        for counter = 1:set2zero
                            cond(counter).duration = zeros(size(cond(counter).onset));
                        end
                        if self.id ==32 && run ==1
                            fprintf('Stays like model 4.\n')
                            save(modelpath,'cond')
                        else
                            fprintf('Cutting conds.\n')
                            %get rid of Cooldown
                            CoolCond = self.nreliefconds(run)+3;
                            cond(CoolCond) = [];
                            %get rid of PainRating at the end
                            RatePainCond = self.nreliefconds(run)+1;
                            cond(RatePainCond).onset(end) = [];
                            cond(RatePainCond).duration(end) = [];
                            save(modelpath,'cond')
                        end
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 7
                    modelname = 'RampOnsetStick_8conds1regr_NoCool'; %based on model_03, onsets are RampOnsets, but we model RampOnsets as sticks, in one vector
                    
                    if ~exist(modelpath) || force == 1
                        if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                        if verbalize == 1
                            fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                        end
                        L               = self.get_log(run);
                        %sort things according to time rather than order of being logged
                        [~,i]           = sort(L(:,1),'ascend');
                        L               = L(i,:);
                        % delete all the events that are after the last scanning..
                        scan_times      = L(find(L(:,2) == 0),1);
                        first_scan_time = min(scan_times);
                        last_scan_time  = max(scan_times);
                        L(L(:,1) < first_scan_time,:) = [];
                        L(L(:,1) > last_scan_time,:)  = [];
                        L(:,1)                        = L(:,1) - first_scan_time; % correct time logging so that it's referring to the first scan as 0. (later first col will be from > dummy trials)
                        
                        % correct timings for dummy scans
                        L(1:6,:) = []; %let dummies go
                        first_real_scan = L(1,1);
                        L(:,end+1) = L(:,1); %store it so we still know the original timing
                        L(:,1) = L(:,1)-first_real_scan;
                        
                        %event types are as follows:
                        %         %Pulse Detection      :     0    info: NaN;
                        %         %Tracker Onset        :     1
                        %         %Cross (tonic) Onset  :     2    info: position
                        %         %Cross (pain) Onset   :     3    info: position
                        %         %Ramp down Onset      :     4    info: ror
                        %         %Treatment Plateau    :     5    info: temp
                        %         %Ramp back onset      :     6    info: ror;
                        %         %Key Presses          :     7    info: keycode;
                        %         %Tracker Offset       :     8    info: NaN;
                        %         %Rate pain Onset		:     9    info: nTrial;
                        %         %Rate pain Offset     :     10   info: nTrial;
                        %         %Rate treat Onset     :     11   info: nTrial;
                        %         %Rate treat Offset    :     12   info: nTrial;
                        %         %Face Onset           :     13   info: dist;
                        %         %Face Offset          :     14   info: dist;
                        %         %FaceStim Fixcross    :     15   info: position(1)
                        %         %Tonic Pain reached   :     16   info: nTrial
                        %         %Fixcross Jump        :     17   info:
                        %         %Cooldown end         :     18   info:
                        %         %dummy fixflip        :     22   info: NaN;
                        %         planned trialstart    :     30   info: NaN
                        %         planned trialend      :     31   info: NaN
                        
                        
                        names = {'VAS','Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PlateauReached','FixJump','CoolDown','TrialStart','TrialEnd'};
                        condnum = [-2 -1 0 1:18 30 31];
                        
                        TR = Project.TR;
                        scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
                        scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
                        
                        %% get onsets
                        
                        conds = [4 11 9 18];% (13 = face is gone) 4 = RampDown, 11 = RateTreat, 9=RatePain 18 = CoolDown at very end.
                        if verbalize ==1
                            fprintf('Considering the following conditions:\n')
                            for c = conds
                                fprintf('%s\n',names{condnum==c})
                            end
                        end
                        
                        
                        all_events          = find(ismember(L(:,2),conds));
                        LL = L(all_events,:);
                        
                        if verbalize ==1
                            fprintf('Overall we have %d events.\n',length(all_events))
                        end
                        
                        %Conds
                        
                        %1 -135
                        %2 -90
                        %3 -45
                        %4  CSP
                        %5 45
                        %6 90
                        %7 135
                        %8 180
                        %9 UCS
                        %10 Null
                        %11 RatePain
                        %12 RateRelief
                        %13 BaselineEnd
                        
                        RampdownOnsets   = LL(LL(:,2)==4,1);
                        RateReliefOnsets = LL(LL(:,2)==11,1);
                        RatePainOnsets   = LL(LL(:,2)==9,1);
                        RampdownFinal    = RampdownOnsets(end);
                        RampdownOnsets(end) = [];
                        
                        cc = 0;
                        cond = [];
                        cc = cc+1;
                        cond_list = L(L(:,2)==13,3);
                        
                        cond(cc).name     = 'RampDown';
                        cond(cc).onset    = RampdownOnsets(cond_list < 500);
                        cond(cc).duration = 0;
                        if verbalize ==1
                            fprintf('Number of onsets for %s: %d.\n',cond(1).name,length(cond(1).onset))
                        end
                        cc = cc+1;
                        cond(cc).name     = 'Null';
                        cond(cc).onset    = RampdownOnsets(cond_list ==3000);
                        cond(cc).duration = 0;
                        
                        if run > 1
                            cc = cc+1;
                            cond(cc).name     = 'UCS';
                            cond(cc).onset    = RampdownOnsets(cond_list == 500);
                            cond(cc).duration = 0;
                        end
                        
                        % Rate Pain
                        cc = cc+1;
                        cond(cc).name   = 'RatePain';
                        cond(cc).onset  = RatePainOnsets;
                        cond(cc).duration = L(L(:,2)==10,1)-RatePainOnsets;
                        if verbalize ==1
                            fprintf('Number of onsets for RatePain: %d.\n',length(cond(cc).onset))
                        end
                        if self.id == 6
                            corrdur = [1.5 2 2 1]; %scanner stopped too early, during last rating.
                            cond(cc).duration(end) = corrdur(run);
                        elseif self.id == 32 && run == 1
                            cond(cc).onset(end-1:end) = [];
                            cond(cc).duration(end-1:end) = [];
                        end
                        
                        % Rate Relief
                        cc = cc+1;
                        cond(cc).name     = 'RateRelief';
                        cond(cc).onset    = RateReliefOnsets;
                        cond(cc).duration = ones(1,length(cond(cc).onset)).*5;
                        if self.id == 32 && run ==1
                            cond(cc).duration(end) = 4; % scanner stopped here
                        end
                        if verbalize ==1
                            fprintf('Number of onsets for RateRelief: %d.\n',length(cond(cc).onset))
                        end
                        
                        %                         % CoolDown phase
                        %                         cc = cc+1;
                        %                         RampDown        = L(L(:,2)==4,1);
                        %                         RampDownEnd     = RampDown(end);
                        %                         cond(cc).name   = '999';
                        %                         cond(cc).onset  = RampDownEnd;
                        %                         cond(cc).duration = RatePainOnsets(end)-L(find(L(:,2)==4,1,'last'),1);
                        %                         if self.id == 4
                        %                             cond(cc).duration = 32.6; %logging problem.
                        %                         elseif self.id == 32 && run == 1 %wrong seq, scanner stopped already
                        %                             cond(cc) = [];
                        %                         end
                        %
                        %                         if verbalize ==1
                        %                             fprintf('Number of onsets for CoolDown %d.\n',length(cond(end).onset))
                        %                         end
                        
                        
                        %% housekeeping
                        for counter = 1:length(cond)
                            cond(counter).onset    = cond(counter).onset./TR;
                            cond(counter).duration = cond(counter).duration./TR;
                            cond(counter).tmod      = 0;
                            cond(counter).pmod      = struct('name',{},'param',{},'poly',{});
                        end
                        
                        %% Tuning information as pmod, VonMises and derivative dVM/dkappa
                        conds = -135:45:180;
                        
                        amp   = 1;
                        kappa = 1;
                        delta = .01;
                        VMlookup = zscore(Tuning.VonMises(conds,amp,kappa,0,0));
                        dVMlookup = -zscore((Tuning.VonMises(conds,amp,kappa+delta,0,0)-Tuning.VonMises(conds,amp,kappa-delta,0,0))./(2*delta)); %central difference formula
                        %
                        condlist = cond_list(cond_list < 500);
                        for ntrial = 1:length(cond(1).onset)
                            deltacsp = condlist(ntrial);
                            ind = Project.compute_deltacsp2ind(deltacsp);
                            VM(ntrial)    = VMlookup(ind);
                            dVM(ntrial)   = dVMlookup(ind);
                        end
                        
                        % linear regression y ~ 1 + vM + dVM
                        % Interaction is done on 2ndlevel
                        cond(1).pmod  = struct('name',{'vM','dvM'},'param',{VM dVM},'poly',{1 1});
                        
                        save(modelpath,'cond');
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 8
                    modelname = 'RampOnsetStick_allconds_interactVMdVM_NoCool'; %based on model_03, onsets are RampOnsets, but we model RampOnsets as stick, in one vector
                    
                    if ~exist(modelpath) || force == 1
                        if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                        if verbalize == 1
                            fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                        end
                        L               = self.get_log(run);
                        %sort things according to time rather than order of being logged
                        [~,i]           = sort(L(:,1),'ascend');
                        L               = L(i,:);
                        % delete all the events that are after the last scanning..
                        scan_times      = L(find(L(:,2) == 0),1);
                        first_scan_time = min(scan_times);
                        last_scan_time  = max(scan_times);
                        L(L(:,1) < first_scan_time,:) = [];
                        L(L(:,1) > last_scan_time,:)  = [];
                        L(:,1)                        = L(:,1) - first_scan_time; % correct time logging so that it's referring to the first scan as 0. (later first col will be from > dummy trials)
                        
                        % correct timings for dummy scans
                        L(1:6,:) = []; %let dummies go
                        first_real_scan = L(1,1);
                        L(:,end+1) = L(:,1); %store it so we still know the original timing
                        L(:,1) = L(:,1)-first_real_scan;
                        
                        %event types are as follows:
                        %         %Pulse Detection      :     0    info: NaN;
                        %         %Tracker Onset        :     1
                        %         %Cross (tonic) Onset  :     2    info: position
                        %         %Cross (pain) Onset   :     3    info: position
                        %         %Ramp down Onset      :     4    info: ror
                        %         %Treatment Plateau    :     5    info: temp
                        %         %Ramp back onset      :     6    info: ror;
                        %         %Key Presses          :     7    info: keycode;
                        %         %Tracker Offset       :     8    info: NaN;
                        %         %Rate pain Onset		:     9    info: nTrial;
                        %         %Rate pain Offset     :     10   info: nTrial;
                        %         %Rate treat Onset     :     11   info: nTrial;
                        %         %Rate treat Offset    :     12   info: nTrial;
                        %         %Face Onset           :     13   info: dist;
                        %         %Face Offset          :     14   info: dist;
                        %         %FaceStim Fixcross    :     15   info: position(1)
                        %         %Tonic Pain reached   :     16   info: nTrial
                        %         %Fixcross Jump        :     17   info:
                        %         %Cooldown end         :     18   info:
                        %         %dummy fixflip        :     22   info: NaN;
                        %         planned trialstart    :     30   info: NaN
                        %         planned trialend      :     31   info: NaN
                        
                        
                        names = {'VAS','Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PlateauReached','FixJump','CoolDown','TrialStart','TrialEnd'};
                        condnum = [-2 -1 0 1:18 30 31];
                        
                        TR = Project.TR;
                        scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
                        scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
                        
                        %% get onsets
                        
                        conds = [4 11 9 18];% (13 = face is gone) 4 = RampDown, 11 = RateTreat, 9=RatePain 18 = CoolDown at very end.
                        if verbalize ==1
                            fprintf('Considering the following conditions:\n')
                            for c = conds
                                fprintf('%s\n',names{condnum==c})
                            end
                        end
                        
                        
                        all_events          = find(ismember(L(:,2),conds));
                        LL = L(all_events,:);
                        
                        if verbalize ==1
                            fprintf('Overall we have %d events.\n',length(all_events))
                        end
                        
                        %Conds
                        
                        %1 -135
                        %2 -90
                        %3 -45
                        %4  CSP
                        %5 45
                        %6 90
                        %7 135
                        %8 180
                        %9 UCS
                        %10 Null
                        %11 RatePain
                        %12 RateRelief
                        %13 BaselineEnd
                        
                        RampdownOnsets   = LL(LL(:,2)==4,1);
                        RateReliefOnsets = LL(LL(:,2)==11,1);
                        RatePainOnsets   = LL(LL(:,2)==9,1);
                        RampdownFinal    = RampdownOnsets(end);
                        RampdownOnsets(end) = [];
                        
                        cc = 0;
                        
                        cc = cc+1;
                        cond_list = L(L(:,2)==13,3);
                        cond = [];
                        cond(cc).name     = 'RampDown';
                        cond(cc).onset    = RampdownOnsets(cond_list < 500);
                        cond(cc).duration = 0;
                        if verbalize ==1
                            fprintf('Number of onsets for %s: %d.\n',cond(1).name,length(cond(1).onset))
                        end
                        cc = cc+1;
                        cond(cc).name     = 'Null';
                        cond(cc).onset    = RampdownOnsets(cond_list ==3000);
                        cond(cc).duration = 0;
                        
                        if run > 1
                            cc = cc+1;
                            cond(cc).name     = 'UCS';
                            cond(cc).onset    = RampdownOnsets(cond_list == 500);
                            cond(cc).duration = 0;
                        end
                        
                        
                        % Rate Pain
                        cc = cc+1;
                        cond(cc).name   = 'RatePain';
                        cond(cc).onset  = RatePainOnsets;
                        cond(cc).duration = L(L(:,2)==10,1)-RatePainOnsets;
                        if verbalize ==1
                            fprintf('Number of onsets for RatePain: %d.\n',length(cond(cc).onset))
                        end
                        if self.id == 6
                            corrdur = [1.5 2 2 1]; %scanner stopped too early, during last rating.
                            cond(cc).duration(end) = corrdur(run);
                        elseif self.id == 32 && run == 1
                            cond(cc).onset(end-1:end) = [];
                            cond(cc).duration(end-1:end) = [];
                        end
                        
                        % Rate Relief
                        cc = cc+1;
                        cond(cc).name     = 'RateRelief';
                        cond(cc).onset    = RateReliefOnsets;
                        cond(cc).duration = ones(1,length(cond(cc).onset)).*5;
                        if self.id == 32 && run ==1
                            cond(cc).duration(end) = 4; % scanner stopped here
                        end
                        if verbalize ==1
                            fprintf('Number of onsets for RateRelief: %d.\n',length(cond(cc).onset))
                        end
                        
                        %                         % CoolDown phase
                        %                         cc = cc+1;
                        %                         RampDown        = L(L(:,2)==4,1);
                        %                         RampDownEnd     = RampDown(end);
                        %                         cond(cc).name   = '999';
                        %                         cond(cc).onset  = RampDownEnd;
                        %                         cond(cc).duration = RatePainOnsets(end)-L(find(L(:,2)==4,1,'last'),1);
                        %                         if self.id == 4
                        %                             cond(cc).duration = 32.6; %logging problem.
                        %                         elseif self.id == 32 && run == 1 %wrong seq, scanner stopped already
                        %                             cond(cc) = [];
                        %                         end
                        %
                        %                         if verbalize ==1
                        %                             fprintf('Number of onsets for CoolDown %d.\n',length(cond(end).onset))
                        %                         end
                        
                        
                        %% housekeeping
                        for counter = 1:length(cond)
                            cond(counter).onset    = cond(counter).onset./TR;
                            cond(counter).duration = cond(counter).duration./TR;
                            cond(counter).tmod      = 0;
                            cond(counter).pmod      = struct('name',{},'param',{},'poly',{});
                        end
                        
                        %% Tuning information as pmod, VonMises and derivative dVM/dkappa
                        conds = -135:45:180;
                        amp   = 1;
                        kappa = 1;
                        delta = .01;
                        VMlookup = zscore(Tuning.VonMises(conds,amp,kappa,0,0));
                        dVMlookup = -zscore((Tuning.VonMises(conds,amp,kappa+delta,0,0)-Tuning.VonMises(conds,amp,kappa-delta,0,0))./(2*delta)); %central difference formula
                        %
                        condlist = cond_list(cond_list < 500);
                        for ntrial = 1:length(cond(1).onset)
                            deltacsp = condlist(ntrial);
                            ind = Project.compute_deltacsp2ind(deltacsp);
                            VM(ntrial)    = VMlookup(ind);
                            dVM(ntrial)   = dVMlookup(ind);
                            VMdVM(ntrial) = VMlookup(ind).*dVMlookup(ind);
                        end
                        
                        % linear regression y ~ 1 + vM + dVM
                        % Interaction is done on 2ndlevel (?)
                        cond(1).pmod  = struct('name',{'vM','dvM','vMxdvM'},'param',{VM dVM VMdVM},'poly',{1 1 1});
                        
                        save(modelpath,'cond');
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 9
                    modelname = 'RampOnsetStick_pmod_rate'; %based on 6, onsets are RampOnsets, but we model RampOnsets as stick (but everything else as box)
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,7)) %model 7 gives onsets as one regressor, On RampDown
                        cond0 = cond;
                        clear cond;
                        %% merge all onsets in one regressor, but keep Null as different (bc no face activity)
                        nreg = {1,[1 3],[1 3],[1 3]}; % only first phase has no UCS at position 3, otherwise we want to merge normal faces and UCS onsets, but keep null as cond(2)
                        [~,ratings] = self.get_rating(run);
                        ratings(1) = []; % take nulltrial out before zscoring
                        ratings = zscore(ratings);
                        onsets = [];
                        for nr = nreg{run}
                            onsets = [onsets(:)' cond0(nr).onset(:)'];
                        end
                        onsets = sort(onsets);
                        cond(1).name = 'RampDown';
                        cond(1).onset = onsets;
                        cond(1).duration = zeros(1,length(onsets));
                        cond(1).tmod = 0;
                        cond(1).pmod = struct('name',{'Relief'},'param',{ratings},'poly',{1});
                        
                        cond(2) = cond0(2); %Nulltrial seperately
                        cond(3) = cond0(end-1); %RatePain
                        cond(4) = cond0(end); %RateRelief
                        save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 10
                    modelname = 'RampOnsetStick_pmod_PR'; %based on 6, onsets are RampOnsets, but we model RampOnsets as stick (but everything else as box) THIS IS FOR FIR, so NO RATINGS
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,7)) %model 7 gives Face onsets as one regressor (called RampDown) and Null, UCS, RateRelief, RatePain as others
                        cond0 = cond;
                        clear cond;
                        %% merge all onsets in one regressor, but keep Null as different (bc no face activity)
                        nreg = {1,[1 3],[1 3],[1 3]}; % only first phase has no UCS at position 3, otherwise we want to merge normal faces and UCS onsets, but keep null as cond(2)
                        
                        PRindex = self.get_pmod_PR(run);
                        PRindex(1) = []; %exclude nulltrial.
                        PRindex = zscore(PRindex);
                        onsets = [];
                        for nr = nreg{run}
                            onsets = [onsets(:)' cond0(nr).onset(:)'];
                        end
                        onsets = sort(onsets);
                        cond(1).name = 'RampDown';
                        cond(1).onset = onsets;
                        cond(1).duration = zeros(1,length(onsets));
                        cond(1).tmod = 0;
                        cond(1).pmod = struct('name',{'PR_ind'},'param',{PRindex},'poly',{1});
                        
                        cond(2) = cond0(2); %Nulltrial seperately
                        cond(3) = cond0(end-1); %RatePain
                        cond(4) = cond0(end); %RateRelief
                        save(modelpath,'cond')
                        %                         %% version for FIR without RatePain or RateRelief
                        %                         %need to merge RampOnset and UCS for that.
                        %                         nreg = 1+sum(ismember(self.condsinphase(run),[500 3000])); %we need to know this to find RateRelief and RatePain behind Onsets, Null, UCS.
                        %                         if run > 1
                        %                             nreg2take = [1 3];
                        %                             onsets = [cond0(1).onset(:)' cond0(3).onset(:)'];
                        %                         else
                        %                             onsets = cond0(1).onset;
                        %                         end
                        %
                        %                         PRindex = self.get_pmod_PR(run);
                        %                         PRindex(1) = []; %exclude nulltrial.
                        %                         onsets = sort(onsets);
                        %                         cond(1).name = 'RampDown';
                        %                         cond(1).onset = onsets;
                        %                         cond(1).duration = zeros(1,length(onsets));
                        %                         cond(1).tmod = 0;
                        %                         cond(1).pmod = struct('name',{'PR_ind'},'param',{PRindex},'poly',{1});
                        %
                        %                         cond(2) = cond0(2);
                        %                         save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                    
                case 11
                    modelname = 'FaceOnset_NoCool'; %based on model_01, onsets are FaceOnsets, we model them as Stick, and kick the Cooldownphase
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,1))
                        set2zero = self.nreliefconds(run);
                        %% set durations to zero
                        for counter = 1:set2zero
                            cond(counter).duration = 0;
                        end
                        if self.id ==32 && run ==1
                            fprintf('Stays like model 1.\n')
                            save(modelpath,'cond')
                        else
                            fprintf('Cutting conds.\n')
                            %get rid of Cooldown
                            CoolCond = self.nreliefconds(run)+3;
                            cond(CoolCond) = [];
                            %get rid of PainRating at the end
                            RatePainCond = self.nreliefconds(run)+1;
                            cond(RatePainCond).onset(end) = [];
                            cond(RatePainCond).duration(end) = [];
                            save(modelpath,'cond')
                        end
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 12
                    modelname = 'FaceOnsetStick_pmod_rate'; %based on 6, onsets are RampOnsets, but we model RampOnsets as stick (but everything else as box)
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,11)) %model 11 has face onsets without Cooldown phase
                        cond0 = cond;
                        clear cond;
                        %% merge all onsets in one regressor, but keep Null as different (bc no face activity)
                        nreg = find(self.condsinphase(run)~=3000); %
                        [~,ratings] = self.get_rating(run);
                        ratings(1) = []; % take nulltrial out before zscoring
                        ratings = zscore(ratings);
                        onsets = [];
                        for nr = nreg(:)'
                            onsets = [onsets(:)' cond0(nr).onset(:)'];
                        end
                        onsets = sort(onsets);
                        cond(1).name = 'RampDown';
                        cond(1).onset = onsets;
                        cond(1).duration = zeros(1,length(onsets));
                        cond(1).tmod = 0;
                        cond(1).pmod = struct('name',{'Relief'},'param',{ratings},'poly',{1});
                        
                        cond(2) = cond0(self.condsinphase(run)==3000); %Nulltrial seperately
                        cond(3) = cond0(end-1); %RatePain
                        cond(4) = cond0(end); %RateRelief
                        save(modelpath,'cond')
                        
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 13
                    modelname = 'FaceOnsetStick_pmod_PR'; %based on 6, onsets are RampOnsets, but we model RampOnsets as stick (but everything else as box)
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,11)) %model 7 gives onsets as one regressor
                        cond0 = cond;
                        clear cond;
                        %% merge normal facetrials and UCS trials in one regressor, but keep Nulltrial seperate
                        nreg = find(self.condsinphase(run)~=3000); %only first phase has no UCS at position 3, otherwise we want to merge normal faces and UCS onsets, but keep null as cond(2)
                        
                        PRindex = self.get_pmod_PR(run);
                        PRindex(1) = []; %exclude nulltrial.
                        PRindex = zscore(PRindex);
                        onsets = [];
                        for nr = nreg(:)'
                            onsets = [onsets(:)' cond0(nr).onset(:)'];
                        end
                        onsets = sort(onsets);
                        cond(1).name = 'RampDown';
                        cond(1).onset = onsets;
                        cond(1).duration = zeros(1,length(onsets));
                        cond(1).tmod = 0;
                        cond(1).pmod = struct('name',{'PR_ind'},'param',{PRindex},'poly',{1});
                        
                        cond(2) = cond0(self.condsinphase(run)==3000); %Nulltrial seperately
                        cond(3) = cond0(end-1); %RatePain
                        cond(4) = cond0(end); %RateRelief
                        save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 15
                    modelname = 'RampDown_6sec_plateau_Rate5'; %based on 6, onsets are RampOnsets
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,6)) %model 7 gives onsets as one regressor
                        cond0 = cond;
                        clear cond;
                        %% merge normal facetrials and UCS trials in one regressor, but keep Nulltrial seperate
                        nreg = find(self.condsinphase(run)~=3000); %only first phase has no UCS at position 3, otherwise we want to merge normal faces and UCS onsets, but keep null as cond(2)
                        pmod = double(self.paradigm{run}.presentation.ucs(2:end));
                        
                        onsets = [];
                        for nr = nreg(:)'
                            onsets = [onsets(:)' cond0(nr).onset(:)'];
                        end
                        onsets = sort(onsets);
                        cond(1).name = 'RampDown';
                        cond(1).onset = onsets;
                        cond(1).duration = 6*ones(1,length(onsets))./self.TR;
                        cond(1).tmod = 0;
                        cond(1).pmod = struct('name',{'PR_ind'},'param',{pmod},'poly',{1});
                        
                        cond(2) = cond0(self.condsinphase(run)==3000); %Nulltrial seperately
                        cond(2).duration = 6./self.TR;
                        cond(3) = cond0(end-1); %RatePain
                        cond(4) = cond0(end); %RateRelief
                        save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 16
                    modelname = 'RampDown_11sec_plateau'; %based on 6, onsets are RampOnsets
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,6)) %model 6 gives onsets as one regressor
                        cond0 = cond;
                        clear cond;
                        %% merge normal facetrials and UCS trials in one regressor, but keep Nulltrial seperate
                        nreg = find(self.condsinphase(run)~=3000); %only first phase has no UCS at position 3, otherwise we want to merge normal faces and UCS onsets, but keep null as cond(2)
                        
                        pmod = double(self.paradigm{run}.presentation.ucs(2:end));
                        
                        onsets = [];
                        for nr = nreg(:)'
                            onsets = [onsets(:)' cond0(nr).onset(:)'];
                        end
                        onsets = sort(onsets);
                        cond(1).name = 'RampDown';
                        cond(1).onset = onsets;
                        cond(1).duration = 11*ones(1,length(onsets))./self.TR;
                        cond(1).tmod = 0;
                        cond(1).pmod = struct('name',{'highlow'},'param',{pmod},'poly',{1});
                        
                        cond(2) = cond0(self.condsinphase(run)==3000); %Nulltrial seperately
                        cond(2).duration = 11./self.TR;
                        cond(3) = cond0(end-1); %RatePain
                        save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 17
                    modelname = 'RampDown_Stick_plus_11sec_plateau'; %based on 6, onsets are RampOnsets
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,6)) %model 7 gives onsets as one regressor
                        cond0 = cond;
                        clear cond;
                        %% merge normal facetrials and UCS trials in one regressor, but keep Nulltrial seperate
                        nreg = find(self.condsinphase(run)~=3000); %only first phase has no UCS at position 3, otherwise we want to merge normal faces and UCS onsets, but keep null as cond(2)
                        
                        pmod =  double(self.paradigm{run}.presentation.ucs(2:end));
                        %self.get_ramptemp(run);
                        %pmod(1) = [];
                        
                        onsets = [];
                        for nr = nreg(:)'
                            onsets = [onsets(:)' cond0(nr).onset(:)'];
                        end
                        onsets = sort(onsets);
                        cond(1).name = 'RampDownStick';
                        cond(1).onset = onsets;
                        cond(1).duration = zeros(1,length(onsets));
                        cond(1).tmod = 0;
                        cond(1).pmod = struct('name',{'highlow'},'param',{pmod},'poly',{1});
                        
                        cond(2).name = 'RampDownPlateau';
                        cond(2).onset = onsets;
                        cond(2).duration = 11*ones(1,length(onsets))./self.TR;
                        cond(2).tmod = 0;
                        cond(2).pmod = struct('name',{},'param',{},'poly',{1});
                        
                        cond(3) = cond0(self.condsinphase(run)==3000); %Nulltrial seperately
                        cond(4) = cond0(end-1); %RatePain
                        save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 18
                    modelname = 'SBS_Schlern'; %RampDown, Boxvcar Plateau, RampUp again..
                    if ~exist(modelpath) || force == 1
                        if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                        if verbalize == 1
                            fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                        end
                        L               = self.get_log(run);
                        %sort things according to time rather than order of being logged
                        [~,i]           = sort(L(:,1),'ascend');
                        L               = L(i,:);
                        % delete all the events that are after the last scanning..
                        scan_times      = L(find(L(:,2) == 0),1);
                        first_scan_time = min(scan_times);
                        last_scan_time  = max(scan_times);
                        L(L(:,1) < first_scan_time,:) = [];
                        L(L(:,1) > last_scan_time,:)  = [];
                        L(:,1)                        = L(:,1) - first_scan_time; % correct time logging so that it's referring to the first scan as 0. (later first col will be from > dummy trials)
                        
                        % correct timings for dummy scans
                        L(1:6,:) = []; %let dummies go
                        first_real_scan = L(1,1);
                        L(:,end+1) = L(:,1); %store it so we still know the original timing
                        L(:,1) = L(:,1)-first_real_scan;
                        
                        %event types are as follows:
                        %         %Pulse Detection      :     0    info: NaN;
                        %         %Ramp down Onset      :     4    info: ror
                        %         %Treatment Plateau    :     5    info: temp
                        %         %Ramp back onset      :     6    info: ror;
                        %         %Rate pain Onset		:     9    info: nTrial;
                        %         %Rate treat Onset     :     11   info: nTrial;
                        %         %Rate treat Offset    :     12   info: nTrial;
                        %         %Tonic Pain reached   :     16   info: nTrial
                        %         %Fixcross Jump        :     17   info:
                        %         %Cooldown end         :     18   info:
                        %         %dummy fixflip        :     22   info: NaN;
                        %         planned trialstart    :     30   info: NaN
                        %         planned trialend      :     31   info: NaN
                        
                        
                        names = {'VAS','Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PlateauReached','FixJump','CoolDown','TrialStart','TrialEnd'};
                        condnum = [-2 -1 0 1:18 30 31];
                        
                        TR = Project.TR;
                        scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
                        scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
                        
                        %% get onsets
                        
                        conds = [4 5 6 9];% (13 = face is gone) 4 = RampDown, 11 = RateTreat, 9=RatePain 18 = CoolDown at very end.
                        if verbalize ==1
                            fprintf('Considering the following conditions:\n')
                            for c = conds
                                fprintf('%s\n',names{condnum==c})
                            end
                        end
                        
                        
                        all_events          = find(ismember(L(:,2),conds));
                        LL = L(all_events,:);
                        
                        if verbalize ==1
                            fprintf('Overall we have %d events.\n',length(all_events))
                        end
                        
                        %Conds
                        
                        %1 -135
                        %2 -90
                        %3 -45
                        %4  CSP
                        %5 45
                        %6 90
                        %7 135
                        %8 180
                        %9 UCS
                        %10 Null
                        %11 RatePain
                        %12 RateRelief
                        %13 BaselineEnd
                        
                        RampdownOnsets   = LL(LL(:,2)==4,1);RampdownOnsets(end) = []; %kick Cooldown already here
                        PlateauReached     = LL(LL(:,2)== 5,1);
                        RampupOnsets     = LL(LL(:,2)== 6,1);
                        RateReliefOnsets = LL(LL(:,2)==11,1);
                        RatePainOnsets   = LL(LL(:,2)==9,1);RatePainOnsets(end) = [];
                        
                        FaceTrials = setdiff(unique(L(L(:,2)==13,3)),3000);
                        %FaceTrials RampDown
                        cond(1) =  struct('name',{'RampDown'},'onset',{RampdownOnsets(2:end)},'duration',{zeros(1,length(RampdownOnsets(2:end)))});
                        %FaceTrials Plateau
                        cond(2) =  struct('name',{'Plateau'},'onset',{RampdownOnsets(2:end)},'duration',{11*ones(1,length(RampdownOnsets(2:end)))});
                        %FaceTrials RampUp
                        cond(3) =  struct('name',{'RampUp'},'onset',{RampupOnsets(2:end)},'duration',{zeros(1,length(RampdownOnsets(2:end)))});
                        
                        %NullTrial RampDown
                        cond(4) =  struct('name',{'RampDown0'},'onset',{RampdownOnsets(1)},'duration',{0});
                        %NullTrial Plateau
                        cond(5) =  struct('name',{'Plateau0'},'onset',{RampdownOnsets(1)},'duration',{11});
                        %NullTrial RampUp
                        cond(6) =  struct('name',{'RampUp0'},'onset',{RampupOnsets(1)},'duration',{0});
                        
                        % Rate Pain
                        RPOffset =  L(L(:,2)==10,1); RPOffset(end) =[];%Cooldown rating, we don't need that
                        RPdur =  RPOffset-RatePainOnsets;
                        cond(7) = struct('name',{'RatePain'},'onset',{RatePainOnsets},'duration',RPdur);
                        
                        if self.id == 6
                            corrdur = [1.5 2 2 1]; %scanner stopped too early, during last rating.
                            cond(end).duration(end) = corrdur(run);
                        elseif self.id == 32 && run == 1
                            cond(end).onset(end-1:end) = [];
                            cond(end).duration(end-1:end) = [];
                        end
                        
                        %                         % Rate Relief
                        %                         cc = cc+1;
                        %                         cond(cc).name     = 'RateRelief';
                        %                         cond(cc).onset    = RateReliefOnsets;
                        %                         cond(cc).duration = ones(1,length(cond(cc).onset)).*5;
                        %                         if self.id == 32 && run ==1
                        %                             cond(cc).duration(end) = 4; % scanner stopped here
                        %                         end
                        %                         if verbalize ==1
                        %                             fprintf('Number of onsets for RateRelief: %d.\n',length(cond(cc).onset))
                        %                         end
                        
                        
                        %% housekeeping
                        for counter = 1:length(cond)
                            cond(counter).onset    = cond(counter).onset./TR;
                            cond(counter).duration = cond(counter).duration./TR;
                            cond(counter).tmod      = 0;
                            cond(counter).pmod      = struct('name',{},'param',{},'poly',{});
                        end
                        
                        % add pmod of ucs or not to RampUp and RampDown.
                        pmodinfo = double(self.paradigm{run}.presentation.ucs(2:end));
                        cond(1).pmod      = struct('name',{'highlow'},'param',{pmodinfo},'poly',{1});
                        cond(3).pmod      = struct('name',{'highlow'},'param',{pmodinfo},'poly',{1});
                        
                        save(modelpath,'cond');
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
                case 19
                    modelname = 'Ramp_Gausspmod_Faces_else'; %based on model_02, Face and Ramp seperately, where Ramp gets 1RegrForAll with pmod Gauss/dGauss
                    
                    if ~exist(modelpath) || force == 1
                        if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                        if verbalize == 1
                            fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                        end
                        L               = self.get_log(run);
                        %sort things according to time rather than order of being logged
                        [~,i]           = sort(L(:,1),'ascend');
                        L               = L(i,:);
                        % delete all the events that are after the last scanning..
                        scan_times      = L(find(L(:,2) == 0),1);
                        first_scan_time = min(scan_times);
                        last_scan_time  = max(scan_times);
                        L(L(:,1) < first_scan_time,:) = [];
                        L(L(:,1) > last_scan_time,:)  = [];
                        L(:,1)                        = L(:,1) - first_scan_time; % correct time logging so that it's referring to the first scan as 0. (later first col will be from > dummy trials)
                        
                        % correct timings for dummy scans
                        L(1:6,:) = []; %let dummies go
                        first_real_scan = L(1,1);
                        L(:,end+1) = L(:,1); %store it so we still know the original timing
                        L(:,1) = L(:,1)-first_real_scan;
                        
                        %event types are as follows:
                        %         %Pulse Detection      :     0    info: NaN;
                        %         %Tracker Onset        :     1
                        %         %Cross (tonic) Onset  :     2    info: position
                        %         %Cross (pain) Onset   :     3    info: position
                        %         %Ramp down Onset      :     4    info: ror
                        %         %Treatment Plateau    :     5    info: temp
                        %         %Ramp back onset      :     6    info: ror;
                        %         %Key Presses          :     7    info: keycode;
                        %         %Tracker Offset       :     8    info: NaN;
                        %         %Rate pain Onset		:     9    info: nTrial;
                        %         %Rate pain Offset     :     10   info: nTrial;
                        %         %Rate treat Onset     :     11   info: nTrial;
                        %         %Rate treat Offset    :     12   info: nTrial;
                        %         %Face Onset           :     13   info: dist;
                        %         %Face Offset          :     14   info: dist;
                        %         %FaceStim Fixcross    :     15   info: position(1)
                        %         %Tonic Pain reached   :     16   info: nTrial
                        %         %Fixcross Jump        :     17   info:
                        %         %Cooldown end         :     18   info:
                        %         %dummy fixflip        :     22   info: NaN;
                        %         planned trialstart    :     30   info: NaN
                        %         planned trialend      :     31   info: NaN
                        
                        
                        names = {'VAS','Text','Pulse','Tracker+','CrossTonic','CrossTreat','RampDown','Plateau','RampBack','Keys','TrackerOff','RatePainOn','RatePainOff','RateReliefOn','RateReliefOff','FaceOn','FaceOff','FaceStimFX','PlateauReached','FixJump','CoolDown','TrialStart','TrialEnd'};
                        condnum = [-2 -1 0 1:18 30 31];
                        
                        TR = Project.TR;
                        scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
                        scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
                        
                        %% get onsets
                        
                        conds = [4 13 11 9 18];% 4 = RampDown, 13 = face cues, 11 = RateTreat, 9=RatePain 18 = CoolDown at very end.
                        if verbalize ==1
                            fprintf('Considering the following conditions:\n')
                            for c = conds
                                fprintf('%s\n',names{condnum==c})
                            end
                        end
                        
                        
                        all_events          = find(ismember(L(:,2),conds));
                        LL = L(all_events,:);
                        
                        if verbalize ==1
                            fprintf('Overall we have %d events.\n',length(all_events))
                        end
                        
                        %Conds
                        
                        %1 -135
                        %2 -90
                        %3 -45
                        %4  CSP
                        %5 45
                        %6 90
                        %7 135
                        %8 180
                        %9 UCS
                        %10 Null
                        %11 RatePain
                        %12 RateRelief
                        %13 BaselineEnd
                        
                        RampdownOnsets   = LL(LL(:,2)==4,1);
                        FaceOnsets    = LL(LL(:,2)==13,1);
                        RateReliefOnsets = LL(LL(:,2)==11,1);
                        RatePainOnsets   = LL(LL(:,2)==9,1);
                        RampdownFinal    = RampdownOnsets(end);
                        RampdownOnsets(end) = [];
                        
                        cc = 0;
                        cond = [];
                        cc = cc+1;
                        cond_list = L(L(:,2)==13,3);
                        
                        cond(cc).name     = 'RampDown';
                        cond(cc).onset    = RampdownOnsets(cond_list < 500);
                        cond(cc).duration = ones(length(cond(cc).onset),1).*6;
                        if verbalize ==1
                            fprintf('Number of onsets for %s: %d.\n',cond(1).name,length(cond(1).onset))
                        end
                        cc = cc+1;
                        cond(cc).name     = 'Null';
                        cond(cc).onset    = RampdownOnsets(cond_list ==3000);
                        cond(cc).duration = 6;
                        
                        if run > 1
                            cc = cc+1;
                            cond(cc).name     = 'UCS';
                            cond(cc).onset    = RampdownOnsets(cond_list == 500);
                            cond(cc).duration = ones(1,length(cond(cc).onset)).*6;
                        end
                        %face cue onsets, single regressors
                        conds_presented = unique(LL(LL(:,2)==13,3));
                       
                        for c = conds_presented(:)';
                            cc = cc+1;
                            trial_ind = find(LL(LL(:,2)==13,3)==c);
                            cond(cc).name     = [mat2str(c) 'Face'];
                            cond(cc).onset    = FaceOnsets(trial_ind);
                            cond(cc).duration = ones(1,length(cond(cc).onset)).*1.5;
                            if verbalize ==1
                                fprintf('Number of onsets for cond %04g: %d.\n',c,length(cond(cc).onset))
                            end
                        end
                        
                        % Rate Pain
                        cc = cc+1;
                        cond(cc).name   = 'RatePain';
                        cond(cc).onset  = RatePainOnsets;
                        cond(cc).duration = L(L(:,2)==10,1)-RatePainOnsets;
                        if verbalize ==1
                            fprintf('Number of onsets for RatePain: %d.\n',length(cond(cc).onset))
                        end
                        if self.id == 6
                            corrdur = [1.5 2 2 1]; %scanner stopped too early, during last rating.
                            cond(cc).duration(end) = corrdur(run);
                        elseif self.id == 32 && run == 1
                            cond(cc).onset(end-1:end) = [];
                            cond(cc).duration(end-1:end) = [];
                        end
                        
                        % Rate Relief
                        cc = cc+1;
                        cond(cc).name     = 'RateRelief';
                        cond(cc).onset    = RateReliefOnsets;
                        cond(cc).duration = ones(1,length(cond(cc).onset)).*5;
                        if self.id == 32 && run ==1
                            cond(cc).duration(end) = 4; % scanner stopped here
                        end
                        if verbalize ==1
                            fprintf('Number of onsets for RateRelief: %d.\n',length(cond(cc).onset))
                        end
                        
                          % CoolDown phase
                        cc = cc+1;
                        RampDown        = L(L(:,2)==4,1);
                        RampDownEnd     = RampDown(end);
                        cond(cc).name   = '999';
                        cond(cc).onset  = RampDownEnd;
                        cond(cc).duration = RatePainOnsets(end)-L(find(L(:,2)==4,1,'last'),1);
                        if self.id == 4
                            cond(cc).duration = 32.6; %logging problem.
                        elseif self.id == 32 && run == 1 %wrong seq, scanner stopped already
                            cond(cc) = [];
                        end
                        
                        %% housekeeping
                        for counter = 1:length(cond)
                            cond(counter).onset    = cond(counter).onset./TR;
                            cond(counter).duration = cond(counter).duration./TR;
                            cond(counter).tmod      = 0;
                            cond(counter).pmod      = struct('name',{},'param',{},'poly',{});
                        end
                        
                        %% Tuning information as pmod, VonMises and derivative dVM/dkappa
                        conds = -135:45:180;
                        
                        amp   = 1;
                        kappa = 1;
                        delta = .01;
                        VMlookup = zscore(Tuning.VonMises(conds,amp,kappa,0,0));
                        dVMlookup = zscore((Tuning.VonMises(conds,amp,kappa+delta,0,0)-Tuning.VonMises(conds,amp,kappa-delta,0,0))./(2*delta)); %central difference formula
                        %
                        condlist = cond_list(cond_list < 500);
                        for ntrial = 1:length(cond(1).onset)
                            deltacsp = condlist(ntrial);
                            ind = Project.compute_deltacsp2ind(deltacsp);
                            VM(ntrial)    = VMlookup(ind);
                            dVM(ntrial)   = dVMlookup(ind);
                        end
                        
                        % linear regression y ~ 1 + vM + dVM
                        % Interaction is done on 2ndlevel
                        cond(1).pmod  = struct('name',{'vM','dvM'},'param',{VM dVM},'poly',{1 1});
                        
                        save(modelpath,'cond');
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
            end
        end
        function [wmcsf, wm, csf] = get_wmcsf_epi(self,nrun)
            force = 1;
            filename = fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_reg_kc%d.mat',self.kickcooldown));
            if ~exist(filename) || force == 1
                
                c2_templ = '^c2meandata.nii';
                c3_templ = '^c3meandata.nii';
                
                segm_dir = strrep(self.path_epi(1,'mean_EPIdartel/'),'data.nii','');
                
                white_matter = spm_select('FPList', segm_dir, c2_templ);
                cs_fluid = spm_select('FPList', segm_dir, c3_templ);
                
                
                WM = spm_vol(white_matter);
                [y1] = spm_read_vols(WM);
                
                CSF = spm_vol(cs_fluid);
                [y2] = spm_read_vols(CSF);
                
                
                rundir = strrep(self.path_epi(nrun),'data.nii','');
                nscans = self.get_total_volumes(nrun);
                if self.kickcooldown
                    nscans = self.get_lastscan_cooldown(nrun);
                end
                files = spm_select('ExtFPList', rundir, '^r.*.nii',1:nscans);
                
                
                V = spm_vol(files);
                y = spm_read_vols(V);
                
                for j = 1:size(files,1)
                    wm(j) = nanmean(nanmean(nanmean(y1.*y(:,:,:,j))));
                    csf(j) = nanmean(nanmean(nanmean(y2.*y(:,:,:,j))));
                end
                
                wmcsf = [wm' csf'];
                savewhere = fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_reg_kc%d.mat',self.kickcooldown));
                save(filename,'wmcsf','wm','csf');
            else
                load(filename);
            end
        end
        function [wmcsf, wm, csf] = get_wmcsf_T1(self,nrun)
            force = 1;
            filename = fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_reg_kc%d.mat',self.kickcooldown));
            if ~exist(filename) || force == 1
                
                c2_templ = '^resize_c2data.nii';
                c3_templ = '^resize_c3data.nii';
                
                segm_dir = strrep(self.path_hr,'data.nii','segm_dartel');
                
                white_matter = spm_select('FPList', segm_dir, c2_templ);
                cs_fluid = spm_select('FPList', segm_dir, c3_templ);
                
                
                WM = spm_vol(white_matter);
                [y1] = spm_read_vols(WM);
                
                CSF = spm_vol(cs_fluid);
                [y2] = spm_read_vols(CSF);
                
                
                rundir = strrep(self.path_epi(nrun),'data.nii','');
                nscans = self.get_total_volumes(nrun);
                if self.kickcooldown
                    nscans = self.get_lastscan_cooldown(nrun);
                end
                files = spm_select('ExtFPList', rundir, '^rdata.*.nii',1:nscans);
                
                
                V = spm_vol(files);
                y = spm_read_vols(V);
                
                %this is different from AT's script on segmented epis, bc. diff resolution in T1 and
                %epis
                for j = 1:size(files,1)
                    %                     wmdummypath  = fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_wm_dummy%d.nii',j));
                    %                     csfdummypath =  fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_csf_dummy%d.nii',j));
                    %                     spm_imcalc([y1,y(:,:,:,j)],wmdummypath,'i1*i2');
                    %
                    %                     spm_imcalc([y2,y(:,:,:,j)],csfdummypath,'i1*i2');
                    wm(j) = nanmean(nanmean(nanmean(y1.*y(:,:,:,j))));
                    csf(j) = nanmean(nanmean(nanmean(y2.*y(:,:,:,j))));
                end
                
                wmcsf = [wm' csf'];
                savewhere = fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_reg_T1_kc%d.mat',self.kickcooldown));
                save(filename,'wmcsf','wm','csf');
            else
                load(filename);
            end
        end
        function resize_wmcsf(self)
            matlabbatch = [];
            matlabbatch{1}.spm.spatial.coreg.write.ref = {self.path_meanepi};%{'/projects/crunchie/treatgen/data/sub004/run001/mrt/meandata.nii,1'};
            matlabbatch{1}.spm.spatial.coreg.write.source = {
                strrep(self.path_data(0,'mrt'),'data.mat','/segm_dartel/c2data.nii')
                strrep(self.path_data(0,'mrt'),'data.mat','/segm_dartel/c3data.nii')
                %                  '/projects/crunchie/treatgen/data/sub004/run000/mrt/segm_dartel/c2data.nii,1'
                %                  '/projects/crunchie/treatgen/data/sub004/run000/mrt/segm_dartel/c3data.nii,1'
                };
            matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'resize_';
            spm_jobman('run', matlabbatch);
        end
        
        function [kickcooldownopts, wmcsfopts,varargout] = get_model_specs(self,modelnum)
            kickcooldownopts = self.kickcooldown;
            wmcsfopts = self.wmcsfregressors;
            tsda_corr = 0;
            switch modelnum
                case 1
                    kickcooldownopts = 1;
                    wmcsfopts = 0;
                case 5
                    kickcooldownopts = 1;
                    wmcsfopts = 0;
                case 6
                    kickcooldownopts = 1;
                    wmcsfopts = 0;
                case 7
                    kickcooldownopts = 1;
                    wmcsfopts = 0;
                case 8
                    kickcooldownopts = 1;
                    wmcsfopts = 0;
                case 9
                    kickcooldownopts = 1;
                    wmcsfopts = 1;
                otherwise
                    kickcooldownopts = 1;
                    wmcsfopts = 0;
                    tsda_corr = 0;
            end
            varargout{1} = tsda_corr;
        end
        function [num_onsetfile] = get_num_onsetfile(self,modelnum)
            switch modelnum
                case 20
                    num_onsetfile = 2;
                 case 21
                    num_onsetfile = 19;   
                otherwise
                    num_onsetfile =3;
            end
            
        end
        function FitFIR(self,nrun,model_num)
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            empty1stlevel    =   1;
            FIRparam         =   Project.orderfir;
            pmod_Rate        =   0;
            pmod_VM          =   0;
            onset_modelnum   =   3;
            All1Regr         =   0;
            wmcsfr           =   0;
            tsda_corr        =   1;
            
            
            if model_num == 3
                onset_modelnum = 10;
                pmod_Rate = 1;
                All1Regr  = 1;
            end
            if model_num == 5
                pmod_VM = 1;
                onset_modelnum = 7;
                All1Regr  = 1;
            end
            if model_num == 6
                pmod_Rate = 1;
                onset_modelnum = 9;
                All1Regr = 1;
            end
            
            spm_dir  = strrep(self.dir_spmmat(nrun(1),model_num),'chrf_00',sprintf('FIR_%02d_10conds_00',FIRparam));
            path_spmmat = fullfile(spm_dir,'SPM.mat');
            
            if ~exist(path_spmmat);mkdir(spm_dir);end
            
            if empty1stlevel == 1
                if exist(spm_dir); system(sprintf('rm -fr %s*',strrep(spm_dir,'//','/')));end %this is AG style.
            end
            clear matlabbatch
            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            %             matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            %             matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            se = 0;
            for session = nrun(:)'
                se = se+1;
                %load files using ...,1, ....,2 format
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans  = cellstr(spm_select('expand',self.path_epi(session,'r')));
                lastscan2include = size(matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans,1);
                if self.kickcooldown
                    lastscan2include = self.get_lastscan_cooldown(session);
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans =  matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans(1:lastscan2include,:);
                
                %load the onsets
                dummy                                              = load(self.path_model(session,onset_modelnum));
                if self.kickcooldown == 0
                    switch session
                        case 1
                            dummy.cond = dummy.cond([1:9 end]); % 8 faces + t0
                        case 2
                            dummy.cond = dummy.cond([1:3 end]);  % CSP, UCS
                        case 3
                            dummy.cond = dummy.cond([1:10 end]); % 8 faces  + UCS + t0
                        case 4
                            dummy.cond = dummy.cond([1:10 end]); % 8 faces + UCS + t0
                    end
                elseif self.kickcooldown == 1
                    if pmod_VM == 1
                        if session == 1
                            dummy.cond = dummy.cond(1:2); %single face regressor with pmod inside, plus nulltrial
                        elseif ismember(session,[3 4])
                            dummy.cond = dummy.cond(1:3);  %single face regressor with pmod inside, plus nulltrial and UCS
                        else
                            fprintf('You entered a session that does probably not make sense. please check.\n')
                        end
                    elseif pmod_Rate == 1
                        %it's all taken care of in PrepareCondFile method 10;
                        
                    else
                        dummy.cond = dummy.cond(1:self.nreliefconds(session)); % 8 faces + UCS + t0 (if there)
                    end
                end
                
                for c = 1:numel(dummy.cond)
                    dummy.cond(c).duration  = 0;
                    dummy.cond(c).onset     = dummy.cond(c).onset -2; %take 2 TRs before first stimulus appears. (onsets are already in unit = TR, so -2 is correct, not -2.*TR).
                    %                     if pmod_Rate == 1
                    %                         if str2double(dummy.cond(c).name)< 900 %don't do it for nulltrial and finalramp, we only have 1 rating here.
                    %                             dummy.cond(c).pmod = self.get_pmod_rating(session,str2double(dummy.cond(c).name));
                    %                         end
                    %                     end
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).cond   = dummy.cond;
                
                %load nuissance parameters
                nuis                                               = self.get_param_motion(session);
                nuis                                               = nuis(1:lastscan2include,:);
                nuis                                               = zscore(nuis);%zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
                for nNuis = 1:size(nuis,2)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis).val   = nuis(:,nNuis);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis).name  = mat2str(nNuis);
                end
                if wmcsfr == 1
                    load(fullfile(self.pathfinder(self.id,session),'midlevel',sprintf('wmcsf_reg_kc%d.mat',self.kickcooldown)),'wm','csf');
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+1).val   = wm';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+1).name  = 'wm';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+2).val   = csf';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+2).name  = 'csf';
                end
                if tsda_corr == 1
                    nui=load(fullfile(self.pathfinder(self.id,session),'midlevel',sprintf('tsdiff_nui.mat')));
                    nui = nui.nui;
                    nnuis_now = size(matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress,2);
                    
                    for tsda_col = 1:size(nui,2)
                        nnuis_now= nnuis_now+1;
                        matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nnuis_now).val   = nui(:,tsda_col);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nnuis_now).name  = ['tsda_' mat2str(tsda_col)];
                    end
                end
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).multi               = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).multi_reg           = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).hpf                 = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                              = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length                  = Project.TR*FIRparam;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order                   = FIRparam;
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                           = -Inf; %.8
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = cellstr(fullfile(self.path_data(0),'mrt','s3_ss_data.nii'));%{''};%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'AR(1)';
            %estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat            = {path_spmmat};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
        end
        %         function CreateModels(self,runs)
        %             %%%%%%%%%%%%%%%%%%%%%%
        %             model_num  = 1;
        %             for run = runs
        %                 model_path = self.path_model(run,model_num);
        %                 if ~exist(fileparts(model_path));mkdir(fileparts(model_path));end
        %                 [scan,id]  = self.StimTime2ScanUnit(run);
        %                 counter    = 0;
        %                 for current_condition = unique(id)
        %                     counter                = counter + 1;
        %                     cond(counter).name     = mat2str(current_condition);
        %                     cond(counter).onset    = scan(id == current_condition);
        %                     cond(counter).duration = zeros(1,length(cond(counter).onset));
        %                     cond(counter).tmod     = 0;
        %                     cond(counter).pmod     = struct('name',{},'param',{},'poly',{});
        %                 end
        %                 save(model_path,'cond');
        %                 %%%%%%%%%%%%%%%%%%%%%%
        %             end
        %         end
        function FitHRF(self,nrun,model_num)
            force_delete = 1;
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            [kickcd,wmcsfr,tsda_corr] = self.get_model_specs(model_num);
            model_num_onsets = self.get_num_onsetfile(model_num);
            
            %set spm dir:
            spm_dir  = self.dir_spmmat(nrun(1),model_num);
            path_spmmat = self.path_spmmat(nrun(1),model_num);
            if ~exist(path_spmmat)
                mkdir(spm_dir)
            else
                if force_delete == 1
                    system(sprintf('rm -fr %s*',strrep(spm_dir,'//','/')));
                end %this is AG style.
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            %             matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            %             matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            se = 0;
            for session = nrun
                se = se+1;
                %load files using ...,1, ....,2 format
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans  = cellstr(spm_select('expand',self.path_epi(session,'r')));%use always the realigned data.
                lastscan2include = size(matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans,1);
                if kickcd == 1
                    lastscan2include = self.get_lastscan_cooldown(session);
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans =  matlabbatch{1}.spm.stats.fmri_spec.sess(se).scans(1:lastscan2include,:);
                %load the onsets
                dummy                                              = load(self.path_model(session,model_num_onsets));
                %kick cooldown onset from cond file
                if strcmp(dummy.cond(end).name,'999')
                    dummy.cond(end)=[];
                end
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).cond   = dummy.cond;
                %load nuissance parameters
                nuis                                               = self.get_param_motion(session);
                nuis                                               = nuis(1:lastscan2include,:);
                nuis                                               = zscore(nuis);%zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
                for nNuis = 1:size(nuis,2)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis).val   = nuis(:,nNuis);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis).name  = mat2str(nNuis);
                end
                if tsda_corr == 1
                    nui=load(fullfile(self.pathfinder(self.id,session),'midlevel',sprintf('tsdiff_nui.mat')));
                    nui = nui.nui;
                    nnuis_now = size(matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress,2);
                    
                    for tsda_col = 1:size(nui,2)
                        nnuis_now= nnuis_now+1;
                        matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nnuis_now).val   = nui(:,tsda_col);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nnuis_now).name  = ['tsda_' mat2str(tsda_col)];
                    end
                end
                if wmcsfr == 1
                    load(fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_reg_kc%d.mat',kickcd)),'wm','csf');
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+1).val   = wm';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+1).name  = 'wm';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+2).val   = csf';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+2).name  = 'csf';
                end
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).multi               = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).multi_reg           = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(se).hpf                 = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                              = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs                  = self.derivatives;%we have [0 0], [ 1 0] or [ 1 1] for 1, 2, or 3 regressors.
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                           = -Inf; %default .8
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = cellstr(spm_select('expand',strrep(self.path_skullstrip,'ss_data.nii','s3_ss_data.nii')));%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'none'; %AR(1)';
            %estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat            = {path_spmmat};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
            %
            %normalize the beta images right away
            %             beta_images = self.path_beta(nrun(1),model_num,'');%'' => with no prefix
            %             self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)
            %             beta_images = self.path_beta(nrun(1),model_num,'wCAT_');%smooth the normalized images.
            %             self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_ww_')
            %             %delete wCAT files, keep s6_wCAT
            %             fprintf('Deleting unsmoothed normalized files..\n');
            %             for n = 1:size(beta_images,1)
            %                 system(sprintf('rm %s',beta_images(n,:)));
            %             end
        end
        function [total_cons] = CreateContrasts(self,nrun,model_num)
            
            path_spm = self.path_spmmat(nrun(1),model_num);
            %             path_con = strrep(path_spm,'SPM.mat',sprintf('con_%04d.nii',con_num));
            %
            matlabbatch = [];
            n = 0;
            name = [];
            vec = [];
            convec = struct([]);
            
            if model_num  == 2 %needs special treatment because both face and relief onset are modelled..
                %
            elseif ismember(model_num,[7 8 9 10 12 13 15 16]) %% models with pmods, so we will skip the main effect / 1st regressor.
                load([self.dir_spmmat(nrun(1),model_num) 'SPM.mat'])
                for pmod = 1:size(SPM.Sess.U(1).P,2) %loop through the pmods we defined
                    
                    vec = zeros(1,self.get_Nbetas(nrun,model_num)-1); %-1 so that phase constant is left out for now.
                    if nrun == 3 && self.id ~= 15
                        vec = zeros(1,self.get_Nbetas(nrun,model_num)./2-1);
                    end
                    n = n+ 1;
                    name{n} = SPM.Sess.U(1).P(pmod).name;
                    vec(1+pmod) = 1; %skip the main effect
                    if nrun == 3  && self.id ~= 15
                        vec = repmat(vec,1,2);
                    end
                    convec{n} = vec;
                end
            elseif ismember(model_num,17)
                regr2take = 1:3;
                regrnames = {'MEdown','PMODdown','Plateau'};
                for regr = regr2take(:)' %loop through the regressors we might want to look at @ 2ndlevel
                    vec = zeros(1,self.get_Nbetas(nrun,model_num)-1); %-1 so that phase constant is left out for now.
                    if nrun == 3 && self.id ~= 15
                        vec = zeros(1,self.get_Nbetas(nrun,model_num)./2-1);
                    end
                    n = n+ 1;
                    name{n} = regrnames{n};
                    vec(regr) = 1; %put 1 to nth regressor
                    if nrun == 3  && self.id ~= 15
                        vec = repmat(vec,1,2);
                    end
                    convec{n} = vec;
                end
            elseif ismember(model_num,18)
                regr2take = 1:5;
                regrnames = {'MEdown','PMODdown','Plateau','MEup','PMODup'};
                for regr = regr2take(:)' %loop through the regressors we might want to look at @ 2ndlevel
                    vec = zeros(1,self.get_Nbetas(nrun,model_num)-1); %-1 so that phase constant is left out for now.
                    if nrun == 3 && self.id ~= 15
                        vec = zeros(1,self.get_Nbetas(nrun,model_num)./2-1);
                    end
                    n = n+ 1;
                    name{n} = regrnames{n};
                    vec(regr) = 1; %put 1 to nth regressor
                    if nrun == 3  && self.id ~= 15
                        vec = repmat(vec,1,2);
                    end
                    convec{n} = vec;
                end
            else
                %% contrast 1: CSdiff
                % CSP > CSN (or UCS > CSN in Cond)
                n = n + 1;
                vec = zeros(1,self.get_Nbetas(nrun,model_num)-1);
                name{n} = 'CSdiff';
                switch nrun
                    case 1
                        vec([4 8])= [1 -1];
                    case 2
                        vec([2 1])= [1 -1];
                    case 3
                        if self.id ~=15
                            vec = zeros(1,self.get_Nbetas(nrun,model_num)./2-1); %replace upper vec, because already double
                            vec([4 8])= [1 -1];
                            vec = repmat(vec,1,2);
                        else %self.id==15
                            vec([4 8])= [1 -1];
                        end
                end
                convec{n} = vec;
                %% Contrasts 2:9 or 2:3: single conditions, esp. used to collect both sessions from run 3/4
                
                thisphaseconds = intersect(self.realconds,unique(self.get_paradigm(nrun).presentation.dist));
                vec0 = zeros(1,self.get_Nbetas(nrun,model_num)-1);
                
                switch nrun
                    case 1
                        for cond = thisphaseconds(:)'
                            vec = vec0; %reset vector, otherwise we accumulate cons [1 0 0 ..] [1 1 0 ..].
                            n = n + 1;
                            name{n} = sprintf('%04d',cond);
                            beta_ind = self.get_beta_index(nrun,model_num,cond);
                            vec(beta_ind) = 1;
                            convec{n} = vec; %for phase 1 and 2 only one session sessions.
                        end
                    case 2
                        % correct this, bc 500 doesnt get spit out above.
                        thisphaseconds = [500 180];
                        for cond = thisphaseconds(:)'
                            vec = vec0; %reset vector, otherwise we accumulate cons [1 0 0 ..] [1 1 0 ..].
                            n = n + 1;
                            name{n} = sprintf('%04d',cond);
                            beta_ind = self.get_beta_index(nrun,model_num,cond);
                            vec(beta_ind) = 1;
                            convec{n} = vec;
                        end
                    case 3
                        if self.id ~=15
                            vec0 = zeros(1,self.get_Nbetas(nrun,model_num)./2-1); %otherwise already double.
                        end
                        for cond = thisphaseconds(:)'
                            vec = vec0; %reset vector, otherwise we accumulate cons [1 0 0 ..] [1 1 0 ..].
                            n = n + 1;
                            name{n} = sprintf('%04d',cond);
                            beta_ind = self.get_beta_index(nrun,model_num,cond);
                            vec(beta_ind) = 1;
                            if self.id ~= 15
                                convec{n} = repmat(vec,1,2);
                            elseif self.id == 15
                                convec{n} = vec;
                            end
                        end
                end
                %% last contrasts: rating. if ever want to look at it...
                n = n+1;
                vec = vec0;
                vec(self.nreliefconds(nrun)+2) = 1; %thats where the RateRelief regressor is.
                name{n} = 'RateRelief';
                if nrun == 3 && self.id ~=15
                    convec{n} = repmat(vec,1,2);
                else
                    convec{n} = vec;
                end
            end
            %% take care that everything adds to 1 / 0 BEFORE creating neg. versions.
            
            for co = 1:n
                if  sum(convec{co}) == 0 % differential, e.g. CSP>CSN [1 -1], so we need to divide those seperately, otherwise we devide by sum=0;
                    convec{co} =convec{co}./sum(convec{co}>0);
                else
                    convec{co} = convec{co}./sum(convec{co});
                end
                %pad with number of session constants
                if self.id == 15 && nrun == 3
                    convec{co} = [convec{co} zeros(1,self.nsessions(nrun)-1)];
                else
                    convec{co} = [convec{co} zeros(1,self.nsessions(nrun))];
                end
            end
            %% add same contrasts' negative version
            for co = 1:n
                convec{co+n} = -1.*convec{co};
                name{co+n}   = [name{co} '_neg'];
            end
            %%
            total_cons = numel(convec);
            % build the batch
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(path_spm);
            for co = 1:numel(name)
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.name    = name{co};
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec  = convec{co};
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
                %sanity check
                if length(convec{co})~= self.get_Nbetas(nrun,model_num)
                    fprintf('Problem with length of convec at contrast %d. Please debug.\n',co);
                    keyboard
                end
            end
            matlabbatch{1}.spm.stats.con.delete = 1;
            A=[];for co = 1:numel(matlabbatch{1}.spm.stats.con.consess);A = [A;matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec];end;
            %             A
            spm_jobman('run',matlabbatch);
        end
        function [total_cons] = CreateContrasts_FIR(self,nrun,model_num,varargin)
            if nargin > 3
                startbin = varargin{1};
                binwin = varargin{2};
            else
                startbin = 4;
                binwin = 4;
            end
            order = self.orderfir;
            nsessions = self.nsessions(nrun); %used to repmat convec later
            if self.id == 15
                nsessions = 1;
            end
            Nnuis = 6;
            if ismember(model_num,[44])
                Nnuis = 8; % with wm and csf regressors
            end
            
            path_spm = [self.path_FIR(nrun,model_num,order,'10conds') 'SPM.mat'];
            load(path_spm);
            matlabbatch = [];
            n = 0;
            name = [];
            convec = struct([]);
            
            Nbetas  = size(SPM.xX.X,2);
            vec0 = zeros(1,(Nbetas./nsessions)-1); %-1 so that phase constant is left out for now.
            
            %
            if ismember(model_num,[5]) %models with PMOD VonMisesonsetsofint(onsetsofint(:,2)==lookup(1,n),1)
                condnames = {'FaceMain','VM','dVM'};
                conds2take = {2:3,[],2:3};
            elseif model_num == 3
                condnames = {'FaceMain','Rating'};%Ratings as pmod included
                conds2take = {2 2 2};
            elseif model_num == 6 %pmod PRind, want to look at both ME and PMOD, so taking 1:2
                condnames = {'FaceMain','PRind'};%Ratings as pmod included
                conds2take = {1:2,1:2,1:2}; %
            else
                if nrun == 2
                    condnames = {'180','500'};
                else
                    condnames = {'-135','-90','-45','0','45','90','135','180','500'};
                end
                conds2take = {1:8,[2 1],1:8};
            end
            
            %% normal single bins collection to prepare full 14bin 2nd level
            % loop through cons and get betas together.
            c = 0;
            for cond_ind = conds2take{nrun}
                c = c+1;
                for bin = 1:order
                    n = n + 1;
                    getcond = (cond_ind-1)*order + bin; %all cons+bins before this con's bin need to be skipped
                    vec = vec0;
                    vec(getcond) = 1;
                    vec = repmat(vec,1,nsessions);
                    vec = padarray(vec,[0 nsessions],0,'post'); %his is done to compare it to Nbetas later. Not necessary per se
                    convec{n} = vec;
                    name{n} = sprintf('bin%02d.%s',bin,condnames{cond_ind});
                end
            end
            % make convecs sum to 1 for positive t-contrasts
            for co = 1:n
                convec{co} = convec{co}./sum(convec{co});
            end
            n1 = n; %first contrast
            
            %% add CSdiff if applicable
            if ismember(model_num,[1 4 44])
                % CSP vs CSN
                cs_ind = [4 8; 2 1; 4 8];
                for bin = 1:order
                    n = n + 1;
                    csp_bin = (cs_ind(nrun,1)-1)*order + bin; %Nth bin beta of CSP
                    csn_bin = (cs_ind(nrun,2)-1)*order + bin; %Nth bin beta of CSN
                    vec = vec0;
                    vec(csp_bin) = 1;
                    vec(csn_bin) = -1;
                    if nrun ==3
                        if self.id ~= 15
                            vec = repmat(vec,1,2);
                        end
                    end
                    if self.id == 15
                        nsessions = 1;
                    end
                    vec = padarray(vec,[0 nsessions],0,'post');
                    convec{n} = vec;
                    name{n} = sprintf('bin%02d.CSdiff',bin);
                end
                % make convecs sum to 1 and -1 for t-contrasts
                for co = (n1+1):n
                    convec{co} = convec{co}./nsessions;
                end
            end
            % get UCS cons
            if nrun == 3
                cond_ind = 9;
                for bin = 1:order
                    n = n + 1;
                    getcond = (cond_ind-1)*order + bin; %all cons+bins before this con's bin need to be skipped
                    vec = vec0;
                    vec(getcond) = 1;
                    vec = repmat(vec,1,nsessions);
                    vec = padarray(vec,[0 nsessions],0,'post'); %his is done to compare it to Nbetas later. Not necessary per se
                    convec{n} = vec./nsessions;
                    name{n} = sprintf('bin%02d.%s',bin,condnames{cond_ind});
                end
            end
            %% bin win thing
            conds2take ={1:8,[2 1],1:9};
            n2 = n;
            if ismember(model_num,[4 44])
                %loop through cons and get betas together.
                %SINGLE CONS ge-bin-ed
                for cond_ind = conds2take{nrun}
                    n = n + 1;
                    binfo1 =  self.findcon_FIR(order,cond_ind,startbin);%all cons+bins before this con's bin need to be skipped
                    binfo2 = self.findcon_FIR(order,cond_ind,startbin+binwin-1);
                    vec = vec0;
                    vec(binfo1:binfo2) = 1;
                    vec = repmat(vec,1,nsessions);
                    vec = padarray(vec,[0 nsessions],0,'post'); %his is done to compare it to Nbetas later. Not necessary per se
                    convec{n} = vec;
                    name{n} = sprintf('bin%02d.win%02d.%s',startbin,binwin,condnames{cond_ind});
                end
                if nrun ==3
                    if self.id ~= 15
                        vec = repmat(vec,1,2);
                    end
                end
                if self.id == 15
                    nsessions = 1;
                end
                %
                cs_ind = [4 8; 2 1; 4 8];
                n = n+1;
                vec = vec0;
                csp_bins =  self.findcon_FIR(order,cs_ind(nrun,1),startbin):self.findcon_FIR(order,cs_ind(nrun,1),startbin+binwin-1); %Nth bin beta of CSP
                csn_bins =  self.findcon_FIR(order,cs_ind(nrun,2),startbin):self.findcon_FIR(order,cs_ind(nrun,2),startbin+binwin-1); %Nth bin beta of CSN
                vec(csp_bins) =  1;
                vec(csn_bins) = -1;
                if nrun ==3
                    if self.id ~= 15
                        vec = repmat(vec,1,2);
                    end
                end
                if self.id == 15
                    nsessions = 1;
                end
                vec = padarray(vec,[0 nsessions],0,'post');
                convec{n} = vec;
                name{n} = sprintf('bin%02d.win%02d.CSdiff',startbin,binwin);
                
                %make convecs sum to 1 and -1 for t-contrasts
                for co = (n2+1):n
                    convec{co} = convec{co}./binwin/nsessions;
                end
            end
            %%
            %             % giant Guassian with certain bins included
            %             if nrun ~= 2
            %                 n = n+1;
            %                 gauss_lookup = spm_Npdf(1:8,4)-mean(spm_Npdf(1:8,4));
            %                 vec = vec0;
            %                 for condcount = 1:8
            %                     take_bins =  self.findcon_FIR(order,condcount,startbin):self.findcon_FIR(order,condcount,startbin+binwin-1); %Nth bin beta of 1st, just as template.
            %                     vec(take_bins) = ones(1,numel(take_bins))*gauss_lookup(condcount);
            %                 end
            %                 if nrun ==3
            %                     if self.id ~= 15
            %                         vec = repmat(vec,1,2);
            %                     end
            %                 end
            %                 if self.id == 15
            %                     nsessions = 1;
            %                 end
            %                 vec = padarray(vec,[0 nsessions],0,'post');
            %                 convec{n} = vec;
            %                 name{n} = sprintf('bin%02d.win%02d.Gauss',startbin,binwin);
            %                 convec{n} = convec{n}./sum(convec{n}(convec{n}>0));
            %             end
            %%
            %visualization and sanity check
            figure(100);
            if nrun == 1
                clf
            end
            subplot(1,3,nrun)
            for n = 1:numel(convec)
                plot(convec{n}+n)
                hold on
            end
            for c = 1:max(conds2take{nrun})
                l=line(repmat((order*(c-1)+1),1,2),ylim);set(l,'Color','k','LineStyle',':')%line on first bin of each condition
            end
            l=line(repmat(max(conds2take{nrun})*order+1,1,2),ylim);set(l,'Color','k','LineWidth',1);%line where contrasts stop being relevant
            l=line(repmat(Nbetas./nsessions-1,1,2),ylim);set(l,'Color','k','LineWidth',2); %end of vec for 1 session
            
            set(gca,'XTick',1:order:order*max(conds2take{nrun}),'XTickLabel',condnames,'YTick',1:order:numel(convec)-1,'YTickLabel',name(1:order:numel(convec)-1))
            title(sprintf('Run %d',nrun));
            %%
            % build the batch
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(path_spm);
            for co = 1:numel(name)
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.name    = name{co};
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec  = convec{co};
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
                %sanity check
                if length(convec{co})~= Nbetas
                    fprintf('Problem with length of convec at contrast %d. Please debug.\n',co);
                    keyboard
                end
            end
            total_cons = numel(name);
            
            matlabbatch{1}.spm.stats.con.delete = 1;
            
            spm_jobman('run',matlabbatch);
            %             for nrun = runs(:)'
            %                 fprintf('\n\nTook %05.2f minutes for 1stLevel contrasts all conds FIR model 0, sub %02d to %02d, run %02d.\n',tocc(nrun)./60,subs(1),subs(end),nrun);
            %             end
            %             fprintf('\n\nTook %05.2f minutes (%05.2f hours) for 1stLevel contrasts all conds FIR model 0, sub %02d to %02d, run %02d to %02d.\n',done./60,done./3600,subs(1),subs(end),runs(1),runs(end));
            
        end
        function [total_cons] = CreateContrasts_FIR_tsda(self,nrun,model_num,varargin)
            if nargin > 3
                startbin = varargin{1};
                binwin = varargin{2};
            else
                startbin = 4;
                binwin = 4;
            end
            order = self.orderfir;
            nsessions = self.nsessions(nrun); %used to repmat convec later
            if self.id == 15
                nsessions = 1;
            end
            
            
            
            
            path_spm = [self.path_FIR(nrun,model_num,order,'10conds') 'SPM.mat'];
            load(path_spm);
            matlabbatch = [];
            n = 0;
            name = [];
            convec = struct([]);
            
            Nbetas  = size(SPM.xX.X,2);
            Nmov= 6;
            Ntsda = load(fullfile(self.path_project,'midlevel','nuisance_N35.mat'));
            Ntsda = Ntsda.all(Ntsda.subs==self.id,:);
            switch nrun
                case 1
                    nconds = 9;
                    vec0 = zeros(1,nconds*order);
                case 2
                    nconds = 3;
                    vec0 = zeros(1,nconds*order);
                case 3
                    nconds = 10;
                    vec0 = zeros(1,nconds*order);
            end
            %
            
            if nrun == 2
                condnames = {'180','500'};
            else
                condnames = {'-135','-90','-45','0','45','90','135','180','500'};
            end
            conds2take = {1:8,[2 1],1:8};
            
            if nrun == 3
                realrun = 3:4;
            else
                realrun = nrun;
            end
            %% normal single bins collection to prepare full 14bin 2nd level
            % loop through cons and get betas together.
            c = 0;
            for cond_ind = conds2take{nrun}
                c = c+1;
                for bin = 1:order
                    n = n + 1;
                    getcond = (cond_ind-1)*order + bin; %all cons+bins before this con's bin need to be skipped
                    vec_pure = vec0;
                    vec_pure(getcond) = 1;
                    %now add all nuis
                    vec = [];
                    for nr = realrun(:)'
                        vec =  [vec vec_pure zeros(1,(Nmov+Ntsda(nr)))];
                    end
                    %add session constants
                    if self.id == 15
                        nsessions = 1;
                    end
                    vec = padarray(vec,[0 nsessions],0,'post'); %his is done to compare it to Nbetas later. Not necessary per se
                    
                    %check Nbeta
                    if ~isequal(Nbetas,length(vec))
                        keyboard;
                    end
                    convec{n} = vec;
                    name{n} = sprintf('bin%02d.%s',bin,condnames{cond_ind});
                end
            end
            % make convecs sum to 1 for positive t-contrasts
            for co = 1:n
                convec{co} = convec{co}./sum(convec{co});
            end
            n1 = n; %first contrast
            
            %% add CSdiff if applicable
            if ismember(model_num,[1 4 44 40])
                % CSP vs CSN
                cs_ind = [4 8; 2 1; 4 8];
                for bin = 1:order
                    n = n + 1;
                    csp_bin = (cs_ind(nrun,1)-1)*order + bin; %Nth bin beta of CSP
                    csn_bin = (cs_ind(nrun,2)-1)*order + bin; %Nth bin beta of CSN
                    vec_pure = vec0;
                    vec_pure(csp_bin) = 1;
                    vec_pure(csn_bin) = -1;
                    
                    %now add all nuis
                    vec = [];
                    for nr = realrun(:)'
                        vec =  [vec vec_pure zeros(1,(Nmov+Ntsda(nr)))];
                    end
                    %add session constants
                    if self.id == 15
                        nsessions = 1;
                    end
                    vec = padarray(vec,[0 nsessions],0,'post'); %his is done to compare it to Nbetas later. Not necessary per se
                    %check Nbeta
                    if ~isequal(Nbetas,length(vec))
                        keyboard;
                    end
                    convec{n} = vec;
                    name{n} = sprintf('bin%02d.CSdiff',bin);
                end
                % make convecs sum to 1 and -1 for t-contrasts
                for co = (n1+1):n
                    convec{co} = convec{co}./nsessions;
                end
            end
            % get UCS cons
            if nrun == 3
                cond_ind = 9;
                for bin = 1:order
                    n = n + 1;
                    getcond = (cond_ind-1)*order + bin; %all cons+bins before this con's bin need to be skipped
                    vec_pure = vec0;
                    vec_pure(getcond) = 1;
                    %now add all nuis
                    vec = [];
                    for nr = realrun(:)'
                        vec =  [vec vec_pure zeros(1,(Nmov+Ntsda(nr)))];
                    end
                    %add session constants
                    if self.id == 15
                        nsessions = 1;
                    end
                    vec = padarray(vec,[0 nsessions],0,'post'); %his is done to compare it to Nbetas later. Not necessary per se
                    %check Nbeta
                    if ~isequal(Nbetas,length(vec))
                        keyboard;
                    end
                    convec{n} = vec./nsessions;
                    name{n} = sprintf('bin%02d.%s',bin,condnames{cond_ind});
                end
            end
            %% bin win thing
            conds2take ={1:8,[2 1],1:9};
            n2 = n;
            if ismember(model_num,[4 44 40])
                %loop through cons and get betas together.
                %SINGLE CONS ge-bin-ed
                for cond_ind = conds2take{nrun}
                    n = n + 1;
                    binfo1 =  self.findcon_FIR(order,cond_ind,startbin);%all cons+bins before this con's bin need to be skipped
                    binfo2 = self.findcon_FIR(order,cond_ind,startbin+binwin-1);
                    vec_pure = vec0;
                    vec_pure(binfo1:binfo2) = 1;
                    %now add all nuis
                    vec = [];
                    for nr = realrun(:)'
                        vec =  [vec vec_pure zeros(1,(Nmov+Ntsda(nr)))];
                    end
                    %add session constants
                    if self.id == 15
                        nsessions = 1;
                    end
                    vec = padarray(vec,[0 nsessions],0,'post'); %his is done to compare it to Nbetas later. Not necessary per se
                    %check Nbeta
                    if ~isequal(Nbetas,length(vec))
                        keyboard;
                    end
                    convec{n} = vec;
                    name{n} = sprintf('bin%02d.win%02d.%s',startbin,binwin,condnames{cond_ind});
                end
                %
                cs_ind = [4 8; 2 1; 4 8];
                n = n+1;
                vec_pure = vec0;
                csp_bins =  self.findcon_FIR(order,cs_ind(nrun,1),startbin):self.findcon_FIR(order,cs_ind(nrun,1),startbin+binwin-1); %Nth bin beta of CSP
                csn_bins =  self.findcon_FIR(order,cs_ind(nrun,2),startbin):self.findcon_FIR(order,cs_ind(nrun,2),startbin+binwin-1); %Nth bin beta of CSN
                vec_pure(csp_bins) =  1;
                vec_pure(csn_bins) = -1;
                %now add all nuis
                vec = [];
                for nr = realrun(:)'
                    vec =  [vec vec_pure zeros(1,(Nmov+Ntsda(nr)))];
                end
                if self.id == 15
                    nsessions = 1;
                end
                vec = padarray(vec,[0 nsessions],0,'post');
                convec{n} = vec;
                name{n} = sprintf('bin%02d.win%02d.CSdiff',startbin,binwin);
                
                %make convecs sum to 1 and -1 for t-contrasts
                for co = (n2+1):n
                    convec{co} = convec{co}./binwin/nsessions;
                end
            end
            
            %visualization and sanity check
            figure(100);
            if nrun == 1
                clf
            end
            subplot(1,3,nrun)
            for n = 1:numel(convec)
                plot(convec{n}+n)
                hold on
            end
            for c = 1:max(conds2take{nrun})
                l=line(repmat((order*(c-1)+1),1,2),ylim);set(l,'Color','k','LineStyle',':')%line on first bin of each condition
            end
            l=line(repmat(max(conds2take{nrun})*order+1,1,2),ylim);set(l,'Color','k','LineWidth',1);%line where contrasts stop being relevant
            l=line(repmat(Nbetas./nsessions-1,1,2),ylim);set(l,'Color','k','LineWidth',2); %end of vec for 1 session
            
            set(gca,'XTick',1:order:order*max(conds2take{nrun}),'XTickLabel',condnames,'YTick',1:order:numel(convec)-1,'YTickLabel',name(1:order:numel(convec)-1))
            title(sprintf('Run %d',nrun));
            %%
            % build the batch
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(path_spm);
            for co = 1:numel(name)
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.name    = name{co};
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec  = convec{co};
                matlabbatch{1}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
                %sanity check
                if length(convec{co})~= Nbetas
                    fprintf('Problem with length of convec at contrast %d. Please debug.\n',co);
                    keyboard
                end
            end
            total_cons = numel(name);
            
            matlabbatch{1}.spm.stats.con.delete = 1;
            
            spm_jobman('run',matlabbatch);
            %             for nrun = runs(:)'
            %                 fprintf('\n\nTook %05.2f minutes for 1stLevel contrasts all conds FIR model 0, sub %02d to %02d, run %02d.\n',tocc(nrun)./60,subs(1),subs(end),nrun);
            %             end
            %             fprintf('\n\nTook %05.2f minutes (%05.2f hours) for 1stLevel contrasts all conds FIR model 0, sub %02d to %02d, run %02d to %02d.\n',done./60,done./3600,subs(1),subs(end),runs(1),runs(end));
            
        end
        function Con1stLevel(self,nrun,model_num)
            n_con = self.CreateContrasts(nrun,model_num);
            fprintf('%d contrasts computed for phase %d, Model %d. \n',n_con,nrun,model_num)
            %prepare paths for Normalization and Smoothing
            path2cons = self.path_con(nrun,model_num,'');
            self.VolumeNormalize(path2cons);%normalize ('wCAT_' and 'wEPI_'' will be added)
            self.VolumeSmooth(self.path_con(nrun,model_num,'wCAT_'));% smooth (s6, or whatever kernel, will be added)
            path2wCAT = self.path_con(nrun,model_num,'wCAT_');
            %remove the unsmoothed image for storage reasons
            self.remove_files(path2wCAT);
        end
        
        function [firstcon,lastcon] = Con1stLevel_FIR(self,nrun,model_num,varargin)
            deletedcons = 1;
            % we need to know which cons are new, so that we can normalize and smooth only them, bc this takes time.
            if deletedcons == 1
                firstcon = 1; %start at the first, go till end.
            else
                load(fullfile(self.path_FIR(nrun,model_num,self.orderfir,'10conds'),'SPM.mat'))
                firstcon = numel(SPM.xCon)+1;
            end
            
            if nargin > 3
                startbin = varargin{1};
                binwin   = varargin{2};
                if model_num == 40
                    n_con = self.CreateContrasts_FIR_tsda(nrun,model_num,startbin,binwin);
                else
                    n_con = self.CreateContrasts_FIR(nrun,model_num,startbin,binwin);
                end
            else
                if model_num == 40
                    n_con = self.CreateContrasts_FIR_tsda(nrun,model_num);
                else
                    n_con = self.CreateContrasts_FIR(nrun,model_num);
                end
            end
            if deletedcons == 1
                lastcon = n_con; %start at the first, go till end.
            else
                lastcon = firstcon+n_con;
            end
            
            fprintf('%d contrasts computed for phase %d, Model %d. \n',n_con,nrun,model_num)
            %prepare paths for Normalization and Smoothing
            for co = firstcon:lastcon
                pathcon  = self.path_FIR(nrun,model_num,self.orderfir,'10conds',sprintf('^con_%04d.nii',co));
                %normalize con image
                self.VolumeNormalize(pathcon);
                % smooth con image
                self.VolumeSmooth(self.path_FIR(nrun,model_num,self.orderfir,'10conds',sprintf('^wCAT_con_%04d.nii',co)));
                %remove the unsmoothed image for storage reasons
                self.remove_files(self.path_FIR(nrun,model_num,self.orderfir,'10conds',sprintf('^wCAT_con_%04d.nii',co)));
            end
        end
        %         function plot_con(nrun,model_num,con_num)
        %             self.CreateContrast(nrun,model_num,con_num);
        %         end
        
        function remove_files(self,allfiles)
            for n = 1:size(allfiles)
                system(sprintf('rm %s',allfiles(n,:)));
            end
        end
        function out = listcond(self,nrun,model_num)
            condfile = load(self.path_model(nrun,model_num));
            cond = condfile.cond;
            out = [];
            for n= 1:length(cond);
                out = char(out,cond(n).name);
            end
            out(1,:) = [];
            out = cellstr(out);
        end
        function t = checkcond(self,nrun,model_num)
            condfile = load(self.path_model(nrun,model_num));
            cond = condfile.cond;
            names = self.listcond(nrun,model_num);
            onset=[];
            dur = [];
            for ncond = 1:length(cond)
                onset(ncond) = length(cond(ncond).onset);
                dur(ncond) = mean(cond(ncond).duration)*self.TR;
            end
            
            t = table(onset',dur','RowNames',names);
            t.Properties.VariableNames = {'onset','dur_secs'};
        end
        
        function log2designmat(self,nrun,model_num)
            load(self.path_spmmat(nrun,model_num))
            pfile = self.get_paradigm(nrun);
            face_onsets_inds = pfile.out.log(:,2)==13;
            timepoints = round(pfile.out.log(face_onsets_inds,1));
            duration = ceil(pfile.out.log(end,1));
            designmat = zeros(duration,self.nreliefconds(nrun));
            
            triallist = pfile.presentation.dist;
            colind = self.compute_deltacsp2ind(triallist);
            colind(1) = 10;
            colind(triallist==500) = 9;
            
            for n =1:10
                designmat(timepoints(colind==n),n)=1;
            end
            designmat(1:(floor(min(timepoints))-10),:)=[];
            designmat(:,sum(designmat)==0)=[];
            figure;
            subplot(1,2,1);imagesc(designmat)
            subplot(1,2,2);imagesc(SPM.xX.X(:,1:self.nreliefconds(nrun)));
        end
        
    end
end
