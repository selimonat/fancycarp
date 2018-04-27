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
    end
    %%
    methods
        function s = Subject(id)%constructor
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
    end
    methods %(behavioral analysis)
        function [out, raw, triallist] = get_rating(self,varargin)
            % this function is mostly to get ratings so that we can fit
            % it.. so it puts it to out.x and out.y, and also only for
            % p.realconds...
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
                    catch
                        warning('Problem while loading paradigm at run %d.',run)
                        raw{run} = NaN;
                        triallist{run} = NaN;
                    end
                elseif run == 5
                    a = self.get_paradigm(3);
                    raw{run} = a.log.ratings.relief(:,3);
                    triallist{run} = a.presentation.dist(:);
                    try
                        a = self.get_paradigm(4);
                        raw{run} = [raw{run}; a.log.ratings.relief(:,3)];
                        triallist{run} = [triallist{run}; a.presentation.dist'];
                    end
                end
            end
            
            for r = 1:length(raw)
                out{r}.x = [];
                out{r}.y = [];
                out{r}.ids = [];
                cc = 0;
                for cond = conds
                    cc = cc+1;
                    ind = triallist{r} == cond;
                    out{r}.x = [out{r}.x repmat(cond,1,sum(ind))];
                    out{r}.y = [out{r}.y raw{r}(ind)'];
                    out{r}.ids = self.id;
                end
            end
            
            if numel(runs) == 1
                out = out{end};
                raw = raw{end};
                triallist = triallist{end};
            end
        end
        
        function out = get_pain(self,varargin)
            if isempty(varargin)
                fprintf('No run selected, collecting all runs.\n')
                for run = 1:4
                    try
                        a = self.get_paradigm(run);
                        dummy = a.log.ratings.pain;
                        dummy = dummy(~isnan(dummy(:,3)),3);
                        if run < 3
                            out(:,run) = [dummy(1:3); NaN];
                        else
                            out(:,run) = dummy(1:4);
                        end
                    catch
                        warning('Problem while loading paradigm at run %d.',run)
                        out(:,run) = nan(4,1);
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
                        dummy = dummy(~isnan(dummy(:,3)),3);
                        if run < 3
                            out = dummy(1:3);
                        else
                            out = dummy(1:4);
                        end
                    catch
                        warning('Problem while loading paradigm at run %d.',run)
                        out = NaN;
                    end
                end
            end
            
        end
        function [M, S] = get_reliefmeans(self,run,varargin)
            conds = self.allconds;
            if nargin > 2
                conds = varargin{1};
            end
            out = self.get_rating(run,conds);
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
        function [out, R] = fit_rating(self,run)
            %will load the rating fit (saved in runXXX/rating) if computed other
            %wise will read the raw ratingdata (saved in runXXX/stimulation)
            %and compute a fit.
            
            fun        = self.selected_fitfun;%vM function
            force      = 0;%repeat the analysis or load from cache
            write_path = sprintf('%s/midlevel/rating_fun_%i.mat',self.pathfinder(self.id,run),fun);
            
            if exist(write_path) && force ==0
                %load directly or
                load(write_path);
                fprintf('Rating Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif force == 1 || ~exist(write_path)
                %compute and save it.
                fprintf('Fitting Rating for run %02d..\n',run)
                R = self.get_rating(run);
                if isempty(R.y)
                    warning('No ratings, so no fit for this person and run.\n');
                    out = [];
                else
                    %create a tuning object and make a single subject fit
                    T                            = Tuning(R);
                    T.SingleSubjectFit(fun);
                    %prepare data for outputting.
                    out.params                   = T.fit_results.params(1,:);
                    out.LL                       = T.fit_results.pval(1,:);
                    out.pval                     = 10.^-T.fit_results.pval(1,:);
                    out.exitflag                 = T.fit_results.ExitFlag(1,:);
                    out.y_hd                     = T.fit_results.y_fitted_HD(1,:);
                    out.x_hd                     = T.fit_results.x_HD(1,:);
                    out.y                        = T.fit_results.y_fitted(1,:);
                    out.x                        = T.fit_results.x(1,:);
                    out.fitfun                   = T.fit_results.fitfun;
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
            
        end
        function plot_ratings(self)
            force = 1;
            savepath = [self.path 'run005/figures/'];
            if ~exist(savepath)
                mkdir(savepath)
            end
            savefigf = [savepath sprintf('ratingplot_sub%02d.fig',self.id)];
            savebmp = [savepath sprintf('ratingplot_sub%02d.bmp',self.id)];
            savepng = [savepath sprintf('ratingplot_sub%02d.png',self.id)];
            if exist(savefigf,'file') && force == 0
                fprintf('Found saved .fig file, opening it as figure %02d\n',self.id)
                openfig(savefigf);
            else
                
                titles = {'Base','Cond','Test1','Test2','Test1/2'};
                f=figure(self.id);
                f.Position = [73 410 1763 476];
                clf;
                
                for run = 1:5
                    [M,S] = self.get_reliefmeans(run,self.allconds);
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
%                 export_fig(gcf,savebmp)
%                 export_fig(gcf,savepng,'-transparent')
            end
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
            %             %% Normalize with CAT12 segmentation
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
            %will load the pmf fit (saved in runXXX/pmf) if computed other
            %wise will read the raw pmf data (saved in runXXX/stimulation)
            %and compute a fit.
            
            if exist(self.path_data(2,'pmf')) && isempty(varargin)
                %load directly or
                load(self.path_data(2,'pmf'));
                %                 fprintf('PMF Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif ~isempty(varargin) || ~exist(self.path_data(2,'pmf'))
                addpath(self.path_palamedes);
                %compute and save it.
                fprintf('Fitting PMF...\n')
                % define a search grid
                searchGrid.alpha  = linspace(0,100,10);    %structure defining grid to
                searchGrid.beta   = 10.^[-1:0.1:1];         %search for initial values
                searchGrid.gamma  = linspace(0,0.5,10);
                searchGrid.lambda = linspace(0,0.1,10);
                paramsFree        = [1 1 1 1];
                PF                = @PAL_Weibull;
                %prepare some variables
                tchain            = size(self.pmf.log.xrounded,3);
                xlevels           = unique(abs(self.pmf.presentation.uniquex));
                NumPos            = NaN(length(xlevels),tchain);
                OutOfNum          = NaN(length(xlevels),tchain);
                sd                = NaN(length(xlevels),tchain);
                %first collapse the two directions (pos/neg differences from
                %csp)
                
                for chain = 1:tchain
                    fprintf('Starting to fit chain %g...\n',chain)
                    %get responses, and resulting PMF from PAL algorithm
                    data = self.pmf.log.xrounded(:,:,chain);
                    rep  = self.pmf.presentation.rep;
                    cl   = 0;
                    for l = xlevels(:)'
                        cl                 = cl+1;
                        ind                = find(abs(self.pmf.presentation.uniquex) == l);
                        collecttrials      = data(ind,1:rep(ind(1)));
                        collecttrials      = collecttrials(:);
                        NumPos(cl,chain)   = sum(collecttrials);% number of "different" responses
                        OutOfNum(cl,chain) = length(collecttrials);%number of presentations at that level
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
                end
                save(self.path_data(2,'pmf'),'out')
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
            verbalize = 1;
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
                        
                        conds = [13 4 11 9 18];% 4 = RampOnset 11 = RateTreat, 9=RatePain
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
                    modelname = 'RampOnsetStick_wmcsf'; %based on 6, onsets are RampOnsets, but we model RampOnsets as stick (but everything else as box)
                    if ~exist(fileparts(modelpath));mkdir(fileparts(modelpath));end
                    if verbalize == 1
                        fprintf('Preparing condition file for phase: ...................... %s.\n',Project.plottitles{run})
                    end
                    if ~exist(modelpath) || force == 1
                        load(self.path_model(run,6))
                        save(modelpath,'cond')
                    else
                        fprintf('Loading cond-mat file from modelpath %s.\n',modelpath)
                        load(modelpath);
                    end
            end
        end
        function out = get_pmod(self,nrun,cond)
            if cond == 999 %turning thermode off in the end -  no relief rating but pain before and after turning off.
                if self.id == 4
                    finalrelief = 0;
                elseif self.id == 15 && nrun == 4;
                    finalrelief = NaN;
                else
                    finalrelief = -diff(self.get_paradigm(nrun).log.ratings.pain(end-1:end,3)); %diff of painratings when thermode is on vs off in the end
                end
                out = struct('name',{'999'},'param',{finalrelief},'poly',1);
            else
                out = struct('name',{num2str(cond)},'param',{self.get_rating(nrun,cond).y},'poly',1);
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
                
                c2_templ = '^c2data.nii';
                c3_templ = '^c3data.nii';
                
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
                files = spm_select('ExtFPList', rundir, '^r.*.nii',1:nscans);
                
                
                V = spm_vol(files);
                y = spm_read_vols(V);
                
                %this is different from AT's script on segmented epis, bc. diff resolution in T1 and
                %epis
                for j = 1%:size(files,1)
                    wmdummypath  = fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_wm_dummy%d.nii',j));
                    csfdummypath =  fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_csf_dummy%d.nii',j));
                    spm_imcalc([y1,y(:,:,:,j)],wmdummypath,'i1*i2');
                    
                    spm_imcalc([y2,y(:,:,:,j)],csfdummypath,'i1*i2');
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
        function [kickcooldownopts, wmcsfopts] = get_model_specs(self,modelnum)
            kickcooldownopts = self.kickcooldown;
            wmcsfopts = self.wmcsfregressors;
            switch modelnum
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
            end
        end
        function FitFIR(self,nrun,model_num)
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            empty1stlevel  =   1;
            FIRparam       =  14;
            addpmod        =   0;
            onset_modelnum =   3;
            
            if model_num == 3
                addpmod =1;
            end
            
            spm_dir  = strrep(self.dir_spmmat(nrun(1),model_num),'chrf_00',sprintf('FIR_%02d_10conds_00',FIRparam));
            path_spmmat = fullfile(spm_dir,'SPM.mat');
            
            if ~exist(path_spmmat);mkdir(spm_dir);end
            
            if empty1stlevel == 1
                if exist(spm_dir); system(sprintf('rm -fr %s*',strrep(spm_dir,'//','/')));end %this is AG style.
            end
            
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
                    dummy.cond = dummy.cond(1:self.nreliefconds(session)); % 8 faces + UCS + t0 (if there)
                end
                
                for c = 1:numel(dummy.cond)
                    dummy.cond(c).duration  = 0;
                    dummy.cond(c).onset     = dummy.cond(c).onset -2; %take 2 TRs before first stimulus appears. (onsets are already in unit = TR, so -2 is correct, not -2.*TR).
                    if addpmod == 1
                        if str2double(dummy.cond(c).name)< 900 %don't do it for nulltrial and finalramp, we only have 1 rating here.
                            dummy.cond(c).pmod = self.get_pmod(session,str2double(dummy.cond(c).name));
                        end
                    end
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
            [kickcd,wmcsfr] = self.get_model_specs(model_num);
          
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
                dummy                                              = load(self.path_model(session,model_num));
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
                    load(fullfile(self.pathfinder(self.id,nrun),'midlevel',sprintf('wmcsf_reg_kc%d.mat',kickcd)),'wm','csf');
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+1).val   = wm';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+1).name  = 'wm';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+2).val   = csf';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(se).regress(nNuis+2).name  = 'csf';
                end
                %
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
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'AR(1)';
            %estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat            = {path_spmmat};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
            %
            %normalize the beta images right away
            beta_images = self.path_beta(nrun(1),model_num,'');%'' => with no prefix
            self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)
            beta_images = self.path_beta(nrun(1),model_num,'wCAT_');%smooth the normalized images.
            self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_ww_')
            %delete wCAT files, keep s6_wCAT
            fprintf('Deleting unsmoothed normalized files..\n');
            for n = 1:size(beta_images,1)
                system(sprintf('rm %s',beta_images(n,:)));
            end
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
%                 Nbetas = self.get_Nbetas(nrun,model_num);
%                 if nrun == 3
%                     if self.id ~= 15
%                         Nbetas = Nbetas./2; %to care for one session first, not necessary for sub 15, as only one session modelled anyway.
%                     end
%                 end
%                 % contrast 1
%                 % all faces main effect
%                 n = n + 1;
%                 name{n} = 'allfaces>else';
%                 face_betas = [];self.path_epi(session,'r')
%                 convec{n} = zeros(1,Nbetas);
%                 for cond = intersect(self.realconds,unique(self.get_paradigm(nrun).presentation.dist))
%                     ind = self.get_beta_index(nrun,model_num,[num2str(cond) 'Face']);
%                     face_betas = [face_betas ind];
%                 end
%                 convec{n}(face_betas) = deal(1);
%                 
%                 % contrast 2
%                 % CSP > CSN (or UCS > CSN in Cond)
%                 % given by FACES
%                 n = n + 1;
%                 convec{n} = zeros(1,Nbetas);
%                 if nrun == 2
%                     convec{n}(self.get_beta_index(nrun,model_num,'500Face'))= 1;
%                     convec{n}(self.get_beta_index(nrun,model_num,'180Face'))= -1;
%                     name{n} = 'UCS>CSN_faces';
%                 else
%                     convec{n}(self.get_beta_index(nrun,model_num,'0Face'))= 1;
%                     convec{n}(self.get_beta_index(nrun,model_num,'180Face'))= -1;
%                     name{n} = 'CSP>CSN_faces';
%                 end
%                 
%                 % contrast 3
%                 % all Relief trials vs else
%                 n = n + 1;
%                 convec{n} = zeros(1,Nbetas);
%                 name{n} = 'allrampdown>else';
%                 thisphaseconds = intersect(self.realconds,unique(self.get_paradigm(nrun).presentation.dist));
%                 convec{n}(self.get_beta_index(nrun, model_num,thisphaseconds)) = 1;
%                 
%                 % contrast 4
%                 % CSP > CSN (or UCS > CSN in Cond)
%                 % given by RAMP
%                 n = n + 1;
%                 convec{n} = zeros(1,Nbetas);
%                 if nrun == 2self.path_epi(session,'r')
%                     convec{n}(self.get_beta_index(nrun,model_num,[500 180]))= [1 -1];
%                     name{n} = 'UCS>CSN_ramp';
%                 else
%                     convec{n}(self.get_beta_index(nrun,model_num,[0 180]))= [1 -1];
%                     name{n} = 'CSP>CSN_ramp';
%                 end
%                 % take care of the two sessions
%                 if nrun ==3 && self.id ~=15
%                     for co = 1:n
%                         convec{co} = repmat(convec{co},1,2);
%                     end
%                 end
                
            elseif ismember(model_num,[7 8])
                load([self.dir_spmmat(1,model_num) 'SPM.mat'])
                for pmod = 1:size(SPM.Sess.U(1).P,2) %loop through the pmods we defined
                    
                    vec = zeros(1,self.get_Nbetas(nrun,model_num)-1);
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
            designmat = zeros(5-,emt[duration,self.nreliefconds(nrun));
            
            triallist = pfile.presentation.dist;
            colind = self.compute_deltacsp2ind(triallist);
            colind(1) = 9;
            
            for n =1:9
                designmat(timepoints(colind==n),n)=1;
            end
            designmat(1:(floor(min(timepoints))-10),:)=[];
            figure;
            subplot(1,2,1);imagesc(designmat)
            subplot(1,2,2);imagesc(SPM.xX.X(:,1:self.nreliefconds(nrun)));
        end
    end
end
