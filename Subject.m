classdef Subject < Project
    
    properties (Hidden)
        paradigm
        default_run       = 1;   
        dicom_serie_id    = [];
        dicom_folders     = [];
        dicom_target_run  = [];
        
        pmf               = []; %raw results of the pmf
        pmf_variablenames = {'alpha_chain01','beta_chain01','gamma_chain01','lapse_chain01','alpha_chain02','beta_chain02','gamma_chain02','lapse_chain02'}
    end
    properties (SetAccess = private)
        id
        path
        csp
        csn
        scr        
        pmf_parameters= [];
        trio_session  = [];
        ratings       = [];
        total_run     = [];
    end
    %%
    methods
        function s = Subject(id)%constructor
            s.id               = id;
            s.path             = s.pathfinder(s.id,[]);
            s.dicom_serie_id   = s.dicom_serie_selector{s.id};
            s.dicom_target_run = s.dicom2run{s.id};
            try
                s.trio_session 	  = s.trio_sessions{s.id};
            end
            
            if exist(s.path)
                for nrun = 1:s.total_run
                    s.paradigm{nrun} = s.load_paradigm(nrun);
                end
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn = s.paradigm{s.default_run}.stim.cs_neg;
                s.scr = SCR(s);                
                
                %                 try
                %                     s.bold = BOLD(s);
                %                 end
            else
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path);
                fprintf('Your path might also be wrong...\n');
            end
        end
    end
    
    methods %(Getters)
        
        function dump_hr(self)
            %will download the latest HR for this subject to default hr
            %path.
            %
            %
            %
            
            %target location for the hr
            hr_target = self.hr_dir;
            %create it if necess.
            if exist(hr_target) == 0
                mkdir(hr_target);
            end
            self.DicomDownload(self.GetDicomHRpath,self.hr_dir);
        end
        function p          = load_paradigm(self,nrun)
            filename = self.path2data(nrun,'stimulation');
            p = [];
            if exist(filename)
                p = load(filename);
                p = p.p;
                %transform id to labels
                if isfield(p,'presentation')
                    %                 p.presentation.stim_label = self.condition_labels(p.presentation.cond_id+1);
                    p.presentation.dist(p.presentation.dist < 500)   = p.presentation.dist(p.presentation.dist < 500) + 180;
                    p.presentation.dist(p.presentation.dist == 500)  = 1001;%ucs
                    p.presentation.dist(p.presentation.dist == 1000) = 1002;%odd
                    p.presentation.dist(isnan(p.presentation.dist))  = 1000;%null
                end
            end
        end 
        function dump_functional(self)
            %Will dump all DICOMS based on Sessions entered in the
            %Project Object. trio_folders are folders in the dicom server,
            %trio2run dictates in which run these folders should be dumped.
            %
            
            %spit out some info for sanity checks
            self.dicomserver_request;
            fprintf('You told me to download the following series: ');
            fprintf('%i,',self.dicom_serie_id);
            fprintf('\nDouble check if everything is fine.\n')
            
            paths               = self.dicomserver_paths;
            self.dicom_folders  = paths(self.dicom_serie_id);
            fprintf('Will now dump series (%s)\n',self.gettime);            
            
            %% save the desired runs to disk
            n = 0;
            for source = self.dicom_folders(:)'
                %
                n 				 = n+1;
                dest             = sprintf('%ssub%03d/run%03d/mrt/',self.path_project,self.id,self.dicom_target_run(n));                
                self.DicomDownload(source{1},dest);
            end
            fprintf('Happily finished dumping...%s\n',self.gettime);
        end
        function rating = get.ratings(self)
            %returns the CS+-aligned ratings for all the runs
            for run = 1:self.total_run;%don't count the first run
                if isfield(self.paradigm{run}.out,'rating')
                    if ~isempty(self.paradigm{run});
                        rating(run).y      = self.paradigm{run}.out.rating';
                        rating(run).y      = circshift(rating.y,[1 4-self.csp ]);
                        rating(run).x      = repmat([-135:45:180],size(self.paradigm{run}.out.rating,2),1);
                        rating(run).i      = repmat(run          ,size(self.paradigm{run}.out.rating,2),size(self.paradigm{run}.out.rating,1));
                        rating(run).y_mean = mean(rating.y);
                    else
                        fprintf('No rating present for this subject and run (%d) \n',nr);
                    end
                end
            end
        end
        function out    = GetSubSCR(self,run,cond)
            if nargin < 3
                cond=1:8;
            end
            conddummy = [-135:45:180 500 1000 3000];
            % s is a subject instance
            out       = [];
            cutnum    = self.scr.findphase(run);
            self.scr.cut(cutnum);
            self.scr.run_ledalab;
            self.scr.plot_tuning_ledalab(cond);
            out.y     = self.scr.fear_tuning;
            out.x     = conddummy(cond);
            out.ind   = cutnum;
        end        
        function [o]    = get.total_run(self)
            %% returns the total number of runs in a folder (except run000)
            o      = length(dir(self.path))-3;%exclude the directories .., ., and run000
        end
        
        function out    = get_log(self,run)
            %loads and plots the Log, here you can put the old experiment
            %plotter.
            
            out        = self.paradigm{run}.out.log;
            %sort things according to time rather than order of being
            %logged
            [~,i]      = sort(out(:,1),'ascend');
            out        = out(i,:);
            % delete all the events that are after the last scanning..
            scan_times = out(out(:,2) == 0,1);
            i          = out(:,1) > max(scan_times);            
        end
        function out    = get.pmf(self)
            %will load the raw pmf data.
            dummy    = load(self.path2data(2,'stimulation'));
            out      = dummy.p.psi;
        end 
        function out    = get.pmf_parameters(self)
            %returns the parameters of the pmf fit (condition x parameter);
            out      = self.fit_pmf;
            out      = [out.params(1,:),out.params(2,:)];
            out      = array2table([out self.id ],'variablenames',[self.pmf_variablenames 'subject_id']);
        end 
    end
    
    methods %(mri, preprocessing))      
        function ConvertDicom(self)
            %% dicom conversion. ATTENTION: dicoms will be deleted
            % and the converted files will be merged to a 4d file nifti
            % file. This file will be named data.nii.
            %
            matlabbatch = [];
            count = 0;
            for nrun = 0:self.tRuns
                files                                                    = spm_select('FPListRec',self.path2data(nrun),'^MR');
                fprintf('Run#%d found %d dicom files.\n',nrun,length(files))
                if ~isempty(files)%only create a batch if there is ^MR files.
                    count = count +1;
                    matlabbatch{count}.spm.util.import.dicom.data             = cellstr(files);
                    matlabbatch{count}.spm.util.import.dicom.root             = 'flat';
                    matlabbatch{count}.spm.util.import.dicom.outdir           = {fileparts(self.path2data(nrun,'mrt'))};
                    matlabbatch{count}.spm.util.import.dicom.protfilter       = '.*';
                    matlabbatch{count}.spm.util.import.dicom.convopts.format  = 'nii';
                    matlabbatch{count}.spm.util.import.dicom.convopts.icedims = 0;
                end
            end
            %don't continue if there is nothing to do...
            if ~isempty(matlabbatch)
                fprintf('Dicom conversion s#%i... (%s)\n',self.id,self.gettime);
                self.RunSPMJob(matlabbatch);
                fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));                                                                
                %% delete the 3d fTRIO files                
                fprintf('Deleting all 3D and DICOM images s#%i... (%s)\n',self.id,self.gettime);
                for nrun = 1:length(matlabbatch)
                    delete(sprintf('%smrt%sMR*',self.path2data(nrun),filesep));
                    delete( matlabbatch{nrun}.spm.util.cat.vols{:} );
                end
                fprintf('Finished... (%s)\n',self.gettime);
            else
                fprintf('No dicom files found for %i\n',self.id)
            end
        end
        function MergeTo4D(self)
            %will create data.nii consisting of all the [f,s]TRIO images
            %merged to 4D.
            
            %% merge to 4D
            fprintf('Merging s#%i...(%s)\n',self.id,self.gettime);
            matlabbatch = [];
            c           = 0;
            for nrun = 0:self.tRuns
                files = spm_select('FPListRec',self.path2data(nrun),'^[f,s]TRIO');
                if ~isempty(files)
                    c                                 = c + 1;
                    matlabbatch{c}.spm.util.cat.vols  = cellstr(files);
                    matlabbatch{c}.spm.util.cat.name  = 'data.nii';
                    matlabbatch{c}.spm.util.cat.dtype = 0;
                end
            end
            if ~isempty(matlabbatch)
                self.RunSPMJob(matlabbatch);
            end
            fprintf('Finished... (%s)\n',self.gettime);
        end        
        function Realign(self)
            % Will realign all runs using spm batch.
            matlabbatch = [];
            for nrun  = 1:self.tRuns
                file = sprintf('%smrt%sdata.nii',self.path2data(nrun),filesep);
            end
            matlabbatch{1}.spm.spatial.realign.estwrite.data             = {cellstr(file)};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight  = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which   = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp  = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask    = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
            fprintf('Realigning s#%i...(%s)\n',self.id,datestr(now,'hh:mm:ss'));
            spm_jobman('run',matlabbatch);
            fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
        end
        function Coreg_Anat2Functional(self)
            %rigid-body wiggles around the anatomical to the mean EPI from
            %the realignment procedure.
            mean_file =sprintf('%smrt%smeandata.nii',self.path2data(1),filesep);%is mean_file always saved to the same place?
            if exist(mean_file) > 0
                matlabbatch{1}.spm.spatial.coreg.estwrite.ref                = { mean_file };%the one that stays constant
                matlabbatch{1}.spm.spatial.coreg.estwrite.source             = { self.hr_path };%anatomical one.
                matlabbatch{1}.spm.spatial.coreg.estwrite.other              = {''};
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun  = 'nmi';
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep       = [4 2];
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol       = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm      = [7 7];
                %
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp    = 6;
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap      = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask      = 0;
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix    = 'r';
                spm_jobman('run',matlabbatch);
            else
                fprintf('I think you need to realign first, the file:\n %s \n is not yet computed...\n',mean_file);
            end
        end        
        function o = MotionParameters(self,run)
            %will load the realignment parameters.
            filename = sprintf('%smrt/rp_data.txt',self.path2data(run));
            if exist(filename)
                o = load(filename);
            else
                fprintf('File:\n %s doesn''t exist.\n Most likely realignment is not yet done.\n');
            end
        end        
        function [scanunit]=StimTime2ScanUnit(self,run)
            %will return stim onsets in units of scan. Usefull for first
            %level.
            
            L          = self.Log(run);
            scan_times = L(L(:,2) == 0,1);%find all scan events
            scan_id    = 1:length(scan_times);
            nstim      = 0;
            %ideally we could divide stim_onsets with TR, but this would
            %create a problem if there are reference measurements during
            %the measurements.
            for stim_times = L(find(L(:,2)==3),1)';%run stim by stim
                nstim          = nstim + 1;
                d              = scan_times - stim_times;
                first_positive = find(d > 0,1);
                decimal        = d(first_positive)./self.TR;
                scanunit(nstim)= (scan_id(first_positive)-1) + decimal;
            end
        end
    end
    methods %fmri path_tools
        function out        = spm_dir(self)
            %returns subject's path to spm folder for run RUN.
            out = sprintf('%smrt/spm/',self.pathfinder(self.id,1));
        end        
        function out        = spm_path(self)
            %returns the path to spm folder for run RUN.
            out = sprintf('%smrt/spm/SPM.mat',self.pathfinder(self.id,1));
        end        
        function out        = hr_dir(self)
            %the directory where hr is located
            out = sprintf('%smrt/',self.pathfinder(self.id,0));
        end        
        function out        = hr_path(self)
            %path to the hr volume
            out = spm_select('ExtFPList',self.hr_dir,'^sTRIO.*.nii$');
        end        
        function [t]        = total_volumes(self,run)
            % will tell you how many volumes are in a 4D image.
            bla = spm_vol_nifti(self.mrt_data(run),1);%simply read the first images header
            t   = bla.private.dat.dim(4);
        end        
        function out        = mrt_path(self,nrun)
            % simply returns the path to the mrt data.
            out = sprintf('%smrt/data.nii',self.pathfinder(self.id,1));
        end        
        function out        = mrt_path_expanded(self,nrun)
            %returns list of filenames of a 4D nii file using comma
            %separated convention (needed for First levels)            
            out = spm_select('ExtFPList',fileparts(self.mrt_data(nrun)),'^rdata.nii');
        end
        function [HRPath]   = GetDicomHRpath(self)
            % finds the dicom path to the latest HR measurement for this
            % subject.
            
            [status2 DicqOutputFull] = system(sprintf('/common/apps/bin/dicq --verbose  --series --exam=%s --folders',self.trio_session));
            %take the latest anatomical scan.
            [status2 HRLine] = system(sprintf('/common/apps/bin/dicq --verbose  --series --exam=%s --folders | grep mprage | tail -n 1',self.trio_session));
            %
            HRPath = [];
            if ~isempty(HRLine);
                HRPath = regexp(HRLine,'/common\S*','match');
                HRPath = HRPath{1};
                %HRPath = GetDicomPath(HRLine);
                fprintf('Dicom Server returns:\n=====\n')
                fprintf(DicqOutputFull);
                fprintf('=====\n');
                fprintf('The latest recorded HR data:\n')
                fprintf(HRLine);
            else
                fprintf('There is no HR data found for this subject.\n Here is the output of the Dicq:\n');
                fprintf(DicqOutputFull);
                
            end
            
        end
        function path2data  = path2data(self,run,varargin)
            % s.path2data(53,4) will return the path to the subject's phase 4
            % s.path2data(53,4,'eye') return the path to the eye data file at the
            % 4th phase.
            
            %will return the path to phase/data_type/
            path2data = self.pathfinder(self.id , run);
            if length(varargin) >= 1
                path2data = sprintf('%s%s%sdata.mat',path2data,varargin{1},filesep);
            end
            if length(varargin) == 2
                path2data = regexprep(path2data,'mat',varargin{2});
            end
        end       
    end
    methods %(fmri analysis)
        function FitFIR(self,nrun,model_num)
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            
            
            spm_dir = sprintf('%s/model_fir_%02d/',self.spm_dir,model_num);
            spm_path= sprintf('%s/model_fir_%02d/SPM.mat',self.spm_dir,model_num);
            
            if ~exist(self.spm_path);mkdir(spm_dir);end
            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            for session = nrun
                %load files using ...,1, ....,2 format
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans  = cellstr(self.mrt_data_expanded(session));
                %load the onsets
                dummy                                                   = load(sprintf('%sdesign/model%02d.mat',self.path2data(session),model_num));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = dummy.cond;
                %load nuissance parameters
                nuis                                                    = self.MotionParameters(nrun);
                nuis                                                    = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
                for nNuis = 1:size(nuis,2)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).val   = nuis(:,nNuis);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).name  = mat2str(nNuis);
                end
                %
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi               = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi_reg           = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).hpf                 = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                              = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length                  = self.TR*15;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order                   = 10;
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                           = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = {''};%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'none';
            %estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat            = {spm_path};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
        end
        function FitHRF(self,nrun,model_num)
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            
            
            spm_dir = sprintf('%s/model_chrf_%02d/',self.spm_dir,model_num);
            spm_path= sprintf('%s/model_chrf_%02d/SPM.mat',self.spm_dir,model_num);
            
            if ~exist(self.spm_path);mkdir(spm_dir);end
            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            for session = nrun
                %load files using ...,1, ....,2 format
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans  = cellstr(self.mrt_data_expanded(session));
                %load the onsets
                dummy                                                   = load(sprintf('%sdesign/model%02d.mat',self.path2data(session),model_num));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = dummy.cond;
                %load nuissance parameters
                nuis                                                    = self.MotionParameters(nrun);
                nuis                                                    = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
                for nNuis = 1:size(nuis,2)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).val   = nuis(:,nNuis);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).name  = mat2str(nNuis);
                end
                %
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi               = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi_reg           = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).hpf                 = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                              = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs                  = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                           = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = {''};%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'none';
            %estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat            = {spm_path};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
        end
    end
    methods %analysis
        function [out] = fit_pmf(self,varargin)
            %will load the pmf fit (saved in runXXX/pmf) if computed other
            %wise will read the raw pmf data (saved in runXXX/stimulation)
            %and compute a fit.            
                
            if exist(self.path2data(2,'pmf')) && isempty(varargin)
                %load directly or 
                load(self.path2data(2,'pmf'));
%                 fprintf('PMF Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif ~isempty(varargin) || ~exist(self.path2data(2,'pmf'))
                addpath('/Users/onat/Documents/Code/Matlab/palamedes1_8_0/Palamedes/');
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
                save(self.path2data(2,'pmf'),'out')
            end
        end
    end
    methods %(plotters)
        function plot_log(self,nrun)
            %will plot the events that are logged during the experiment.
            L       = self.get_log(nrun);
            tevents = size(L,1);
            figure;
            plot(L(1:tevents,1) - L(1,1),L(1:tevents,2),'o','markersize',10);
            ylim([-2 8]);
            set(gca,'ytick',[-2:8],'yticklabel',{'Rating On','Text','Pulse','Tracker+','Cross+','Stim+','CrossMov','UCS','Stim-','Key+','Tracker-'});
            grid on
            drawnow;
        end       
        function plot_pmf(self)
            
            try
                addpath(self.palamedes_path);
            catch
                fprintf('Add Palamedes to your path.\n')
                return
            end
            %figure            
            out            = self.fit_pmf;
            tchain         = size(self.pmf.presentation.x,2);
            xlevels        = unique(abs(self.pmf.presentation.uniquex));
            StimLevelsFine = [min(xlevels):(max(xlevels) - min(xlevels))./1000:max(xlevels)];
            
            %plot the Fit
            ccc = rand(1,3);
            for chain = 1:tchain
                subplot(tchain,1,chain)
                Fit = out.PF(out.params(chain,:),StimLevelsFine);
                errorbar(xlevels,out.PropCorrectData(chain,:),out.sd(chain,:),'k.','Markersize',40);
                set(gca,'Fontsize',12);
                hold on;
                plot(StimLevelsFine,Fit,'-','Linewidth',3,'color',ccc);
                legend('data point','Fit','location','southeast');
                legend boxoff
                box off
                xlabel('Delta Degree');
                ylabel('p(diff)');
            end
            %
            subplot(tchain,1,1)
            title(sprintf('Sub %d, CSP (%d), estimated alpha = %g, LL = %g',self.id,self.csp, out.params(1,1),out.LL(1)));
            subplot(tchain,1,2)
            title(sprintf('CSN (%d), estimated alpha = %g, LL = %g',self.csn,out.params(2,1),out.LL(2)));
        end
        
    end
end
