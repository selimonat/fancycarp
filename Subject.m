classdef Subject < Project
    %   Getter Methods are used to get stuff e.g. ratings, epis, stimulation
    %   log, behavioral data etc. It also implements downloading data from the
    %   dicom server.
    %
    %   path_tools Methods are used to return the path to different data
    %   types, not really interesting, used internally.
    %
    %   plotter Methods plot subject specific data.
    %
    %   fmri analysis methods contain first-level modelling.
    %
    % 
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
            fprintf('Subject Constructor for id:%i is called:\n',id);
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
            %path (which is run000)
            %
            %
            %
            
            fprintf('dump_hr:\nWill now dump the latest HR (%s)\n',self.gettime);            
            %target location for the hr
            hr_target = self.hr_dir;
            %create it if necess.
            if exist(hr_target) == 0
                mkdir(hr_target);
            end
%             self.DicomDownload(self.GetDicomHRpath,self.hr_dir);
            self.DicomTo4D(self.hr_dir);
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
                dest             = self.epi_dir(self.dicom_target_run(n));                
                self.DicomDownload(source{1},dest);
            end
            self.DicomTo4D(dest);
        end
        function rating = get.ratings(self)
            %returns the CS+-aligned ratings for all the runs
            for run = 1:self.total_run;%don't count the first run
                if isfield(self.paradigm{run}.out,'rating')
                    if ~isempty(self.paradigm{run});
                        rating(run).y      = self.paradigm{run}.out.rating';
                        rating(run).y      = circshift(rating.y,[1 4-self.csp ]);
                        rating(run).x      = repmat([-135:45:180],size(self.paradigm{run}.out.rating,2),1);
                        rating(run).ids    = repmat(self.id,size(self.paradigm{run}.out.rating,2),size(self.paradigm{run}.out.rating,1));
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
        function o      = GetMotionParameters(self,run)
            %will load the realignment parameters, of course you have to
            %realign the EPIs first.
            filename = sprintf('%smrt/rp_data.txt',self.path2data(run));
            if exist(filename)
                o = load(filename);
            else
                fprintf('File:\n %s doesn''t exist.\n Most likely realignment is not yet done.\n');
            end            
        end
    end
    
    methods %(mri, preprocessing))      
        
        function preprocess_pipeline(self,run)
            %meta method to run all the required steps for hr
            %preprocessing.
            self.segment;
            self.SkullStrip;
            self.Dartel;
            self.Re_Coreg(run);
        end
        function segment(self)
            %take run000/mrt/data.nii and produces {r}c{1,2}data.nii and
            %data_seg8.mat files.
            
            matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(self.hr_path);
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[self.tpm_dir 'TPM.nii,1']};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[self.tpm_dir 'TPM.nii,2']};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[self.tpm_dir 'TPM.nii,3']};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[self.tpm_dir 'TPM.nii,4']};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[self.tpm_dir 'TPM.nii,5']};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[self.tpm_dir 'TPM.nii,6']};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
            %
            self.RunSPMJob(matlabbatch);
        end        
        function SkullStrip(self)
            %needs results of segment, will produce a skullstripped version
            %of hr.
            
            c1         = regexprep(self.hr_path,'mrt/data','mrt/c1data');
            c2         = regexprep(self.hr_path,'mrt/data','mrt/c2data');                
                        
            if exist(c1) && exist(c2)
                matlabbatch{1}.spm.util.imcalc.input            = {self.hr_path,c1,c2};
                matlabbatch{1}.spm.util.imcalc.output           = self.skullstrip;
                matlabbatch{1}.spm.util.imcalc.outdir           = {self.hr_dir};
                matlabbatch{1}.spm.util.imcalc.expression       = 'i1.*((i2+i3)>0.2)';
                matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask     = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp   = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype    = 4;
                self.RunSPMJob(matlabbatch);
            else
                fprintf('Need to run segment first...\n')
            end
        end        
        function Dartel(self)
            %will write u_rc1 file.
            
            rc1         = regexprep(self.hr_path,'mrt/data','mrt/rc1data');
            rc2         = regexprep(self.hr_path,'mrt/data','mrt/rc2data');
            if exist(rc1) && exist(rc2)
                matlabbatch{1}.spm.tools.dartel.warp1.images = {{rc1},{rc2}}';
                matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {self.dartel_templates(1)};
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {self.dartel_templates(2)};
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {self.dartel_templates(3)};
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {self.dartel_templates(4)};
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {self.dartel_templates(5)};
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {self.dartel_templates(6)};
                matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
                matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;
                self.RunSPMJob(matlabbatch);
            else
                fprintf('Need to run segment first...\n')
            end
        end        
        function Re_Coreg(self,run)
            %will realign and coregister. Right now it cannot deal with
            %multiple runs simultaneously.
                                    
            if exist(self.epi_path(run))
                
                mean_epi    = regexprep( self.epi_path(run),'mtr/data','mrt/meandata');
                
                %double-pass realign EPIs and reslice the mean image only.
                matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cellstr(self.epi_path(run));
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
                
                %%coregister to skullstrip (only the affine matrix is modified)
                matlabbatch{2}.spm.spatial.coreg.estimate.ref    = cellstr(self.skullstrip);
                matlabbatch{2}.spm.spatial.coreg.estimate.source = cellstr(mean_epi);
                matlabbatch{2}.spm.spatial.coreg.estimate.other  = cellstr(self.epi_path(run));
                matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                
                %%write EPIs
                matlabbatch{3}.spm.spatial.realign.write.data            = cellstr(self.epi_path(run));
                matlabbatch{3}.spm.spatial.realign.write.roptions.which  = [2 1];%all images as well as the mean image.
                matlabbatch{3}.spm.spatial.realign.write.roptions.interp = 4;
                matlabbatch{3}.spm.spatial.realign.write.roptions.wrap   = [0 0 0];
                matlabbatch{3}.spm.spatial.realign.write.roptions.mask   = 1;
                matlabbatch{3}.spm.spatial.realign.write.roptions.prefix = 'r';
                self.RunSPMJob(matlabbatch(2));
                self.RunSPMJob(matlabbatch(3));
            else
                fprintf('EPI is not here...\n')
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
    methods %fmri path_tools which are related to the subject              
        function out = skullstrip(self)
            out = regexprep(self.hr_path,'mrt/data','mrt/data_skullstrip');
        end
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
            %the directory where hr is located
            out = sprintf('%smrt/data.nii',self.pathfinder(self.id,0));
        end                        
        function [t]        = total_volumes(self,run)
            % will tell you how many volumes are in a 4D image.
            bla = spm_vol_nifti(self.mrt_data(run),1);%simply read the first images header
            t   = bla.private.dat.dim(4);
        end        
        function out        = epi_path(self,nrun)
            % simply returns the path to the mrt data.
            if nargin == 2
                out = sprintf('%smrt/data.nii',self.pathfinder(self.id,nrun));                
            else
                fprintf('Need to give an input...\n')
                return
            end
        end                
        function out        = epi_dir(self,nrun)
            % simply returns the path to the mrt data.
            
            if nargin == 2                
                out = sprintf('%smrt/',self.pathfinder(self.id,nrun));                
            else
                fprintf('Need to give an input...\n')
                return
            end
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
                path2data = regexprep(path2data,'mat',varargin{2});
            end
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
                addpath(self.palamedes_path);
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
end
