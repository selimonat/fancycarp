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
        derivatives       = [0 0];%specifies expansion degree of the cHRF when running models.
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
                try
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn = s.paradigm{s.default_run}.stim.cs_neg;
                end
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
        
            
            fprintf('dump_hr:\nWill now dump the latest HR (%s)\n',self.current_time);            
            %target location for the hr: self.hr_dir;            
            %create it if necess.
            if exist(self.hr_dir) == 0
                mkdir(self.hr_dir);
            end            
            self.DicomDownload(self.GetDicomHRpath,self.hr_dir);
            self.ConvertDicom(self.hr_dir);
            files       = spm_select('FPListRec',self.hr_dir,'^sTRIO');
            if ~isempty(files)
                movefile(files,regexprep(files,'/sTRIO.*$','/data.nii'));%rename it to data.nii
            end
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
                dest             = self.epi_dir(n);                                
                self.DicomDownload(source{1},dest);                
                self.DicomTo4D(dest);
            end
            
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
        function out = GetTotalVolume(self,run)
            %Returns the number of volumes present in a 4D image
            out = size(spm_select('expand',self.epi_path(run)),1);            
        end
        function out = GetTotalVolumeLogged(self,run)
            %returns number of pulses logged in stimulus computer during the experiment
            L   = self.get_log(run);
            out = sum(L(:,2) == 0);
        end
    end
    
    methods %(mri, preprocessing))      
        
        function preprocess_pipeline(self,runs)
            %meta method to run all the required steps for hr
            %preprocessing. RUNS specifies the functional runs, make it a
            %vector if needed.
            self.SegmentSurface;
            self.SkullStrip;            
            self.Re_Coreg(runs);
        end
           
        
        function SkullStrip(self)
            %needs results of SegmentSurface, will produce a skullstripped
            %version of hr (filename: ss_data.nii). It will also
            %automatically create a normalized version as well
            %(w_ss_data.nii).
            %c
            c1         = regexprep(self.hr_path,'mrt/data','mrt/mri/p1data');
            c2         = regexprep(self.hr_path,'mrt/data','mrt/mri/p2data');
            if exist(c1) && exist(c2)
                matlabbatch{1}.spm.util.imcalc.input            = cellstr(strvcat(self.hr_path,c1,c2));
                matlabbatch{1}.spm.util.imcalc.output           = self.skullstrip;
                matlabbatch{1}.spm.util.imcalc.outdir           = {self.hr_dir};
                matlabbatch{1}.spm.util.imcalc.expression       = 'i1.*((i2+i3)>0.2)';
                matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask     = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp   = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype    = 4;
                self.RunSPMJob(matlabbatch);
                
                self.VolumeNormalize(self.skullstrip);
            else
                fprintf('Need to run segment first...\n')
            end
        end        
        
        function Re_Coreg(self,runs)
            %will realign and coregister. 
              
            %% collect all the EPIs as a cell array of cellstr
            c = 0;
            for nr = runs
                if exist(self.epi_path(nr))
                    c = c +1;
                    epi_run{c} = cellstr(spm_select('expand',self.epi_path(nr)));
                end
            end
            %%
                mean_epi    = regexprep( self.epi_path(1),'mrt/data','mrt/meandata');%mean EPI will always be saved here.
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
                matlabbatch{2}.spm.spatial.coreg.estimate.ref    = cellstr(self.skullstrip);
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
        function SegmentSurface(self)            
            %runs CAT12 Segment Surface routine.
            matlabbatch{1}.spm.tools.cat.estwrite.data = {spm_select('expand',self.hr_path)};
            matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
            matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {sprintf('%s/TPM.nii',self.tpm_dir)};
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
        function VolumeNormalize(self,path2image)
            %SegmentSurface writes deformation fields (y_*), which are here used
            %to normalize the native hr images. Adds a prefix w- to
            %resampled images. path2image is the image to be resampled.
            %Example: 
            %
            %s.VolumeNormalize(s.beta_path(1,1))
            %s.VolumeNormalize(s.skullstrip);
            
            for nf = 1:size(path2image,1)
                matlabbatch{nf}.spm.spatial.normalise.write.subj.def      = cellstr(regexprep(self.hr_path,'data.nii','mri/y_data.nii'));
                matlabbatch{nf}.spm.spatial.normalise.write.subj.resample = {path2image(nf,:)};
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70
                                                          78 76 85];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.vox    = [Inf Inf Inf];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.interp = 4;
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.prefix = 'w_';
            end
            %
            self.RunSPMJob(matlabbatch);
        end
        
    end
    methods %fmri path_tools which are related to the subject              
        function out        = skullstrip(self,varargin)
            %returns filename for the skull stripped hr. Use VARARGIN to
            %add a prefix to the output, such as 'w' for example.
            if nargin == 1
                out = sprintf('%s%s',self.hr_dir,'ss_data.nii');
            elseif nargin == 2
                out = sprintf('%s%s_%s',self.hr_dir,varargin{1},'ss_data.nii');
            end
        end        
        function out = spmmat_dir(self,nrun,model_num)
            %Returns the path to SPM folder in a given NRUN responsible for
            %the model MODEL_NUM. VARARGIN is used for the derivatives.            
            out = sprintf('%s/spm/model_%02d_chrf_%d%d/',self.path2data(nrun),model_num,self.derivatives(1),self.derivatives(2));
        end
        
        function out        = spmmat_path(self,nrun,model_num)
            %returns the path to spm folder for run RUN.
            dummy = self.spmmat_dir(nrun,model_num,self.derivatives);
            out   = sprintf('%s/SPM.mat',dummy);
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
            bla = spm_vol_nifti(self.epi_path(run),1);%simply read the first images header
            t   = bla.private.dat.dim(4);
        end        
        function out        = epi_path(self,nrun,varargin)
            % simply returns the path to the mrt data. use VARARGIN to add
            % prefixes.
            if nargin == 2
                out = sprintf('%smrt/data.nii',self.pathfinder(self.id,nrun));                
            elseif nargin == 3
                out = sprintf('%smrt/%sdata.nii',self.pathfinder(self.id,nrun),varargin{1});
            else                
                fprintf('Need to give an input...\n')
                return
            end
        end                        
        function out        = epi_dir(self,nrun)
            % simply returns the path to the mrt data.
            
            if nargin == 2                
                out = sprintf('%smrt/',self.pathfinder(self.id,self.dicom_target_run(nrun)));                
            else
                fprintf('Need to give an input...\n')
                return
            end
        end   
        function out = beta_path(self,nrun,model_num,prefix,varargin)
            %returns the path for beta images computed in NRUN for
            %MODEL_NUM. Use VARARGIN to select a subset by indexing.
            %Actually spm_select is not even necessary here.
           
            out = self.spmmat_dir(nrun,model_num);
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
        
        function [HRPath]   = GetDicomHRpath(self)
            % finds the dicom path to the latest HR measurement for this
            % subject.
             
            HRPath = [];
            if ~ismac & ~ispc
                [status2 DicqOutputFull] = system(sprintf('/common/apps/bin/dicq --verbose  --series --exam=%s --folders',self.trio_session));
                %take the latest anatomical scan.
                [status2 HRLine] = system(sprintf('/common/apps/bin/dicq --verbose  --series --exam=%s --folders | grep mprage | tail -n 1',self.trio_session));
                %                
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
            else
                fprintf('GetDicomHRpath: To use dicom query you have to use one of the institute''s linux boxes\n');
            end
        end
        function path2data  = path2data(self,run,varargin)
            % s.path2data(4) will return the path to the subject's phase 4
            % s.path2data(4,'eye') return the path to the eye data file at the
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
        function plot_motionparams(self,nrun)
            dummy = self.GetMotionParameters(nrun);
            subplot(2,1,1);
            plot(dummy(:,1:3));           
            legend({'x','y' 'z'})
            legend boxoff;ylabel('mm');box off
            axis tight
            ylim([-5 5])
            set(gca,'ygrid','on')
            %
            subplot(2,1,2)
            plot(dummy(:,4:6));
            legend({'pitch' 'roll' 'yaw'});
            legend boxoff;ylabel('degrees');box off;
            xlabel('volumes')
            axis tight
            ylim([-5 5].*10.^-2)
            set(gca,'ygrid','on')
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
        function [stim_scanunit,stim_ids]=StimTime2ScanUnit(self,run)
            %will return stim onsets in units of scan. Will check the Log
            %for stim onsets, it will discard those trials occuring outside
            %the first and last scans based on their time-stamps. In
            %FearAmy, we have the reference measurements, i.e. scanner
            %stops for a while during the experiment. Those onsets are also
            %excluded (if the next scan unit is more than the TR far away).
            
            L               = self.get_log(run);
            scan_times      = L(L(:,2) == 0,1);%find all scan events and get their times            
            scan_id         = 1:length(scan_times);%label pulses with increasing numbers
            last_scan_time  = max(scan_times);%time of the last
            first_scan_time = min(scan_times);
            trial           = 0;
            %collect info on stim onsets and discard those not occurring
            %during scanning.
            stim_times      = L(find(L(:,2)==3),1)';
            valid           = stim_times<last_scan_time & stim_times>first_scan_time;%i.e. during scanning
            %sanity check: check if valid onsets matches to tTRIAL
            if sum(valid) ~= self.paradigm{run}.presentation.tTrial
                fprintf('Number of stimulus onsets found in the log are different than tTRIAL.\n');
                keyboard
            end            
            stim_times = stim_times(valid);            
            %
            trial = 0;
            for stim_time = stim_times;%run stim by stim                
                d                        = scan_times - stim_time;
                first_positive           = find(d > 0,1);%find the first positive value
                decimal                  = d(first_positive)./self.TR;
                if decimal < 1%if stimuli are shown but the scanner is not running
                    trial                = trial + 1;
                    stim_scanunit(trial) = (scan_id(first_positive)-1) + decimal;
                    stim_ids(trial)       = self.paradigm{run}.presentation.dist(trial);
                end
            end            
        end
        %
        function CreateModels(self,run)
            %%%%%%%%%%%%%%%%%%%%%%
            model_num  = 1;
            model_path = sprintf('%sdesign/model%02d/data.mat',self.path2data(run),model_num);
            if ~exist(fileparts(model_path));mkdir(fileparts(model_path));end
            [scan,id]  = self.StimTime2ScanUnit(run);
            counter    = 0;
            for current_condition = unique(id)
                counter                = counter + 1;
                cond(counter).name     = mat2str(current_condition);
                cond(counter).onset    = scan(id == current_condition);
                cond(counter).duration = zeros(1,length(cond(counter).onset));
                cond(counter).tmod     = 0;
                cond(counter).pmod     = struct('name',{},'param',{},'poly',{});
            end
            save(model_path,'cond');
            %%%%%%%%%%%%%%%%%%%%%%
        end
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
            %run the model MODEL_NUM for data in NRUN. NRUN can be a vector.            
                        
            
            %set spm dir: saves always to run1
            spm_dir  = self.spmmat_dir(nrun,model_num);
            spm_path = self.spmmat_path(nrun,model_num);
            if ~exist(self.spm_path);mkdir(spm_dir);end
            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            for session = nrun
                %load files using ...,1, ....,2 format
                out                                                     = self.epi_path(session,'r');
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans  = cellstr(spm_select('expand',out));
                %load the onsets                
                dummy                                                   = load(sprintf('%sdesign/model%02d/data.mat',self.path2data(session),model_num));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = dummy.cond;
                %load nuissance parameters
                nuis                                                    = self.GetMotionParameters(session);
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
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs                  = self.derivatives;%we have [0 0], [ 1 0] or [ 1 1] for 1, 2, or 3 regressors.
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                           = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = {''};%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'none';
            %estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat            = {spm_path};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
            %
            %normalize the beta images right away
            beta_images = self.beta_path(nrun,model_num,'');%'' => with no prefix
            self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)
            self.VolumeSmooth(self,beta_images);%('s_' will be added)
        end
     end
end
