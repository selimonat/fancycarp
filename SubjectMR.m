classdef SubjectMR < ProjectMR
    properties (Hidden)
        paradigm
        default_run   = 1;
        dicom_folder  = [];
    end
    properties (SetAccess = private)
        id
        path
        csp
        csn
        scr    
        pmf
		trio_session  = [];
        trio_id       = [];
    end
    methods
        function s = SubjectMR(id)%constructor
            s.id              = id;
            s.path            = s.pathfinder(s.id,[]);
            try
                s.trio_session 	  = s.trio_sessions{s.id};
                s.dicom_folder 	  = s.dicom_folders{s.id};
                s.trio_id         = s.trio_ids{s.id};
            end
            
			if exist(s.path)
                for nrun = 1:5
                    s.paradigm{nrun} = s.load_paradigm(nrun);
                end
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn = s.paradigm{s.default_run}.stim.cs_neg;
                s.scr = SCR(s);
                try
                    dummy = s.fitPMF;
                    s.pmf = dummy.params;
                end
%                 try
%                     s.bold = BOLD(s);
%                 end
            else
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path)
            end
        end
    end
    
    methods
        function GetLatestHR(self)
            %will download the latest HR for this subject to default hr
            %path.
            
            %target location for the hr
            hr_target = self.hr_dir;
            %create it if necess.
            if exist(hr_target) == 0
                mkdir(hr_target)
            end
            self.DicomDownload(self.GetDicomHRpath,self.hr_dir);
        end
        function [HRPath]=GetDicomHRpath(self)
            % finds the dicom path to the latest HR measurement for this
            % subject.
            [status2 DicqOutputFull] = system(sprintf('/common/apps/bin/dicq --verbose  --series --patient=%s --folders',self.trio_id));
            %take the latest anatomical scan.
            [status2 HRLine] = system(sprintf('/common/apps/bin/dicq --verbose  --series --patient=%s --folders | grep mprage | tail -n 1',self.trio_id));
            %
            HRPath = [];
            if ~isempty(HRLine);
                HRPath = regexp(HRLine,'/common\S*','match');
                HRPath = HRPath{1};
                %HRPath = GetDicomPath(HRLine);
                fprintf('All recorded HR data:\n')
                fprintf(DicqOutputFull);
                fprintf('The latest one recorded HR data:\n')
                fprintf(HRLine);               
            else
                fprintf('There is no HR data found for this subject.\n Here is the output of the Dicq:\n');
                fprintf(DicqOutputFull);
                keyboard
            end
            
        end
        function ConvertDicom(self)
            %% dicom conversion. ATTENTION: dicoms will be deleted
            % and the converted files will be merged to a 4d file nifti
            % file. This file will be named data.nii.
            %             
            matlabbatch = [];
            count = 0;
            for nrun = 0:self.tRuns
                files                                                    = cellstr(spm_select('FPListRec',self.path2data(self.id,nrun),'^MR'));
                if ~isempty(files)%only create a batch if there is ^MR files.
                count = count +1;
                matlabbatch{count}.spm.util.import.dicom.data             = files;
                matlabbatch{count}.spm.util.import.dicom.root             = 'flat';
                matlabbatch{count}.spm.util.import.dicom.outdir           = {fileparts(self.path2data(self.id,nrun,'mrt'))};
                matlabbatch{count}.spm.util.import.dicom.protfilter       = '.*';
                matlabbatch{count}.spm.util.import.dicom.convopts.format  = 'nii';
                matlabbatch{count}.spm.util.import.dicom.convopts.icedims = 0;
                end
            end
            fprintf('Dicom conversion s#%i... (%s)\n',self.id,datestr(now,'hh:mm:ss'));            
            spm_jobman('run', matlabbatch);
            fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
            %% delete the dicom files            
            fprintf('Cleaning s#%i... (%s)\n',self.id,datestr(now,'hh:mm:ss'));
            for nrun = 0:self.tRuns
                delete(sprintf('%smrt%sMR*',self.path2data(self.id,nrun),filesep));
            end            
            fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
            %% merge to 4D
            fprintf('Merging s#%i...(%s)\n',self.id,datestr(now,'hh:mm:ss'));
            matlabbatch = [];
            c           = 0;
            for nrun = 1:self.tRuns
                files = spm_select('FPListRec',self.path2data(self.id,nrun),'^fTRIO')
                if ~isempty(files)
                    c     = c + 1;
                    matlabbatch{c}.spm.util.cat.vols  = cellstr(files);
                    matlabbatch{c}.spm.util.cat.name  = 'data.nii';
                    matlabbatch{c}.spm.util.cat.dtype = 0;
                end
            end           
            if ~isempty(matlabbatch)
                spm_jobman('run', matlabbatch);
            end
            fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
            %% delete the 3d fTRIO files            
            fprintf('Deleting the fTRIO images s#%i... (%s)\n',self.id,datestr(now,'hh:mm:ss'));
            for nrun = 1:self.tRuns
                delete(sprintf('%smrt%sfTRIO_*',self.path2data(self.id,nrun),filesep));
            end            
            fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
        end
        function GetDicom(self)            
			%Will dump all the Dicoms based on Sessions entered in the
			%Project Object. trio_folders are folders in the dicom server, trio2run dictates in which run these folders should be dumped.
			fprintf('Making a dicom query, sometimes this might take longer...\n')
			[status paths]   = system(['/common/apps/bin/dicq -t --series --exam=' self.trio_name]);
			paths            = strsplit(paths,'\n');%split paths            
            [status result]  = system(['/common/apps/bin/dicq --series --exam=' self.trio_name]);
            fprintf('This is what I found for you:\n')
			result
			fprintf('Will now dump the following series:\n');
			fprintf('Series %i ---> \n',self.trio_folder);
			fprintf('to these runs:\n'); 
			fprintf('run%03d <--- \n',self.trio2run{self.id})
			fprintf('\n');
			%% save the desired runs to disk
            n = 0;
            for f = self.trio_folder(:)'
                %
                n 				 = n+1;
                dest             = sprintf('%ssub%03d/run%03d/mrt/',self.path_project,self.id,self.trio2run{self.id}(n));
                source           = paths{f};                
                self.DicomDownload(source,dest);   
            end
			fprintf('Happily finished dumping...\n');
        end
        function Realign(self) 
            % Will realign all runs using spm batch.
            matlabbatch = [];
            for nrun  = 1:self.tRuns
                file = sprintf('%smrt%sdata.nii',self.path2data(self.id,nrun),filesep);
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
            mean_file =sprintf('%smrt%smeandata.nii',self.path2data(self.id,1),filesep);%is mean_file always saved to the same place?            
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
        function out = Log(self,run)
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
            out(i,:)   = [];            
        end
        function LogPlot(self,nrun)
            L = self.Log(nrun);
            tevents = size(L,1);
            figure;
            plot(L(1:tevents,1) - L(1,1),L(1:tevents,2),'o','markersize',10);        
            ylim([-2 8]);
            set(gca,'ytick',[-2:8],'yticklabel',{'Rating On','Text','Pulse','Tracker+','Cross+','Stim+','CrossMov','UCS','Stim-','Key+','Tracker-'});
            grid on
            drawnow;
        end
        function o = MotionParameters(self,run)
            %will load the realignment parameters.
            filename = sprintf('%smrt/rp_data.txt',self.path2data(self.id,run));
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
        
        function out = spm_dir(self,run)
            %returns the path to spm folder for run RUN.
            out = sprintf('%smrt/spm/',self.pathfinder(self.id,1));
        end
        
        function out = spm_path(self,run)
            %returns the path to spm folder for run RUN.
            out = sprintf('%smrt/spm/SPM.mat',self.pathfinder(self.id,1));
        end
        
        function out = hr_dir(self)
            %the directory where hr is located
            out = sprintf('%smrt/',self.pathfinder(self.id,0));            
        end
        
        function out = hr_path(self)
            %path to the hr volume            
            out = spm_select('ExtFPList',self.hr_dir,'^sTRIO.*.nii$');
        end
        
        function [t]=total_volumes(self,run)            
            % will tell you how many volumes are in a 4D image.
            bla = spm_vol_nifti(self.mrt_data(run),1);%simply read the first images header
            t   = bla.private.dat.dim(4);
        end
        
        function out = mrt_data(self,nrun)
            % simply returns the path to the mrt data.
            out = sprintf('%smrt/data.nii',self.pathfinder(self.id,nrun));
        end
        
        function out = mrt_data_expanded(self,nrun)
            %returns list of filenames of a 4D nii file using comma
            %separated convention (needed for First levels)
            
            out = spm_select('ExtFPList',fileparts(self.mrt_data(nrun)),'^rdata.nii');
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
                dummy                                                   = load(sprintf('%sdesign/model%02d.mat',self.path2data(self.id,session),model_num));
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
                dummy                                                   = load(sprintf('%sdesign/model%02d.mat',self.path2data(self.id,session),model_num));
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
        
        function out = getPMFraw(self)
            dummy = load(self.path2data(self.id,2,'stimulation'));
            out = dummy.p.psi;
        end
        
        function plotPMFraw(self)
            figure
            tchain = size(self.pmf.presentation.x,2);
            for sub = 1:tchain
                subplot(2,1,sub)
            x = self.pmf.presentation.uniquex;
            y = nanmean(self.pmf.log.xrounded(:,:,sub),2);
            e = nanstd(self.pmf.log.xrounded(:,:,sub),0,2);
            errorbar(x,y,e,'k.','MarkerSize',20)
            end
        end 
        function [out] = fitPMF(self,varargin)
            
            if exist(self.path2data(self.id,2,'pmf')) && isempty(varargin)
                
                load(self.path2data(self.id,2,'pmf'));
                ('PMF Fit found and loaded successfully, will just plot it...')
                
            else
                if ~isempty(varargin)
                    fprintf('Forcefit detected.. will now fit it again...\n')
                elseif ~exist(self.path2data(self.id,2,'pmf'))
                    fprintf('No Fit found yet, will fit it now...\n')
                end
                
                % define a search grid
                searchGrid.alpha = linspace(0,100,100);    %structure defining grid to
                searchGrid.beta  = 10.^[-1:0.1:1];         %search for initial values
                searchGrid.gamma = linspace(0,0.5,100);
                searchGrid.lambda = linspace(0,0.1,100);
                paramsFree = [1 1 1 1];
                PF         = @PAL_Weibull;
                %prepare some variables
                tchain   = size(self.pmf.log.xrounded,3);
                xlevels  = unique(abs(self.pmf.presentation.uniquex));
                NumPos   = NaN(length(xlevels),tchain);
                OutOfNum = NaN(length(xlevels),tchain);
                %first collapse the two directions (pos/neg differences from
                %csp)
                
                for chain = 1:tchain
                    fprintf('Starting to fit chain %g...',chain)
                    %get responses, and resulting PMF from PAL algorithm
                    data = self.pmf.log.xrounded(:,:,chain);
                    rep  = self.pmf.presentation.rep;
                    cl = 0;
                    for l = xlevels(:)'
                        cl = cl+1;
                        ind = find(abs(self.pmf.presentation.uniquex) == l);
                        collecttrials = data(ind,1:rep(ind(1)));
                        collecttrials = collecttrials(:);
                        NumPos(cl,chain)   = sum(collecttrials);% number of "different" responses
                        OutOfNum(cl,chain) = length(collecttrials);%number of presentations at that level
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
                        NumPos(:,chain), OutOfNum(:,chain), searchGrid, paramsFree, PF);
                    fprintf('%s . \n',output.message )
                    fprintf('\n')
                    
                    
                    out.NumPos(chain,:) = NumPos(:,chain);
                    out.OutOfNum(chain,:) = OutOfNum(:,chain);
                    out.PropCorrectData(chain,:) = NumPos(:,chain)./OutOfNum(:,chain);
                    out.params(chain,:) = paramsValues;
                    out.LL(chain,:) = LL;
                    out.exitflag(chain,:) = exitflag;
                    out.PF = PF;
                    out.xlevels = xlevels;
                end
                save(self.path2data(self.id,2,'pmf'),'out')
            end
            
            figure
            tchain = size(out.params,1);
            %plot the Fit
            for chain = 1:tchain
                subplot(tchain,1,chain)
                StimLevelsFine = [min(out.xlevels):(max(out.xlevels)- ...
                    min(out.xlevels))./1000:max(out.xlevels)];
                Fit = out.PF(out.params(chain,:),StimLevelsFine);
                plot(out.xlevels,out.PropCorrectData(chain,:),'k.','Markersize',40);
                set(gca,'Fontsize',12);
                hold on;
                plot(StimLevelsFine,Fit,'g-','Linewidth',3);
                legend('data point','Fit')
                title(sprintf('estimated alpha = %g, LL = %g',out.params(chain,1),out.LL(chain)));
            end
            
        end
        
        function p         = load_paradigm(self,nrun)
            %HAST TO GO TO THE PROJECT ACTUALLY TOGETHER WITH
            %CONDTION_>COLOR DESCRIPTION
            filename = self.path2data(self.id,nrun,'stimulation');
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
        function degree    = stimulus2degree(self,stim_id)
            %will transform condition indices to distances in degrees from
            %the csp face. stim_id is a cell array. This is a subject
            %method as it depends on the subject specific CSP face.
            
            ind_valid     = find(cellfun(@(x) ~isempty(x),regexp(stim_id,'[0-9]')));
            degree        = stim_id;
            for i = ind_valid(:)'
                degree{i} = mat2str(MinimumAngle( 0 , (stim_id{i}-self.csp)*45 ));
            end
        end
        function color     = condition2color(self,cond_id)
            cond_id/45+4;
        end
        function rating    = GetRating(self,run,align)
            %align optional, default = 1;
            if nargin < 3
                align = 1;
            end
            % s is a subject instance
            rating = [];
            if ~isempty(self.paradigm{run})
                rating.y      = self.paradigm{run}.out.rating';
                if align
                    rating.y  = circshift(rating.y,[1 4-self.csp ]);
                end
                rating.x      = repmat([-135:45:180],size(self.paradigm{run}.out.rating,2),1);
                rating.i      = repmat(run          ,size(self.paradigm{run}.out.rating,2),size(self.paradigm{run}.out.rating,1));
                rating.y_mean = mean(rating.y);
            else
                fprintf('no rating present for this subject and run (%d) \n',run);
            end            
        end
        function out    = GetSubSCR(self,run,cond)
            if nargin < 3
                cond=1:8;
            end
            conddummy=[-135:45:180 500 1000 3000];
            % s is a subject instance
            out = [];
            cutnum = self.scr.findphase(run);
            self.scr.cut(cutnum);
            self.scr.run_ledalab;
            self.scr.plot_tuning_ledalab(cond);
            out.y = self.scr.fear_tuning;
            out.x = conddummy(cond);
            out.ind = cutnum;
        end        
        function [o]=tRuns(self)
            %% returns the total number of runs in a folder
            a      = dir(self.path);
            %include all runs except run000
            o      = cellfun( @(x) ~isempty(x), regexp({a(:).name},'run[0-9][0-9][0-9]')).*cellfun( @(x) isempty(x), regexp({a(:).name},'run000'));
            o      = sum(o);
        end
        
    end
end
