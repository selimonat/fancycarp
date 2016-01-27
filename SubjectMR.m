classdef SubjectMR < ProjectMR
    properties (Hidden)
        paradigm
        default_run  = 1;
    end
    properties (SetAccess = private)
        id
        path
        csp
        csn
        scr    
        pmf
		trio_name    =[];
		trio_folder  =[];
        
    end
    methods
        function s = SubjectMR(id)%constructor
            s.id              = id;
            s.path            = s.pathfinder(s.id,[]);
			s.trio_name 	  = s.trio_names{s.id};
            s.trio_folder 	  = s.trio_folders{s.id};
			if exist(s.path)
                for nrun = 1:5
                    s.paradigm{nrun} = s.load_paradigm(nrun);
                end
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn = s.paradigm{s.default_run}.stim.cs_neg;
                s.scr = SCR(s);
                try
                    s.pmf = s.getPMF;
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
        
        function ConvertDicom(self)
            %% dicom conversion. ATTENTION: dicoms will be deleted
            % and the converted files will be merged to a 4d file nifti
            % file. This file will be named data.nii.
            %             
            matlabbatch = [];
            for nrun = 1:self.tRuns
                files                                                    = cellstr(spm_select('FPListRec',self.path2data(self.id,nrun),'^MR'));
                matlabbatch{nrun}.spm.util.import.dicom.data             = files;
                matlabbatch{nrun}.spm.util.import.dicom.root             = 'flat';
                matlabbatch{nrun}.spm.util.import.dicom.outdir           = {fileparts(self.path2data(self.id,nrun,'mrt'))};
                matlabbatch{nrun}.spm.util.import.dicom.protfilter       = '.*';
                matlabbatch{nrun}.spm.util.import.dicom.convopts.format  = 'nii';
                matlabbatch{nrun}.spm.util.import.dicom.convopts.icedims = 0;
            end
            fprintf('Dicom conversion s#%i... (%s)\n',self.id,datestr(now,'hh:mm:ss'));            
            spm_jobman('run', matlabbatch);
            fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
            %% delete the dicom files            
            fprintf('Cleaning s#%i... (%s)\n',self.id,datestr(now,'hh:mm:ss'));
            for nrun = 1:self.tRuns
                delete(sprintf('%smrt%sMR*',self.path2data(self.id,nrun),filesep));
            end            
            fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
            %% merge to 4D
            fprintf('Merging s#%i...(%s)\n',self.id,datestr(now,'hh:mm:ss'));
            matlabbatch = [];
            c           = 0;
            for nrun = 1:self.tRuns
                files = spm_select('FPListRec',self.path2data(self.id,nrun),'^fTRIO');
                if ~isempty(files)
                    c     = c + 1;
                    matlabbatch{c}.spm.util.cat.vols  = cellstr(files);
                    matlabbatch{c}.spm.util.cat.name  = 'data.nii';
                    matlabbatch{c}.spm.util.cat.dtype = 0;
                end
            end           
            spm_jobman('run', matlabbatch);
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
                n 				 = n +1;
                dest             = sprintf('%ssub%03d/run%03d/mrt/',self.path_project,self.id,self.trio2run{self.id}(n))
                if exist(dest) == 0
					fprintf('The folder %s doesn''t exist yet, will create it...\n',dest)
					mkdir(dest)
				end
				fprintf('Calling system''s COPY function to dump the data...\n')
				[a b]            = system(sprintf('cp -vr %s/* %s',paths{f},dest));
            	a = 0;
				if a ~= 0
					fprintf('There was a problem while dumping...\n');
					keyboard
				end
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
        
        function out = getPMF(self)
            load(sprintf('%smidlevel%sweibull%sdata.mat',Project.path_project,filesep,filesep));
            out = data(self.id);
            %lines in the output are 
            %1/ CS+ before
            %2/ CS- before
            %3/ CS+ after
            %4/ CS- after
            out.subject_alpha = mean(out.params1(1:2,1),1);
            out.subject_beta  = mean(out.params1(1:2,2),1);
        end        
        function pmfplot(self)
            plotpath = sprintf('%s%sp05%sfigures%sfearcloud_FitPMFs_RE.fig',self.path,filesep,filesep,filesep);
            if exist(plotpath)
            openfig(plotpath);
            else
                fprintf('no figure found!')
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
