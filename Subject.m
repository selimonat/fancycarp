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
        trio_session  = [];
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
                    s.paradigm{nrun} = s.get_paradigm(nrun);
                end
                try
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn = s.paradigm{s.default_run}.stim.cs_neg;
                end
                s.scr = SCR(s);                       
                
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
            files       = spm_select('FPListRec',self.dir_hr,'^sTRIO');
            if ~isempty(files)
                movefile(files,regexprep(files,sprintf('%ssTRIO.*$',filesep),sprintf('%sdata.nii',filesep)));%rename it to data.nii
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
					keyboard
					fprintf('Stopped here as a sanity check\nIt seems the destination folder doesn''t exist.')
				end
            end
            
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
        function o      = get_param_motion(self,run)
            %will load the realignment parameters, of course you have to
            %realign the EPIs first.

            filename = sprintf('%smrt%srp_data.txt',self.path_data(run),filesep);
            if exist(filename)
                o = load(filename);
            else
                fprintf('File:\n %s doesn''t exist.\n Most likely realignment is not yet done.\n');
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
        function out    = get_totalvolumelogged(self,run)
            %returns number of pulses logged in stimulus computer during the experiment
            L   = self.get_log(run);
            out = sum(L(:,2) == 0);
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
    end
    
    methods %(mri, preprocessing))              
        function preprocess_pipeline(self,runs)
            %meta method to run all the required steps for hr
            %preprocessing. RUNS specifies the functional runs, make it a
            %vector if needed.
            self.SegmentSurface;
            self.SkullStrip;%removes non-neural voxels
%             self.MNI2Native;%brings the atlas to native space
            self.Re_Coreg(runs);            

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
            mean_epi    = strrep( self.path_epi(1),sprintf('mrt%sdata',filesep),sprintf('mrt%smeandata',filesep));
            
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
        function SegmentSurface(self)            
            %runs CAT12 Segment Surface routine.
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
        function VolumeNormalize(self,path2image)
            %SegmentSurface writes deformation fields (y_*), which are here used
            %to normalize the native hr images. Adds a prefix w- to
            %resampled images. path2image is the image to be resampled.
            %Example: 
            %
            %s.VolumeNormalize(s.path_beta(1,1))
            %s.VolumeNormalize(s.path_skullstrip);
            
            for nf = 1:size(path2image,1)

                matlabbatch{nf}.spm.spatial.normalise.write.subj.def      = cellstr(strrep(self.path_hr,'data.nii',sprintf('mri%sy_data.nii',filesep)));
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
        function MNI2Native(self)
            %brings the atlas to MNI space and saves it in run000/atlas.
            %Same as VolumeNormalize but uses the inverse deformation
            %fields but same batch. Currently functions only with the 120th
            %volume (which is right amygdala).
            
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
        end  
        function meanNorm(self)
            VolumeNormalize(self,cellstr(strrep(self.path_hr,sprintf('run000%smrt%sdata.nii',filesep,filesep),sprintf('run001%smrt%smeandata.nii',filesep,filesep))))
        end
%         function meanEPInorm(self)
%         %Create warped T1 and mean EPI
%             matlabbatch = [];
%             strrep(self.path_hr,sprintf('mrt%sdata',filesep),sprintf('mrt%smri%sp1data',filesep,filesep));
%                 st_dir      = [base_dir name '/fMRI' filesep 'mprage' filesep];
%                 u_rc1_file  = spm_select('FPList', st_dir, u_rc1_templ);
%                 strip_file  = spm_select('FPList', st_dir, skullstrip_templ);
% 
%                 f_dir       = [base_dir name '/fMRI' filesep 's1' ];
%                 mean_file   = spm_select('FPList', f_dir, mean_func_templ);
% 
%                 matlabbatch{g}.spm.tools.dartel.crt_warped.flowfields = cellstr(strvcat(u_rc1_file,u_rc1_file));
%                 matlabbatch{g}.spm.tools.dartel.crt_warped.images = {cellstr(strvcat(mean_file,strip_file))};
%                 matlabbatch{g}.spm.tools.dartel.crt_warped.jactransf = 0;
%                 matlabbatch{g}.spm.tools.dartel.crt_warped.K = 6;
%                 matlabbatch{g}.spm.tools.dartel.crt_warped.interp = 1;
%             
%             self.RunSPMJob(matlabbatch);
%         end
    end
    methods %path_tools which are related to the subject              
        function out        = path_skullstrip(self,varargin)
            %returns filename for the skull stripped hr. Use VARARGIN to
            %add a prefix to the output, such as 'w' for example.
            if nargin == 1
                out = sprintf('%s%s',self.dir_hr,'ss_data.nii');
            elseif nargin == 2

                out = sprintf('%s%s_%s',self.dir_hr,varargin{1},'ss_data.nii');
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
                keyboard%sanity check
            end
            %select if VARARGIN provided
            if nargin > 4
                selector        = varargin{1};
                out             = out(selector,:);
            end
        end
        function [HRPath]   = path_hr_dicom(self)
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
            dummy = self.get_param_motion(nrun);
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
            catch_times      = L(find(L(:,2)==11),1)';
            rate_times      = L(find(L(:,2)==7),1)';
            valid_s           = stim_times<last_scan_time & stim_times>first_scan_time;%i.e. during scanning
            valid_c           = catch_times<last_scan_time & catch_times>first_scan_time;%i.e. during scanning
            valid_r           = rate_times<last_scan_time & rate_times>first_scan_time;%i.e. during scanning
            %sanity check: check if valid onsets matches to tTRIAL
            if sum(valid_s) ~= self.paradigm{run}.presentation.tTrial
                fprintf('Number of stimulus onsets found in the log are different than tTRIAL.\n');
                keyboard
            end            
            stim_times = stim_times(valid_s); 
            catch_times = catch_times(valid_c); 
            rate_times = rate_times(valid_r); 
            %
            trial = 0;
            for stim_time = stim_times;%run stim by stim                
                d                        = scan_times - stim_time;
                first_positive           = find(d > 0,1);%find the first positive value
                decimal                  = d(first_positive)./self.TR;
                if decimal < 1%if stimuli are shown but the scanner is not running
                    trial                = trial + 1;
                    stim_scanunit(trial) = (scan_id(first_positive)-1) + decimal;
                    stim_ids(trial)       = self.paradigm{run}.presentation.con_id(trial);
                end
            end  
            for catch_times = catch_times;%run stim by stim                
                d                        = scan_times - catch_times;
                first_positive           = find(d > 0,1);%find the first positive value
                decimal                  = d(first_positive)./self.TR;
                if decimal < 1%if stimuli are shown but the scanner is not running
                    trial                = trial + 1;
                    stim_scanunit(trial) = (scan_id(first_positive)-1) + decimal;
                    stim_ids(trial)       = max(self.paradigm{run}.presentation.con_id)+1;
                end
            end
            for rate_times = rate_times;%run stim by stim                
                d                        = scan_times - rate_times;
                first_positive           = find(d > 0,1);%find the first positive value
                decimal                  = d(first_positive)./self.TR;
                if decimal < 1%if stimuli are shown but the scanner is not running
                    trial                = trial + 1;
                    stim_scanunit(trial) = (scan_id(first_positive)-1) + decimal;
                    stim_ids(trial)       = max(self.paradigm{run}.presentation.con_id)+2;
                end
            end 
        end
        function CreateModels(self,runs)
            %%%%%%%%%%%%%%%%%%%%%
            dummy=load([self.path_raw, sprintf('%02d',self.sub_list(self.id)), '\stim_onset.mat']);
            for run = runs
                model_num  = 1;
                model_path = self.path_model(run,model_num);
                if ~exist(fileparts(model_path));mkdir(fileparts(model_path));end
                for current_condition =1:size(dummy.stim_onset,2)
                    cond(current_condition).name     = mat2str(current_condition);
                    cond(current_condition).onset    = dummy.stim_onset{run,current_condition}';
                    cond(current_condition).duration = self.duration(current_condition)*ones(1,length(cond(current_condition).onset));
                    cond(current_condition).tmod     = 0;
                    cond(current_condition).pmod     = struct('name',{},'param',{},'poly',{});
                end
                current_condition = current_condition+1;
                cond(current_condition).name     = mat2str(current_condition);
                cond(current_condition).onset    = dummy.catch_onset{run}';
                cond(current_condition).duration = self.duration(current_condition)*ones(1,length(cond(current_condition).onset));
                cond(current_condition).tmod     = 0;
                cond(current_condition).pmod     = struct('name',{},'param',{},'poly',{});
                current_condition = current_condition+1;
                cond(current_condition).name     = mat2str(current_condition);
                cond(current_condition).onset    = dummy.rate_onset{run}';
                cond(current_condition).duration = self.duration(current_condition)*ones(1,length(cond(current_condition).onset));
                cond(current_condition).tmod     = 0;
                cond(current_condition).pmod     = struct('name',{},'param',{},'poly',{});
                save(model_path,'cond');
                %%%%%%%%%%%%%%%%%%%%%%
            end
        end        
        function FitHRF(self,nrun,model_num)

            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            
            %set spm dir: saves always to run1
            spm_dir  = self.dir_spmmat(nrun(1),model_num);
            path_spm = self.path_spmmat(nrun(1),model_num);%stuff is always saved to the first run.
            if ~exist(self.path_spm);mkdir(spm_dir);end

            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            dummy                                                   = load(self.path_model(1,model_num));
            cond_use   = [1:length(dummy.cond_name)];
            duration = 2;
            dummyScans = 6;
            for session = nrun
                %load files using ...,1, ....,2 format

                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans  = cellstr(spm_select('expand',self.path_epi(session,'r')));%use always the realigned data.
                %load the onsets
                for conds = 1:length(cond_use);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds).name     = dummy.cond_name{cond_use(conds)};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds).onset = dummyScans+dummy.stim_onset{session,conds};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds).duration = duration;
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds+1).name     = 'rating';
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds+1).onset = dummyScans+dummy.rate_onset{session};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds+1).duration = duration;
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds+2).name     = 'catch';
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds+2).onset = dummyScans+dummy.catch_onset{session};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond(conds+2).duration = 9;
                all_nuis{session}                                       = self.get_param_motion(session);
                n_nuis         = size(all_nuis{session},2);
                for nuis = 1:n_nuis
                      matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nuis) = struct('name', cellstr(num2str(nuis)), 'val', all_nuis{session}(:,nuis));
                end
%                 matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
%                 matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = dummy.cond;
                %load nuissance parameters
%                 nuis                                                    = self.get_param_motion(session);
%                 nuis                                                    = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
%                 for nNuis = 1:size(nuis,2)
%                     matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).val   = nuis(:,nNuis);
%                     matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).name  = mat2str(nNuis);
%                 end
                %
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi               = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi_reg           = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).hpf                 = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                              = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs                  = self.derivatives;%we have [0 0], [ 1 0] or [ 1 1] for 1, 2, or 3 regressors.
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
%             matlabbatch{1}.spm.stats.fmri_spec.mthresh
%             = -Inf;   %it was uncommented in Selim version
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = {''};%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'AR(1)';%'none' was selim setting;
            %estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat            = {path_spm};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
            %
            %normalize the beta images right away
            beta_images = self.path_beta(nrun(1),model_num,'');%'' => with no prefix
            self.VolumeSmooth(beta_images);%('s_' will be added)       
            self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)

        end
        function estContrast(self,nrun,model_num)
            nuis                  = self.get_param_motion(1);
            nuis                  = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
            n_nuis    = size(nuis,2);
            n_sess    = length(nrun);
            path_spm = self.path_spmmat(nrun(1),model_num);
            
            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {path_spm};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

            matlabbatch{2}.spm.stats.con.spmmat = {path_spm};
            matlabbatch{2}.spm.stats.con.delete = 1;

            
            matlabbatch{2}.spm.stats.con.consess{1}.fcon.name   = 'eff_of_int';
            matlabbatch{2}.spm.stats.con.consess{1}.fcon.convec = {[repmat([eye(self.cond_use) zeros(self.cond_use,n_nuis+2)],1,n_sess) zeros(self.cond_use,n_sess)]};

            for co=2:length(self.con_name)+1
                matlabbatch{2}.spm.stats.con.consess{co}.tcon.name    = self.con_name{co-1};
                matlabbatch{2}.spm.stats.con.consess{co}.tcon.convec  = [repmat([self.con_value(co-1,:) zeros(1,n_nuis)],1,n_sess) zeros(1,n_sess)];
                matlabbatch{2}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
            end
            spm_jobman('initcfg');
            spm('defaults', 'FMRI');
            spm_jobman('run',matlabbatch);
            con_images = self.path_con(nrun(1),model_num,'');%'' => with no prefix
            self.VolumeNormalize(con_images);%normalize them ('w_' will be added)
            self.VolumeSmooth(con_images);%('s_' will be added)   
        end
        function analysis_spm_firstlevel(self,nrun,model_num)            
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            
            %set spm dir: saves always to run1
            spm_dir  = self.dir_spmmat(nrun(1),model_num);
            path_spm = self.path_spmmat(nrun(1),model_num);%stuff is always saved to the first run.
            if exist(spm_dir)
                rmdir(spm_dir,'s');
            end
            mkdir(spm_dir);               
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'secs'; %'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            for session = nrun
                %load files using ...,1, ....,2 format
                
                dummy_scan = cellstr(spm_select('expand',self.path_epi(session,'r')));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans  = dummy_scan(7:end); %use always the realigned data.
                %load the onsets
                dummy                                                   = load(self.path_model(session,model_num));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = dummy.cond;
                %load nuissance parameters
                nuis                                                    = self.get_param_motion(session);
                nuis                                                    = nuis(7:end,:);
%                 nuis                                                    = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
                for nNuis = 1:size(nuis,2)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).val   = nuis(:,nNuis);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).name  = mat2str(nNuis);
                end
                %
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi               = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi_reg           = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).hpf                 = self.HParam;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                              = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs                  = [0 0]; % self.derivatives;%we have [0 0], [ 1 0] or [ 1 1] for 1, 2, or 3 regressors.
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
%             matlabbatch{1}.spm.stats.fmri_spec.mthresh
%             = -Inf;  % was not commented %myself
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = {''};%add a proper mask here.;%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'none';%'AR(1)';%'none';%'AR(1)'; % selim was 'none';  %myself
            spm_jobman('run', matlabbatch);%create SPM file first
            % now adapt for session effects.
%             spm_fmri_concatenate(path_spm, [910 895 self.get_total_volumes(nrun)-910-895]);
            
            matlabbatch = [];
            %estimation
            matlabbatch{1}.spm.stats.fmri_est.spmmat            = {path_spm};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
            %                        
%             matlabbatch = [];
%             beta_images = self.path_beta(nrun(1),model_num,'');%'' => with no prefix
%             %normalize the beta images right away            
%             self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)
% %             self.VolumeSmooth(beta_images);%smooth the native images ('s_' will be added, resulting in 's_')
%             beta_images = self.path_beta(nrun(1),model_num,'w_');%smooth the normalized images too.
%             self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_w_')
            %%  
        end
        function image2surface(self,nrun,model_num)
            beta_images = self.path_beta(nrun(1),model_num,'');
            matlabbatch = [];
            for co = 1:size(beta_images,1)
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.data_vol = cellstr(beta_images(co,:));
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.data_mesh_lh = cellstr(strrep(self.path_hr,'data.nii',sprintf('surf%slh.central.data.gii',filesep)));
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.sample = {'maxabs'};
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.interp = {'linear'};
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.datafieldname = '';
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class = 'GM';
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint = 0;
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize = 0.1;
                matlabbatch{co}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint = 1;
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.data_vol = cellstr(beta_images(co,:));
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.data_mesh_lh = cellstr(strrep(self.dartel_templates(1),sprintf('templates_1.50mm%sTemplate_1_IXI555_MNI152.nii',filesep),sprintf('templates_surfaces%slh.central.Template_T1_IXI555_MNI152.gii',filesep)));
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.sample = {'maxabs'};
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.interp = {'linear'};
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.datafieldname = 'intensity';
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.mapping.rel_mapping.class = 'GM';
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.mapping.rel_mapping.startpoint = 0;
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.mapping.rel_mapping.stepsize = 0.1;
%                 matlabbatch{co}.spm.tools.cat.stools.vol2surftemp.mapping.rel_mapping.endpoint = 1;
            end
            spm_jobman('run',matlabbatch);
            mapped_surf  = spm_select('FPList', strrep(self.path_hr,'data.nii',sprintf('surf%s',filesep)), '^[lr]h.intensity.beta_*'); 
            for i=1:size(mapped_surf,1)
                [p,n,e] = fileparts(mapped_surf(i,:));
                new_name = [p filesep strrep(n,'intensity_beta','beta') '.data'];
%                 new_name = strrep(new_name,'.intensity.Template_T1_IXI555_MNI152_w_','.');
                movefile(mapped_surf(i,:),new_name);
            end
        end
%         function MoveFiles(self,nrun,model_num)
%             mapped_surf = spm_select('FPList', self.dir_spmmat(nrun(1),model_num) , '^[lr]h.beta_*');
%             for i=1:size(mapped_surf,1)
%                 [new_name e p] = fileparts(mapped_surf(i,:));
%                 new_name = strrep(new_name,'fancycarp',sprintf('fancycarp%ssurfaces',filesep));
%                 if ~exist(new_name);   mkdir(new_name);  end
%                 movefile(mapped_surf(i,:),[new_name filesep e p]);
%             end
%         end
        function surface_normalize_smooth(self)
            matlabbatch = [];
            mapped_surf  = spm_select('FPList', strrep(self.path_hr,'data.nii',sprintf('surf%s',filesep)) , '^[lr]h.beta_.*');%cave only 10 and up!!
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.data_surf = cellstr(mapped_surf);
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.fwhm = 6;
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.nproc = 0;
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.lazy = 0;
            spm_jobman('run',matlabbatch);
        end
        
    end
end
