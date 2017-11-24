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
        function p      = get_paradigm(self,nrun)
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
			for source = self.dicom_folders(:)'
				%
                n 				 = n+1
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
            o      = length(dir(self.path))-3;%exclude the directories .., ., and run000
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
            L(L(:,1) > last_scan_time,:)  = [];
            L(:,1)          = L(:,1) - first_scan_time;
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
            %vector if needed. RUNS will be used with Re_Coreg.
            if nargin > 1
	    		self.SegmentSurface_HR;%cat12 segmentation
            	self.SkullStrip;%removes non-neural voxels
            	self.MNI2Native;%brings the atlas (if present) to native space
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
        function Segment_meanEPI(self)            
            % Runs the new segment of SPM12 on the mean EPI image.
			% Will write to the disk:
			% meandata_seg8.mat
			% iy_meandata.nii
            % c{1-5}meandata.nii
			% y_meandata.nii
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
            
            for nf = 1:size(path2image,1)
                matlabbatch{nf}.spm.spatial.normalise.write.subj.def      = cellstr(regexprep(self.path_meanepi,'meandata','y_meandata'));
                matlabbatch{nf}.spm.spatial.normalise.write.subj.resample = {path2image(nf,:)};
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70
                                                                              78 76 85];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.vox    = [Inf Inf Inf];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.interp = 4;
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.prefix = 'wEPI_';
            end
            self.RunSPMJob(matlabbatch);
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
    end
    methods %path_tools which are related to the subject              
        function out        = path_skullstrip_meanepi(self)
				out = regexprep(self.path_meanepi,'meandata','ss_meandata');
            end
        function out        = path_meanepi(self)
            %returns the path to the meanepi (result of realignment).
            %returns empty if non-existent. Assumes that the first run contains the mean epi.
            first_run = self.dicom_target_run(1);
	    out       = strrep( self.path_epi(first_run),sprintf('mrt%sdata',filesep),sprintf('mrt%smeandata',filesep));
        end
		function out        = path_tpm(self,n)
			%return the path to the Nth TPM image from the spm
			out = sprintf('%s/TPM.nii,%i',self.tpm_dir,n);
		end
		function out        = path_meanepi_segmented(self,num)
			%Returns the path to the output of Segment_meanEPI, N can be a vector.
			mean_epi = self.path_meanepi;
			out      = '';
			for n = num
				out      = strvcat(out,regexprep(mean_epi,'meandata',sprintf('c%dmeandata',n)))
			end
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
					HRLine = regexp(HRLine,'Series.*','match');
					HRLine = HRLine{1};
                    HRPath = regexp(HRLine,'/common.*','match');
                    HRPath = HRPath{1};
                    HRPath = regexprep(HRPath,'\n','');
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
        function analysis_spm_firstlevel(self,nrun,model_num)            
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
            
            for session = nrun
                %load files using ...,1, ....,2 format
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans  = cellstr(spm_select('expand',self.path_epi(session,'r')));%use always the realigned data.
                %load the onsets
                dummy                                                   = load(self.path_model(session,model_num));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond   = dummy.cond;
                %load nuissance parameters
                nuis                                                    = self.get_param_motion(session);
                nuis                                                    = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
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
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs                  = self.derivatives;%we have [0 0], [ 1 0] or [ 1 1] for 1, 2, or 3 regressors.
            matlabbatch{1}.spm.stats.fmri_spec.volt                              = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                            = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                           = -Inf;%
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = {self.path_skullstrip_meanepi};
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'none';
            spm_jobman('run', matlabbatch);%create SPM file first            
            %% normalize and smooth beta images right away.            
            beta_images          = self.path_beta(nrun(1),model_num,'');%'' => with no prefix            
            %
            %normalize them ('w_EPI' and 'w_CAT' will be added)
            self.VolumeNormalize(beta_images);
            %now smooth these normalized images
            beta_images          = self.path_beta(nrun(1),model_num,'wEPI_');%smooth the normalized images too.
            self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_w_')
            beta_images          = self.path_beta(nrun(1),model_num,'wCAT_');%smooth the normalized images too.
            self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_w_')
            %%                        
        end
     end
end
