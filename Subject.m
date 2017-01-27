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
        default_model_name= [];
    end
    properties (SetAccess = private)
        id
        path
        pupil
        csp
        csn
        scr
        eye
        pmf_param     = [];
        rating_param  = [];
        scr_param     = [];
        trio_session  = [];
        rating        = [];
        total_run     = [];
        pmf        
        detected_oddballs;
        detected_face;%
    end
    %%
    methods
        function s = Subject(id,varargin)%constructor
            fprintf('Subject Constructor for id:%i is called:\n',id);
            s.id               = id;
            s.path             = s.pathfinder(s.id,[]);
            s.dicom_serie_id   = s.dicom_serie_selector{s.id};
            s.dicom_target_run = s.dicom2run;
            try
                s.trio_session 	  = s.trio_sessions{s.id};
            end
            
            if exist(s.path)
                for nrun = 1:s.total_run
                    s.paradigm{nrun} = s.get_paradigm(nrun);
                end
                try
                s.csp   = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn   = s.paradigm{s.default_run}.stim.cs_neg;                
                end
                if any(cellfun(@(i) strcmp(i,'scr'),varargin))
                    s.scr   = SCR(s);
                end
                if any(cellfun(@(i) strcmp(i,'pupil'),varargin))
                    s.pupil = Pupil(s.id,1);
                end
%                 try
%                     s.eye = Fixmat(id,0);
%                 end
            else
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path);
                fprintf('Your path might also be wrong...\n');
            end
            s.default_model_name = sprintf('chrf_%d_%d',s.derivatives(1),s.derivatives(2));%specifies the default model path name.
        end
    end
    
    methods %(Getters)        
        function              get_hr(self)
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
        function              get_epi(self)
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
            
            % save the desired runs to disk into folders specified in dicom2run
            n = 0;
            for source = self.dicom_folders(:)'
                %
                n 				 = n+1;
                dest             = self.dir_epi(self.dicom2run(n));
                if exist(dest)
					self.DicomDownload(source{1},dest);                
				else
					keyboard
					fprintf('Stopped here as a sanity check\nIt seems the destination folder doesn''t exist.')
				end
            end
            % merge the data into 4D
            for nrun = unique(self.dicom2run(:))'
                self.DicomTo4D(self.dir_epi(nrun));
            end
                
        end        
        function [o]        = get.total_run(self)
            % returns the total number of ALL (epi+others) runs in a folder.
            o      = length(dir(sprintf('%s/run*',self.path)));%
        end        
        function L          = get_log(self,run)
            % Loads the ptb log.s Removes all events before/after the
            % first/last scan and defines zero as the first scan.
            %
            % Fearamy specificity: appends the microblock index of events
            % to end+1th column. This removes some repeated code in the
            % object.
            
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
            L(:,1)                        = L(:,1) - first_scan_time;            
            % add the microblock index of events
            % identify each event with its microblock id.
            tevents        = size(L,1);
            mbi            = nan(tevents,1);
            current_mb     = 1;
            microb_events  = L(:,2) ~= 9;
            for n = 1:tevents
                if microb_events(n)%is 0 when it is microblock event.
                    mbi(n,1) = current_mb;
                else%if ti is micro_block event use its info value to update the current_mb
                    current_mb = L(n,3);
                    mbi(n,1)   = current_mb;
                end
            end
            L(:,end+1)     = mbi;
        end
        function L          = get_physio2log(self)            
            %returns logged events from the physio computer in the same
            %format as the log file. Events are aligned to the first valid
            %scan pulse.
            %
            
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
        function o          = get_param_motion(self,run)
            %will load the realignment parameters, of course you have to
            %realign the EPIs first.

            filename = sprintf('%smrt%srp_data.txt',self.path_data(run),filesep);
            if exist(filename)
                o = load(filename);
            else
                fprintf('File:\n %s doesn''t exist.\n Most likely realignment is not yet done.\n');
            end            
        end                
        function cond       = get_modelonsets(self,nrun,model_num)
            %returns stimulus onsets for NRUN defined by model specified by
            %MODEL_NUM
            cond = [];
            load(self.path_model(nrun,model_num));
        end
        function N          = get_param_nuissance(self,nrun)
            %Computes nuissances parameters from MotionParameters.
            N = self.get_param_motion(nrun);
            N = zscore([N [zeros(1,size(N,2));diff(N)] N.^2 [zeros(1,size(N,2));diff(N)].^2 ]);
        end
        function [t]        = get_total_volumes(self,run)
            % will tell you how many volumes are in a 4D image.
            bla = spm_vol_nifti(self.path_epi(run),1);%simply read the first images header
            t   = bla.private.dat.dim(4);
        end        
        function out        = get_totalvolumelogged(self,run)
            %returns number of pulses logged in stimulus computer during the experiment
            L   = self.get_log(run);
            out = sum(L(:,2) == 0);
        end        
        function XYZmm      = get_XYZmmNative(self,mask_id)
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
        function path2mask  = get_NativeMaskPath(self,mask_id)
            %Creates a binary mask image from the 4D native atlas. Selects
            %3D volumes from the 4D native atlas based on MASK_ID, binarize
            %the images individually, merges them to 3D and saves the final
            %volume in run000/mask. Binarization is based on the atlas2mask_threshold
            %property.                        
            %
            V = [];
            for nv = mask_id(:)'
                mask_handle = spm_vol(self.path_native_atlas(nv));%read the mask
                ROI         = spm_read_vols(mask_handle);
                if isempty(V)
                    V = logical(zeros(size(ROI)));                    
                end
                V           = V | (ROI > self.atlas2mask_threshold);%threshold it            
            end          
            %
            filename             =  sprintf('%.0f_',mask_id);
            folder               =  fileparts(self.path_data(0,'mask'));
            if ~exist(folder);mkdir(folder);end
            path2mask            = sprintf('%s/mask_%s.nii',folder,filename(1:end-1));
            mask_handle          = spm_vol(self.path_native_atlas(1));
            mask_handle.fname    = path2mask;            
            spm_write_vol(mask_handle,V);        
        end        
        function D          = get_data(self,file,mask_id)
            %will read the data specified in FILE 
            %FILE is the absolute path to a 3/4D .nii file.            
            %
            %MASK_ID is used to select voxels in the native space, so this
            %only makes sense when the file is also in the native space.
            %
            %This relies on get_XYZmmNative and get_mm2vox methods.
            %                        
            
            vh      = spm_vol(file);
            if spm_check_orientations(vh)
                XYZmm   = self.get_XYZmmNative(mask_id);
                XYZvox  = self.get_mm2vox(XYZmm,vh(1));%in EPI voxel space.
                D       = spm_get_data(vh,XYZvox);
            else
                fprintf('The data in\n %s\n doesn''t have same orientations...',file);
            end            
        end
        function [K]        = get_highpassfilter(self)
            % Get high-pass filter a la SPM.                         
            run_borders              = [[0 910 910+895]+1;[910 910+895  self.get_total_volumes(1)]];
            %
            K(1:size(run_borders,2)) = struct('HParam', self.HParam, 'row',    [] , 'RT',     self.TR ,'X0',[]);
            c = 0;
            for b = run_borders
                c        = c + 1;
                K(c).row = b(1):b(2);
                K(c)     = spm_filter(K(c));
                K(c).X0  = [ones(length(K(c).row),1)*std(K(c).X0(:)) K(c).X0];
            end
        end
        function C          = get_constant_terms(self)
            %return the constant terms for a GLM.
            total_volume = self.get_total_volumes(1);
            C            = [[ones(910,1);zeros(total_volume-910,1)] [zeros(910,1);ones(895,1);zeros(total_volume-910-895,1)] [zeros(910+895,1);ones(total_volume-895-910,1)]];
        end
        function [X,nbasis] = get_designmatrix(self,name,varargin)
            %%
            %will return the design matrix used by spm for a given model at
            %a given run. If nargin == 2, then get_modelonsets is called,
            %one can also directly feed in a cond structure. This would be
            %useful in the case of a mumfordian analysis. NAME is the type
            %of basis functions. NBASIS is the number of basis functions in
            %the current set. For HRF is 1, for Fourier it depends on the
            %order.
            
            if nargin == 4
                fprintf('Will read model %i from run %i...\n',varargin{2},varargin{1});                
                cond                  = self.get_modelonsets(varargin{1},varargin{2});            
            elseif nargin == 3
                fprintf('Cond structure is directly fed into...\n');
                cond                  = varargin{1};
            end
                
            fMRI_T                = 16;
            fMRI_T0               = 1;
            xBF.T                 = fMRI_T;
            xBF.T0                = fMRI_T0;
            xBF.dt                = self.TR/xBF.T;
            xBF.UNITS             = 'scans';
            xBF.Volterra          = 1;            
            xBF.name              = name;
            nbasis = 1;
                        
            if strcmp(name,'Fourier set (Hanning)')                
                xBF.order                  = 4;
                xBF.length                 = self.TR*20;
                self.default_model_name    = sprintf('fourier_%d_%d',4,20);
                nbasis                     = xBF.order*2+1;
            end            
            xBF                   = spm_get_bf(xBF);
            %
            for i = 1:length(cond);%one regressor for each condition
                Sess.U(i).dt     = xBF.dt;%- time bin (seconds)                
                Sess.U(i).ons    = cond(i).onset;%- onsets    (in SPM.xBF.UNITS)
                Sess.U(i).name   = {sprintf('%02d',i)};%- cell of names for each input or cause                
                %no parametric modulation here
                Sess.U(i).dur    =  repmat(0,length(Sess.U(i).ons),1);%- durations (in SPM.xBF.UNITS)
                Sess.U(i).P.name =  'none';
                Sess.U(i).P.P    =  'none';
                Sess.U(i).P.h    =  0;%- order of polynomial expansion
                Sess.U(i).P.i    =  1;%- sub-indices of u pertaining to P
            end
            %
            k                       = self.get_total_volumes(1);
            SPM.xBF                 = xBF;
            SPM.nscan               = k;
            SPM.Sess                = Sess;
            SPM.Sess.U              = spm_get_ons(SPM,1);            
            %
            % Convolve stimulus functions with basis functions
            [X,Xn,Fc]               = spm_Volterra(SPM.Sess.U,SPM.xBF.bf,SPM.xBF.Volterra);
            % Resample regressors at acquisition times (32 bin offset)
            X                       = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
% % %             %sanity checks, the difference to SPM.xX.X should be 0.            
% % %             load(self.path_spmmat(nrun,model_num));
% % %             d = abs(SPM.xX.X(:,1:size(X,2)) - X);
% % %             fprintf('SanityCheck: Difference to SPM''s own DM is %i\n',sum(d(:)));
        end
        %behavioral data getters
        function rating     = get.rating(self)
            %returns the CS+-aligned rating for all the runs.
            for run = unique(self.dicom2run(:))';%don't count the first run
                if isfield(self.paradigm{run}.out,'rating')
                    if ~isempty(self.paradigm{run});
                        rating(run).y      = self.paradigm{run}.out.rating';
                        rating(run).y      = circshift(rating.y,[1 4-self.csp ]);
                        rating(run).x      = repmat([-135:45:180],size(self.paradigm{run}.out.rating,2),1);
                        rating(run).ids    = repmat(self.id,size(self.paradigm{run}.out.rating,2),size(self.paradigm{run}.out.rating,1));
                        rating(run).y_mean = mean(rating.y);
                        rating(run).y_sem  = std(rating.y)./sqrt(2);                                                
                    else
                        fprintf('No rating present for this subject and run (%d) \n',nr);
                    end
                end
            end
        end
        function [out]      = get.rating_param(self)
            %returns the parameters of the pmf fit (1 x parameter);
            out      = [self.fit_rating.params(1,:)];
        end
        function [out]      = get.scr_param(self)
            %returns the parameters of the pmf fit (1 x parameter);
            out      = [self.fit_scr.params(1,:)];
        end
        function out        = get.scr(self)
            %returns scr after cleaning bad microblocks
            if ~isempty(self.scr)                
                self.scr.run_ledalab;                
                y                                                   = self.scr.ledalab.y;%[time trial]
                y                                                   = mean(y);%remove tiem dimension
                                
                y                                                   = reshape(y,[65,9]);%[time mbi cond]
                y([Project.mbi_oddball Project.mbi_transition],:)   = [];
                y(self.ucs_vector == 1,:)                           = [];
                
                out.y                                               = y(:,1:8);
                out.x                                               = repmat([-135:45:180],size(out.y,1),1);
                out.ids                                             = repmat(self.id,size(out.y,1),size(out.y,2));
                out.y_mean                                          = mean(out.y);
                out.y_sem                                           = std(out.y)./sqrt(size(out.y,1));                
            end
        end
        function out        = get_scr_timecourse(self)
            %returns the time-course of scr trials in out.y [time mbi
            %cond], use get_scr to obtain timeless amplitudes
            
           if ~isempty(self.scr)
               
               self.scr.run_ledalab;
               
               y                                                   = self.scr.ledalab.y;%[time trial]
               y                                                   = reshape(y,[size(y,1) 65,9]);%[time mbi cond]
               y(:,[Project.mbi_oddball Project.mbi_transition],:) = [];
               y(:,self.ucs_vector == 1,:)                         = [];
               
               out.y = y;
               out.x = repmat(self.scr.ledalab.x(:,1),[1 size(data.y,2) size(data.y,3)]);
           end
        end
        function f          = get_facecircle_fixmat(self)
            %out        = get_facemaps(self)
            %
            %   generates a fixmat object after aligning the fixation
            %   points during the facecircle task.            
            partition = 1;
            out       = self.get_facecircle(partition);            
            x         = out.raw(1,:);%x
            y         = out.raw(2,:);%y
            wedge     = out.raw(6,:);%in line with ptb rects 
            rects     = self.paradigm{1}.stim.circle_rect;
            for w = unique(wedge);
                i    = (wedge == w);%here we subtract from each fixation the [0 0] coordinate the stimuli it was directed at...
                x(i) = x(i) - rects(w,1);% + rects(2,1);
                y(i) = y(i) - rects(w,2);
            end
            %% start a new fixmat here
            f             = Fixmat([],[]);            
            t             = length(x);
            f.x           = round(x);
            f.y           = round(y);
            f.rect        = [0 0 500 500];
            f.subject     = repmat(self.id,1,t);
            f.phase       = repmat(3,1,t);
            f.deltacsp    = out.raw(8,:);
            f.realcond    = f.deltacsp;
            f.selection   = ~(f.x < f.rect(2) | f.x > (f.rect(2)+f.rect(4)-1) | f.y < f.rect(1) | f.y > (f.rect(1)+f.rect(3)-1) );                        
            f.trialid     = ones(1,length(x));
            f.start       = out.raw(11,:);
            f.stop        = out.raw(15,:);
            f.weight      = out.raw(16,:);
            f.fix         = out.raw(10,:);
            %% remove fixations outside of the image border
            f.selection = ~(f.x < f.rect(2) | f.x > (f.rect(2)+f.rect(4)-1) | f.y < f.rect(1) | f.y > (f.rect(1)+f.rect(3)-1) );
            f.ApplySelection;
            %% round coordinates to pixels
            f.x = round(f.x);
            f.y = round(f.y);
        end
        function out        = get_facecircle(self,partition)
            %
            % fixation data are partitioned into PARTITION many partitions,
            % this is done based on the start time of fixations. If you
            % want to get the grand-mean average, use PARTITION = 1. This
            % method deals with aligning fixation points from PTB space to
            % a CS+ aligned space. Fixations are both spatially and
            % temporally cleaned. Each fixation is assigned to a wedge 
            % index. Based on the wedge angle the fixation points are
            % rotated for alignment. The .RAW field contains for each
            % single fixation different information, like position on the
            % screen etc, see the list below.
            %                     .RAW field is [ 
            %                            1/X coor;
            %                            2/Y coor
            %                            3/duration of fixation
            %                            4/Angle of the fixation
            %                            5/distance to center
            %                            6/Wedge Index, meaning position that PTB used to draw
            %                            7/Wedge Angle, the angle in PTB definition,
            %                            basically rounded fixation angles.
            %                            8/Face angle, this is the angular distance to CS+
            %                            9/Aligned Face Index computed based on 8, CS+:4
            %                            10/Rank; this is the order of the fixation along the 30 seconds
            %                            11/START time
            %                            12/angle of fixation maps aligned to CS+.
            %                            13/old amp.
            %                            14/subject index.
            %                            15/STOP time
            %                            16/Weight for fixations based on their wedge index.
            %
            %
            % mind that the .RAW field interacts with PARTITION, if you
            % need to have all fixations in the .RAW fields you have to
            % have PARTITION set to 1.
            
            out                               = [];
            filename                          = sprintf('%sfacecircle_%02d.mat',self.path_midlevel(3),partition);%facecircle is recorded in run 3.
            force = 0;
            %if you like to recache the Project.screen_size and
            %Fixmat.window has to be readjusted to original values, cache
            %it and change back to aligned settings.
            if exist(filename) == 0 || force;
                
                % get the masks to count the fixations as a f(cs+ distance)
                res                               = self.screen_resolution;%resolution of the screen needed to get the center
                center                            = res./2;%[y x]
                wedge_angle                       = 135:-45:-180;%this is the midpoint of the wedge. This is in PTB register, i.e. WedgeIndex of 1 means 135 degrees.
                %
                fix                               = Fixmat(self.id,3);%get the eye data from that face
                valid                             = (fix.start >= 0)&(fix.start < 30000);%take only fixations which are within the limits.
                x                                 = fix.x(valid)-center(2);%center X and Y
                y                                 = res(1)-fix.y(valid)-center(1);%image to cartesian coordinates for the y-axis.
                W                                 = double(abs(fix.stop(valid) - fix.start(valid)));%duration
                rank                              = double(fix.fix(valid));%rank of the fixation
                start                             = double(fix.start(valid));%start of the fixation point.
                stop                              = double(fix.stop(valid));%start of the fixation point.
                %extract parameter from x and y on the ptb coordinates
                [theta, rho]                      = cart2pol(x,y);%transform X and Y to
                [theta]                           = rad2deg(theta);%angle of the fixation point in degrees
                shift_angle                       = -mod(theta+22.5-45,360);%shift the angles so that we can get the wedge index.
                wedge_index                       = mod(ceil(shift_angle./45)+3-1,8)+1;%moves together the same as PTB drawing
                face_angle                        = self.paradigm{1}.stim.circle_order(wedge_index)';%drawn condition in the wedge
                face_index                        = face_angle./45+4;
                out.raw                           = [x+center(2);-y+res(1)-center(1);double(W);theta;rho;wedge_index;wedge_angle(wedge_index);face_angle;face_index;rank;start];
                %remove fixations wrt the annulus (spatial cleaning)
                valid                             = (out.raw(5,:) > 190)&(out.raw(5,:) < 380);%take only those fixations which are on an annulus containing face stimuli.
                out.raw(:,~valid)                 = [];%remove again
                %sort fixation based on their drawing order
                [~,i]                             = sort(out.raw(6,:));%sort according to wedge_index i.e. drawing order.
                out.raw                           = out.raw(:,i);
                %get aligned angles
                [out.raw(12,:),out.raw(13,:)]     = pol2cart(deg2rad(out.raw(8,:)+(out.raw(4,:)-out.raw(7,:))),out.raw(5,:));%(new theta, old amp)
                %save the subject's index.
                out.raw(14,:)                     = self.id;
                %add the stop later to conserve compatbility of indices 
                stop(~valid)                      = [];
                out.raw(15,:)                     = stop(i);
                %
                %get the weight business
                W1                                = ones(1,8);%weight for each wedge indices.
                try
                    filename2     = sprintf('%smidlevel/run003/facecircle_position_weights_limit_%i.mat',self.path_project,60);
                    load(filename2);
                    W2           = W;
                catch
                    cprintf([1 0 0],'ouups!!! weights cannot be loaded correctly\n');
                    W2 = W1;
                end
                out.raw(16,:) = W2(out.raw(6,:));%store the weights, this will be used for the fixmat
                %compute two counts based on weights
                tfix          = double(max(fix.fix));
                limits        = linspace(0,tfix,partition+1)+1;
                [a b]         = histc(out.raw(10,:),limits);
                for P   = 1:partition
                    i                         = b == P;
                    out.count(P,:)    = accumarray(out.raw(9,i)',1,[8 1] );%no weight
                    out.countw(P,:)   = accumarray(out.raw(9,i)',W2(out.raw(6,i)'),[8 1] );%posiiton weights
                    out.x(P,:)        = [-135:45:180];
                    out.ids(P,:)      = repmat(P,1,8);
                end
                %add these weights to the raw data matrix also.
                cprintf([0 1 0],'facecircle for partition %03d is now cached...\n',partition);
                save(filename,'out');
            else
                fprintf('facecircle for partition %03d loading from cache...\n',partition);
                load(filename);
            end
        end        
        function [detected] = get.detected_oddballs(self)
            %returns number of oddbals that are detected.
                        
            keypresses  = find(self.paradigm{1}.out.response);
            oddballs    = find(self.paradigm{1}.presentation.oddball);
            detected    = ismember(oddballs,keypresses);
%             detected    = sum(detected);
            
        end
        function selected   = get.detected_face(self)
            %returns the face that is selected as paired with UCS in aligned coordinates.            
            selected   = self.paradigm{1}.out.selectedface;                            
        end
        function out        = get.pmf(self)
            %will load the raw pmf data.
            
            %first a double check for hiwi-fakups
            if 1 || self.csp == (self.paradigm{2}.stim.cs_plus)./45+1               
                out               = self.paradigm{2}.psi;
                %create a third chain by pooling responses of thefirst 2
                %chains
                out.log.xrounded = cat(3,[out.log.xrounded(:,:,1) nan(15,30)],[out.log.xrounded(:,:,2) nan(15,30)],[out.log.xrounded(:,:,1) out.log.xrounded(:,:,2)]);
                %prepare some variables
                tchain            = size(out.log.xrounded,3);
                xlevels           = unique(abs(out.presentation.uniquex));
                out.NumPos        = NaN(length(xlevels),tchain);
                out.OutOfNum      = NaN(length(xlevels),tchain);
                out.sd            = NaN(length(xlevels),tchain);
                out.y_mean        = NaN(length(xlevels),tchain);
                out.x             = xlevels;
                %first collapse the two directions (pos/neg differences from
                %csp)                
                for chain = 1:tchain
                    %get responses, and resulting PMF from PAL algorithm
                    data          = out.log.xrounded(:,:,chain);%responses 1:yes, 0:no
                    rep           = sum(~isnan(data),2);
                    cl            = 0;
                    for l = xlevels(:)'
                        cl                     = cl+1;                        
                        valid_row              = find(abs(out.presentation.uniquex) == l);
                        D                      = data(valid_row,:);
                        valid_index            = find(~isnan(D));
                        out.NumPos(cl,chain)   = sum(D(valid_index));
                        out.OutOfNum(cl,chain) = length(valid_index);                                                
                        out.y_mean(cl,chain)   = out.NumPos(cl,chain)./out.OutOfNum(cl,chain);
                        out.y_var(cl,chain)    = (out.y_mean(cl,chain)*(1-out.y_mean(cl,chain)))./out.OutOfNum(cl,chain);%var of binomial distr. (np(1-p))                        
                    end
                end
            else
                cprintf([1 0 0],'!!!! Hi-wi fakup detected, will discard this subject''s pmf data...\n')
                out = [];                
            end
        end 
        function out        = get.pmf_param(self)
            %returns the parameters of the pmf fit (chain x parameter);            
            
            if ~isempty(self.pmf)
                out      = [self.fit_pmf.params];
            else
                out = [];
                cprintf([1 0 0],'No pmf data for this subject.\n')
            end
        end          
        function out        = get_bold_fourier_spacetime(self,sk)
            
            out      = [];
            filename = sprintf('%sget_bold_fourier_spacetime_%02d.mat',self.path_midlevel(1),sk);
            %
            if exist(filename) == 0
                fprintf('Collecting bold timexcondition maps: SK:%02d\n',sk);
                xBF                     = [];
                TR                      = self.TR;
                xBF.T                   = 16;
                xBF.T0                  = 1;
                xBF.dt                  = TR/xBF.T;
                xBF.UNITS               = 'scans';
                xBF.Volterra            = 1;
                xBF.name                = 'Fourier set (Hanning)';
                xBF.order               = 4;
                xBF.length              = 20*TR;
                fir_xBF                 = spm_get_bf(xBF);
                %             self.default_model_name = [self.default_model_name '_mumfordian'];
                self.default_model_name = 'fourier_4_20_mumfordian';
                beta_files              = self.path_beta(1,0,'w_');
                vol                     = spm_vol(beta_files);                
                vol                     = spm_smoothto16bit(vol,sk);                
                %
                for nroi = 1:length(self.roi)
                    XYZ                    = vol(1).mat\[self.roi(nroi).xyz 1]';
                    data                   = spm_get_data(vol,XYZ);
                    data                   = reshape(data,[65,11 length(vol)/65/11]);
                    data                   = reshape(data,65*11,9);
                    expansion              = fir_xBF.bf*data';
                    expansion              = expansion';
                    %                 expansion              = mean(expansion);
                    data                   = reshape(expansion,65,11,321);                    
                    out(:,:,:,nroi)        = data(self.mbi_valid,1:8,:);
                end
                fprintf('Saving...\n');
                save(filename,'out');
            else
                fprintf('Loading from cache...\n');
                load(filename);
            end
        end                
        function out        = get_bold_spacetime(self,nroi,sk)
            %out        = get_bold_spacetime(self,nroi,sk)
            %
            %returns the BOLD response in time x condition format with
            %smoothing kernel of SK at roi NROI.
            %
            %this part has to be adapted to the ROI object convention
                       
            run      = 1;
            out      = [];
            filename = sprintf('%sget_bold_spacetime_roi_%03d_sk_%02d.mat',self.path_midlevel(run),nroi,sk);
            if exist(filename) == 0                
                fprintf('Collecting bold timexcondition maps: SK:%02d, ROI:%03d, ROIthreshold:%03d\n',sk,nroi,self.atlas2mask_threshold);
                self.default_model_name = 'chrf_0_0_mumfordian';
                beta_files              = self.path_beta(run,0,'w_');
                vol                     = spm_vol(beta_files);
                vol                     = spm_smoothto16bit(vol,sk);
                XYZ_mm                  = self.get_XYZmmNormalized(nroi);
                XYZ_vx                  = self.get_mm2vox(XYZ_mm,vol(1));                                
                data                    = spm_get_data(vol,XYZ_vx);%[volume x voxels]
                data                    = reshape(data,[65,11 size(data,2)]);%distribute to condition and mbi
                out(:,:,:)              = data(self.mbi_valid,1:8,:);%discard bad mbi                
                fprintf('Saving...\n');
                save(filename,'out');
            else
                fprintf('Loading from cache...\n');
                load(filename);
            end
        end        
        function out        = get_pupil_spacetime(self)
                fprintf('Collecting pupil timexcondition maps\n');
                dummy  = self.pupil.get_singletrials(self.id,'clean');
                out    = dummy(self.id).mean(:,1:8);
        end        
        function out        = get_all_spacetime(self,sk)
            out        = [];            
            out        = cat(4,out,self.get_bold_spacetime(sk));                   
            out        = cat(4,out,self.get_pupil_spacetime);
        end        
        
    end
    
    methods %(preprocessing))              
        function preprocess_pipeline(self,runs)
            %meta method to run all the required steps for hr
            %preprocessing. RUNS specifies the functional runs, make it a
            %vector if needed.
            self.SegmentSurface;
            self.SkullStrip;%removes non-neural voxels
            self.MNI2Native;%brings the atlas to native space
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
              
            % collect all the EPIs as a cell array of cellstr
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
            
            matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(strrep(self.path_hr,'data.nii',sprintf('mri%sy_data.nii',filesep)));
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(path2image);
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70
                                                                         78   76  85];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [Inf Inf Inf];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_';
            
            %
            self.RunSPMJob(matlabbatch);
        end                        
        function MNI2Native(self)
            %brings the atlas to native space and saves it in run000/atlas.
            %Same as VolumeNormalize but uses the inverse deformation
            %fields but same batch. Currently functions only with the 120th
            %volume (which is right amygdala).
            
            %copy the atlas to subject's folder.
            copyfile(fileparts(self.path_atlas),[self.path_data(0) 'atlas/']);
            %
            nf = 0;
            for roi = [48 58 120 182]%left, right, bilateral amygdala and visual cortex
                nf                                                          = nf + 1;
                filename                                                    = self.path_native_atlas(roi);%for test purposes and to gain speed I focus on amygdala right now.                
                matlabbatch{nf}.spm.spatial.normalise.write.subj.def        = cellstr(regexprep(self.path_hr,'data.nii','mri/iy_data.nii'));%iy_ not y_!
                matlabbatch{nf}.spm.spatial.normalise.write.subj.resample   = {filename};
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.bb     = [-78 -112 -70
                    78 76 85];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.vox    = [Inf Inf Inf];
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.interp = 4;
                matlabbatch{nf}.spm.spatial.normalise.write.woptions.prefix = 'w';%inverse warped
                %                
            end
            self.RunSPMJob(matlabbatch);
            target_file                                                     = regexprep(self.path_native_atlas,'data.nii','wdata.nii');%created by the above batch;
            movefile(target_file,self.path_native_atlas);
        end
    end
    methods %(fmri analysis)                        
        function beta       = analysis_firstlevel(self,nrun,model_num,mask_id)
            %will compute beta weights "manually" without calling SPM, and
            %will make a sanity check comparing both betas values.
            
            
            X         = self.get_designmatrix(nrun,model_num);%returns the Design Matrix, Nuissance Matrix, and High-pass Filtering Matrix                        
            % Nuissance and high-pass filter Parameters
            N         = self.get_param_nuissance(nrun);
            K         = self.get_highpassfilter;
            %            
            Y         = self.get_data(self.path_epi(nrun,'r'),mask_id);%return the realigned data.
            %this is a little bit spm unorthodox: it will make the beta
            %value to be interpretable in units of z-score.
            Y         = (Y-mean(Y(:)))./std(Y(:));                        
            %
            Y         = spm_filter(K,Y);%high-pass filtering.
            %            
            DM        = [X N ones(size(X,1),1)];%append together Onsets, Nuissances and a constant
            DM        = spm_filter(K,DM);%filter also the design matrix
            DM        = spm_sp('Set',DM);
            DM        = spm_sp('x-',DM);% projector;
            beta      = DM*Y;
            beta      = beta';%(voxels x betas);
%             %sanity check with real betas.
%             path_beta = self.path_beta(nrun,model_num,'',1:8);
%             Y2        = self.get_data(path_beta,mask_id);%return the realigned data.
            
        end                
        function beta       = analysis_mumfordian(self,design_name)
            %This will run a GLM on all single trials one by one while
            %keeping all the other trials in another regressor. DESIGN_NAME
            %goes directly to the xBF.name. 'Fourier set (Hanning)' or
            %'hrf' are your choices.
                                    
            nrun      = 1;            
            % get the nuissance and highpass filtering matrices.
            fprintf('Getting design ingredients...\n')
            N         = self.get_param_nuissance(nrun);            
            K         = self.get_highpassfilter;
            C         = self.get_constant_terms;            
            % get the data and highpass filter right away.
            fprintf('Loading the data and global correction...\n')            
            load(self.path_spmmat(nrun,2));%any SPM model will contain the VY and the global mean normalized values;
            VY        = SPM.xY.VY;%these are all realigned data, but not yet normalized
            clear SPM;       
            %% now we run across all trials and get a specific design matrix for each trial.
            L                            = self.get_log(nrun);
            [onsets,stim_id,mbi]         = self.analysis_StimTime2ScanUnit(L);
            stim_id(stim_id<500)         = stim_id(stim_id<500)./45+4;
            stim_id(stim_id == 1000)     = 9;
            stim_id(stim_id == 1001)     = 10;
            stim_id(stim_id == 1002)     = 11;
%             plot(stim_id,mbi,'o');
            %
            cond                         = [];
            ttrial                       = length(stim_id);
            for ntrial = 1:ttrial
                cond{ntrial}(1).name     = num2str(stim_id(ntrial));
                cond{ntrial}(1).onset    = onsets(ntrial);
                cond{ntrial}(1).duration = 0;
                cond{ntrial}(1).tmod     = [];
                cond{ntrial}(1).pmod     = struct('name',{},'param',{},'poly',{});
                %
                cond{ntrial}(2).name     = 'all_other_stims';%exclude Null trials, UCS and oddball
                cond{ntrial}(2).onset    = setdiff(onsets(stim_id < 9),onsets(ntrial));%everything else
                cond{ntrial}(2).duration = 0;
                cond{ntrial}(2).tmod     = [];
                cond{ntrial}(2).pmod     = struct('name',{},'param',{},'poly',{});
                %
                cond{ntrial}(3).name     = 'ucs_oddballs';
                cond{ntrial}(3).onset    = setdiff(onsets(stim_id > 9),onsets(ntrial));%everything else;%everything else
                cond{ntrial}(3).duration = 0;
                cond{ntrial}(3).tmod     = [];
                cond{ntrial}(3).pmod     = struct('name',{},'param',{},'poly',{});
            end
            %%
            %we will go slice by slice otherwise single subject data is
            %about 12 GB, it stuffes the computer            
            [~,nbasis] = self.get_designmatrix(design_name,cond{1});%run this once to get the number of basis functions
            beta       = nan(VY(1).dim(1),VY(1).dim(2),VY(1).dim(3),nbasis,max(mbi),11,'single');
            tslice     = VY(1).dim(3);
            for nslice = 1:tslice                                
                Y         = squeeze(spm_data_read(VY,'slice',nslice));                
                Y         = reshape(Y,size(Y,1)*size(Y,2),size(Y,3))';
                Ys        = spm_filter(K,Y);%high-pass filtering.
                for ntrial = 1:ttrial
                    fprintf('Fitting Subject %03d''s %03dth trial (stim_id:%02d, onset:%3.5g) of %i, slice %i of slice %i (%2.3g percent)\n',self.id,ntrial,stim_id(ntrial),onsets(ntrial),ttrial,nslice,tslice,ntrial./ttrial*(nslice./tslice)*100);                    
                    [X]   = self.get_designmatrix(design_name,cond{ntrial});%returns the Design Matrix, Nuissance Matrix, and High-pass Filtering Matrix                    
                    %
                    DM        = [X N C];%append together Onsets, Nuissances and a constant
                    DM        = spm_filter(K,DM);%filter also the design matrix
                    DM        = spm_sp('Set',DM);
                    DM        = spm_sp('x-',DM);% projector;
                    dummy     = DM*Ys;%(voxels x betas);
                    beta(:,:,nslice,:,mbi(ntrial),stim_id(ntrial))    = reshape(dummy(1:nbasis,:)',VY(1).dim(1),VY(1).dim(2),nbasis);
                    %(x,y,z,microblock,stim,voxel)
                end
            end
            %%            
            self.default_model_name = [self.default_model_name '_mumfordian'];
            write_folder            = sprintf('%s',fileparts(self.path_spmmat(1,0)));
            write_folder            = regexprep(write_folder,'/projects','/Volumes/VeryBigHardDick');
            if exist(write_folder) == 0
                mkdir(write_folder);
            end
            %%
            
            Vo                      = VY(1);
            Vo.dt                   = [16 0];
            Vo.pinfo                = [1 0 352]';
            c = 0;
            for nbetas = 1:size(beta,4)
                for y = 1:size(beta,6)%stim_id
                    for x = 1:size(beta,5)%mbi
                        c               = c + 1;
                        Vo.fname        = sprintf('%sbeta_%04d.nii',write_folder,c);
                        spm_write_vol(Vo,beta(:,:,:,nbetas,x,y));
                    end
                end
            end
            
        end
        function epi        = analysis_mumfordian_rsa(self)
            
            force    = 0;
            filename = sprintf('%sanalysis_mumfordian_rsa.mat',self.path_midlevel(1));
            if exist(filename) == 0 | force
                %get the beta, these are mumford analysis in native space
                self.default_model_name  ='chrf_0_0_mumfordian';
                betas                    = self.path_beta(1,0,'w_');%(run, model)
                sk_counter               = 0;
                for sk = self.smoothing_factor
                    sk_counter          = sk_counter  +1;
                    epi(sk_counter).sk  = sk;
                    %load and smooth the data.
                    vol                 = spm_vol(betas);
                    vol                 = spm_smoothto16bit(vol,sk);%smooth on the fly, slow but disk efficient.
                    roi_counter         = 0;
                    for nroi = [64 103 120 101 102 126 165 182 163 164]%more or less where the univariate peaks
                        roi_counter                                                          = roi_counter + 1;
                        %load the current roi and detect voxels in steps of 10th percentile
                        d                                                                    = spm_read_vols(spm_vol(self.path_atlas(nroi)));
                        i                                                                    = d>0;
                        borders                                                              = prctile((d(i(:))),linspace(10,100,11));
                        threshold_counter                                                    = 0;
                        for threshold = borders(1:end-1);
                            threshold_counter                                                = threshold_counter +1;
                            fprintf('Subject: %02d, Smooth: %02d, Roi: %03d, Threshold: %03d\n',self.id,sk,nroi,threshold);
                            epi(sk_counter).roi(roi_counter,threshold_counter).name          = self.get_atlasROIname(nroi);
                            epi(sk_counter).roi(roi_counter,threshold_counter).threshold     = threshold;
                            %get the coordinates for ROI and THRESHOLD
                            self.atlas2mask_threshold                                        = threshold;
                            XYZmm                                                            = self.get_XYZmmNormalized(nroi);
                            XYZvox                                                           = self.get_mm2vox(XYZmm,vol(1));%in EPI voxel space.
                            %get the BOLD at these coordinates.
                            D                                                                = spm_get_data(vol,XYZvox);
                            %sanity check: f there a column full of NaN it must
                            %be outside of the FOV, kill it;
                            invalid                                                          = sum(isnan(D)) == size(D,1);
                            D(:,invalid)                                                     = [];
                            %get the similarity metric.
                            D                                                                = reshape(D,[65 11 size(D,2)]);%[mbi condition voxel]
                            D                                                                = D(self.mbi_valid,1:8,:);
                            for nmbi = 1:size(D,1)
                                RRR = corrcoef(squeeze(D(nmbi,:,:))');
%                                 if any(isnan(RRR(:)));
%                                     keyboard
%                                 end
                                epi(sk_counter).roi(roi_counter,threshold_counter).rsa(:,:,nmbi) = RRR;
                            end
                        end
                    end
                    save(filename,'epi');
                end
            else
                load(filename);
            end
        end
        function epi        = analysis_rsa(self)
            
            force    = 0;
            filename = sprintf('%sanalysis_rsa.mat',self.path_midlevel(1));
            if exist(filename) == 0 | force
                %get the beta, these are mumford analysis in native space                
                betas                    = self.path_beta(1,2,'w_',1:8);%(run, model)
                sk_counter               = 0;
                for sk = self.smoothing_factor
                    sk_counter          = sk_counter  +1;
                    epi(sk_counter).sk  = sk;
                    %load and smooth the data.
                    vol                 = spm_vol(betas);
                    vol                 = spm_smoothto16bit(vol,sk);%smooth on the fly, slow but disk efficient.
                    roi_counter         = 0;
                    for nroi = [64 103 120 101 102 126 165 182 163 164]%more or less where the univariate peaks
                        roi_counter                                                          = roi_counter + 1;
                        %load the current roi and detect voxels in steps of 10th percentile
                        d                                                                    = spm_read_vols(spm_vol(self.path_atlas(nroi)));
                        i                                                                    = d>0;
                        borders                                                              = prctile((d(i(:))),linspace(10,100,11));
                        threshold_counter                                                    = 0;
                        for threshold = borders(1:end-1);
                            threshold_counter                                                = threshold_counter +1;
                            fprintf('Subject: %02d, Smooth: %02d, Roi: %03d, Threshold: %03d\n',self.id,sk,nroi,threshold);
                            epi(sk_counter).roi(roi_counter,threshold_counter).name          = self.get_atlasROIname(nroi);
                            epi(sk_counter).roi(roi_counter,threshold_counter).threshold     = threshold;
                            %get the coordinates for ROI and THRESHOLD
                            self.atlas2mask_threshold                                        = threshold;
                            XYZmm                                                            = self.get_XYZmmNormalized(nroi);
                            XYZvox                                                           = self.get_mm2vox(XYZmm,vol(1));%in EPI voxel space.
                            %get the BOLD at these coordinates.
                            D                                                                = spm_get_data(vol,XYZvox);                                                                                    
                            D                                                                = D';
                            
                            RRR                                                              = corrcoef(D);
                            imagesc(RRR);drawnow;title(sprintf('Subject: %02d, Smooth: %02d, Roi: %03d, Threshold: %03d\n',self.id,sk,nroi,threshold));
                            if any(isnan(RRR(:)));
                                keyboard
                            end
                            epi(sk_counter).roi(roi_counter,threshold_counter).rsa = RRR;
                            
                        end
                    end
                    save(filename,'epi');
                end
            else
                load(filename);
            end
        end
    end
    methods %sanity checks
        function sanity_realignement(self)
            %will load both the realigned and original data and compute the
            %correlation between consecutive frames.            
            %%
            V = [];c = 0;
            for nrun = unique(self.dicom2run(:))'
                vol  = spm_select('expand',self.path_epi(nrun));
                volh = spm_vol(vol);
                volh = volh;
                tvol = length(volh);
                %
                for n = 1:tvol
                    c = c+1;
                    fprintf('Reading volumes %3d percent',round((n/tvol)*100));
                    V(:,:,:,c)    = spm_read_vols(volh(n));
                    if n < tvol;fprintf(repmat('\b',1,27));end
                end
                fprintf('\n')
            end
            %
            V        = reshape(V,[size(V,1)*size(V,2)*size(V,3) size(V,4)]);
            R(:,:,1) = CancelDiagonals(corrcoef(V),NaN);
            C(:,1)   = diag(R(:,:,1),-1);
            clear V;
            % load the realigned data 
            V = [];c = 0;
            for nrun = unique(self.dicom2run(:))'
                vol  = spm_select('expand',self.path_epi(nrun,'r'));                
                volh = spm_vol(vol);
                volh = volh;
                tvol = length(volh);
                %
                for n = 1:tvol
                    c = c+1;
                    fprintf('Reading volumes %3d percent',round((n/tvol)*100));
                    V(:,:,:,c)    = spm_read_vols(volh(n));
                    if n < tvol;fprintf(repmat('\b',1,27));end
                end
                fprintf('\n')
            end
            V        = reshape(V,[size(V,1)*size(V,2)*size(V,3) size(V,4)]);
            R(:,:,2) = CancelDiagonals(corrcoef(V),NaN);
            C(:,2)   = diag(R(:,:,2),-1);
            clear V
            %%            
            ffigure(1);imagesc([R(:,:,1) R(:,:,2)]);axis image;colorbar;
            SaveFigure(sprintf('%ssanity_realignment01.png',self.path_midlevel(1)));
            
            ffigure(2);
            plot(C,'o-');            
            SaveFigure(sprintf('%ssanity_realignment02.png',self.path_midlevel(1)));
            %%            
        end
        function sanity_saveslice(self)
            %save a slice from the mean epi resulting from the realignement
            vol  = self.path_meanepi;
            V    = spm_read_vols(spm_vol(vol));
            imagesc(V(:,:,13));
            axis image;
            colormap('gray');
            grid on;
            colorbar;                    
            SaveFigure(sprintf('%ssanity_saveslice.png',self.path_midlevel(1)));
        end
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
        function out        = path_model(self,run,model_num)
            %returns the path to a model specified by MODEL_NUM in run RUN.
            out = sprintf('%sdesign/model%02d/data.mat',self.path_data(run),model_num);
        end
        function out        = path_spmmat(self,nrun,model_num)
            %returns the path to spm folder for run RUN.
            dummy = self.dir_spmmat(nrun(1),model_num);
            out   = sprintf('%s%sSPM.mat',dummy,filesep);
        end        
        function out        = path_native_atlas(self,varargin)
            %path to subjects native atlas, use VARARGIN to slice out a
            %given 3D volume. VARARGIN can be a vector to select more than
            %one volume.
            prefix = sprintf('%satlas/data.nii',self.path_data(0));
            if nargin > 1                
                for nv = 1:length(varargin{1})
                    out{nv} = sprintf('%s,%d',prefix,varargin{1}(nv));
                end
                out = strvcat(out{:});
            else
                out = prefix
            end
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
        function out        = path_beta(self,nrun,model_num,prefix,varargin)
            %returns the path for beta images computed in NRUN for
            %MODEL_NUM. Use VARARGIN to select a subset by indexing.
            %Actually spm_select is not even necessary here.
           
            out1 = self.dir_spmmat(nrun,model_num);
            fprintf('Searching for beta images in:\n%s\n',out1)
            out = spm_select('FPList',out1,sprintf('^%sbeta_*',prefix'));
            if isempty(out)
                cprintf([1 0 0],'No beta images found, probably wrong prefix/run/etc is entered...\n');                
                fprintf('Here:\n%s\n',out1);
                keyboard%sanity check
            end
            %select if VARARGIN provided
            if nargin > 4                
                selector        = varargin{1};
                out             = out(selector,:);
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
                [status2 DicqOutputFull] = system(sprintf('/common/apps/bin/dicq --verbose  --series --exam=%s --folders',self.trio_session));
                %get the patient_id from exam
                a         = regexp(DicqOutputFull,'Patient: V[0-9]*','match');
                PatientID = str2num(a{1}(regexp(a{1},'\d')));
                %take the latest anatomical scan.
                [status2 HRLine] = system(sprintf('/common/apps/bin/dicq --verbose  --series --patient=V%i --folders | grep mprage | tail -n 1',PatientID));
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
        function out        = path_midlevel(self,run)
            out = sprintf('%s%smidlevel%s',self.path_data(run),filesep,filesep);
        end
        function out        = path_meanepi(self)
            %path to meanepi.
           out = strrep( self.path_epi(1),sprintf('mrt%sdata',filesep),sprintf('mrt%smeandata',filesep));
        end
        function out        = dir_epi(self,nrun)
            % simply returns the path to the mrt data.
            
            if ismember(nrun,self.dicom2run)
                fprintf('Requested a run which doesn''t exist\n');
                keyboard%sanity check.
            end
            
            if nargin == 2                
                out = sprintf('%smrt%s',self.pathfinder(self.id,self.dicom_target_run(nrun)),filesep);                
            else
                fprintf('Need to give an input...\n')
                return
            end            
        end
        function out        = dir_spmmat(self,nrun,model_num)
            %Returns the path to SPM folder in a given NRUN responsible for
            %the model MODEL_NUM. VARARGIN is used for the derivatives.            
            out = sprintf('%sspm%smodel_%02d_%s%s',self.path_data(nrun),filesep,model_num,self.default_model_name,filesep);

        end
        function out        = dir_hr(self)
            %the directory where hr is located
            out = sprintf('%smrt%s',self.pathfinder(self.id,0),filesep);
        end
        function out        = path_scr(self)
            %the directory where SCR is located
            out = sprintf('%sscr%sdata.smr',self.pathfinder(self.id,1),filesep);
        end
    end      
    methods %(plotters)
        function plot_log(self,nrun)
            %will plot the events that are logged during the experiment.
            L           = self.get_log(nrun);
            tevents     = size(L,1);            
            scan_events = find(L(:,2) == 0);                        
            stim_events = find(L(:,2) == 3);
            scan_times  = L(scan_events,1);            
            stim_types  = L(stim_events,3);
            stim_times  = L(stim_events,1);
            t_scan      = length(scan_times);            
            %
            plot(L(1:tevents,1),L(1:tevents,2),'o','markersize',10);%plot events as dots.
            % plot lines for pulses
            hold on;                        
            plot([scan_times(5:5:end) scan_times(5:5:end)],ylim,'k','linewidth',.1);%plot every 5 th pulse event as a line
            % text pulse indices for each line as well.            
            text(scan_times(5:5:end),repmat(0,length(5:5:t_scan),1),num2str([5:5:t_scan]'),'color','r');
            % mark with a star missing pulses (if any)
            miss        = find(diff(scan_times) > self.TR*1.1);
            if ~isempty(miss)
                plot(scan_times(miss)+self.TR,0,'mp','markersize',40);
            end
            % text condition ids on dots.            
            text(stim_times,repmat(3,length(stim_times),1),num2str(stim_types),'color','r');            
%             text(stim_times,repmat(3.5,length(stim_times),1),num2str(onsets),'color','r');            
            %
            hold off;            
            set(gca,'ytick',[-2:8],'yticklabel',{'Rating On','Text','Pulse','Tracker+','Cross+','Stim+','CrossMov','UCS','Stim-','Key+','Tracker-'});
            grid off;box off;
            ylim([-2 10]);
            set(gca,'Position',[0.05 .05 1-.1 1]);
            xlim([-100 max(scan_times)+100]);
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
        function plot_activitymap(self,files)
            %plots the data in a volume as average intensity projection, as
            %a 3xN panel, where N is the number of volumes. The masking
            %operates in the native space, so please specify a native
            %image.
            
            %%            
            %%                    
            files      = vertcat(files{:});
            vh         = spm_vol(files);%volume handle
            data_mat   = spm_read_vols(vh);%read the activity data, bonus: it checks also for orientation.
            [XYZmm]    = self.get_XYZmmNative(120);%world space;
            XYZvox     = self.get_mm2vox(XYZmm,vh(1));%get indices for the selected voxels.
            s          = size(data_mat);if length(s) == 3;s = [s 1];end
            mask       = logical(zeros(s(1:3)));%create an empty mask.
            i          = sub2ind(size(data_mat),XYZvox(1,:),XYZvox(2,:),XYZvox(3,:));
            mask(i)    = true;              
            %%
            tcond      = size(data_mat,4);
            %mask the data
            data_mat   = data_mat.*repmat(mask,[1 1 1 tcond]);           
            %
            x          = [min(find(sum(squeeze(sum(mask,3)),2) ~= 0)) max(find(sum(squeeze(sum(mask,3)),2) ~= 0))];
            y          = [min(find(sum(squeeze(sum(mask,3))) ~= 0)) max(find(sum(squeeze(sum(mask,3))) ~= 0))];
            z          = [min(find(sum(squeeze(sum(mask))) ~= 0)) max(find(sum(squeeze(sum(mask))) ~= 0))];
                        
            %vectorize it so that we can get the contribution of all voxels to the
            %colormap computation
            
            %for the first COL conditions
            col                                    = tcond;            
            dummy                                  = reshape(data_mat,prod(s(1:3)),s(4));
            dummy                                  = mean(dummy(mask(:),:),2);
            sigma_cmap = 4;
            [roi.cmap(1:tcond,1) roi.cmap(1:tcond,2)]  = GetColorMapLimits(dummy(:),sigma_cmap);
            %            
            for n = 1:tcond
                current          = data_mat(x(1):x(2),y(1):y(2),z(1):z(2),n);                
                current_mask     = mask(x(1):x(2),y(1):y(2),z(1):z(2));                
                % get the data
                roi.x(:,:,n)     = squeeze(nanmean(current,1));
                roi.y(:,:,n)     = squeeze(nanmean(current,2));
                roi.z(:,:,n)     = squeeze(nanmean(current,3));
                
                % get the alpha masks
                roi.x_alpha      = squeeze(mean(double(current_mask),1) ~=  0);
                roi.y_alpha      = squeeze(mean(double(current_mask),2) ~=  0);
                roi.z_alpha      = squeeze(mean(double(current_mask),3) ~=  0);
                %
            end            
%             [roi.xcmap(1:tcond,1) roi.xcmap(1:tcond,2)]   = GetColorMapLimits(Vectorize(roi.x(:,:,1:tcond)),sigma_cmap);
%             [roi.ycmap(1:tcond,1) roi.ycmap(1:tcond,2)]   = GetColorMapLimits(Vectorize(roi.y(:,:,1:tcond)),sigma_cmap);
%             [roi.zcmap(1:tcond,1) roi.zcmap(1:tcond,2)]   = GetColorMapLimits(Vectorize(roi.z(:,:,1:tcond)),sigma_cmap);
            %%
            fields_alpha = {'x_alpha' 'y_alpha' 'z_alpha' };
%             fields_cmap  = {'xcmap' 'ycmap' 'zcmap' };
            fields_data  = {'x' 'y' 'z'};
            ffigure(1);clf;
            for ncol = 1:tcond;%run over conditions
                for nrow = 1:3%run over rows
                    subplot(3,tcond, ncol+(tcond*(nrow-1)) );                    
                    imagesc(roi.(fields_data{nrow})(:,:,ncol),[roi.cmap(1,1) roi.cmap(1,2)]);
                    alpha(double(roi.(fields_alpha{nrow})));
                    box off  
                    axis off
                    thincolorbar('horizontal');            
                    set(gca,'xtick',[],'ytick',[],'linewidth',1);
                    axis image;
                    drawnow;
                end
            end
            
        end        
        function plot_rating(self)
            
            %plot subjects rating as a bar plot.
            self.plot_bar(mean(self.rating.x), self.rating.y_mean,self.rating.y_std);%plot the data            
            
            %if the fit is better than flat line, paint accordingly.
            if  (self.fit_rating.LL < -log10(.05))%plot simply a blue line if the fit is not significant.
                PlotTransparentLine(self.fit_rating.x_hd,repmat(mean(self.rating.y(:)),100,1),.5,'k','linewidth',2.5);
            else
                PlotTransparentLine(self.fit_rating.x_hd(:),self.fit_rating.y_hd(:),.5,'k','linewidth',2.5);%this is the fit.
            end
            hold off;
            ylim([0 10]);
            title(sprintf('id:%02d (+:%d)',self.id,self.csp),'fontsize',12);%subject and face id
        end
        function plot_pmf(self,chains)
            % plot the fits
            
            if ~isempty(self.pmf)
                colors = {'r' 'c' 'k'};
                for chain = chains(:)'
                    plot(self.fit_pmf.x,self.fit_pmf.y(chain,:),'color',colors{chain},'linewidth',3);
                    hold on;                                
                    errorbar(self.pmf.x,self.pmf.y_mean(:,chain),self.pmf.y_var(:,chain),'o','markersize',8,'color',colors{chain});
                    plot([self.fit_pmf.params(chain,1) self.fit_pmf.params(chain,1)],[0 1],'color',colors{chain});
                end
            end            
            axis tight;box off;axis square;ylim([-0.1 1.2]);xlim([0 135]);drawnow;
            title(sprintf('id:%02d (+:%d)',self.id,self.csp),'fontsize',12);%subject and face id            
            hold on;plot(xlim,[0 0 ],'k-');plot(xlim,[0.5 0.5 ],'k:');plot(xlim,[1 1 ],'k-');hold off;%plot grid lines            
        end
        function plot_scr(self)
            %plot_scr(self);
            if ~isempty(self.scr)
                out = self.scr;
                self.plot_bar(out.x(1,:),out.y_mean,out.y_sem);
                
                %if the fit is better than flat line, paint accordingly.
                if  (self.fit_scr.LL < -log10(.05))%plot simply a blue line if the fit is not significant.
                    PlotTransparentLine(self.fit_scr.x_hd,repmat(mean(self.scr.y(:)),100,1),.5,'k','linewidth',2.5);
                else
                    PlotTransparentLine(self.fit_scr.x_hd(:),self.fit_scr.y_hd(:),.5,'k','linewidth',2.5);%this is the fit.
                end
                
                
                title(sprintf('id:%02d (+:%d)',self.id,self.csp),'fontsize',12);%subject and face id
            end
        end
        function plot_facecircle(self,partition,fun)
            %plot subjects face_circle performance as a bar plot.
            %%
            out  = self.get_facecircle(partition);
            Y = out.countw;
            Y    = self.circconv2(out.countw,[1 1]/2);
            self.plot_bar(Y');%plot the data            
            hold on;
            %if the fit is better than flat line, paint accordingly.
            for P = 1:partition
                if  (self.fit_facecircle(partition).LL(P) < -log10(.05))%plot simply a blue line if the fit is not significant.
                    PlotTransparentLine(linspace(1,8,100)+9*(P-1),repmat(mean(out.countw(P,:)),100,1),.35,'b','linewidth',2.5);
                else
                    PlotTransparentLine(linspace(1,8,100)+9*(P-1),self.fit_facecircle(partition).y(P,:)',.35,'k','linewidth',2.5);%this is the fit.
                end
            end
            hold off;
            axis tight;
            title(sprintf('id:%02d (+:%d)',self.id,self.csp),'fontsize',12);%subject and face id                        
        end
        function plot_facecircle_ptb(self)
            %plots the ptb screen during the face circle "task".
            
            %%            
            %A rect in the Psychophysics Toolbox is a set of rectangular
            %coordinates [x1 y1 x2 y2] specifying the upper left (x1, y1)
            %and lower right (x2, y2) coordinates of a rectangle. The
            %origin (0, 0) of the screen is in the upper left, so y
            %increases from the top to the bottom of the screen, and x
            %increases from left to right.        
            figure(7);clf;
            C = round([self.paradigm{1}.stim.circle_rect(:,1:2) self.paradigm{1}.stim.circle_rect(:,3)-self.paradigm{1}.stim.circle_rect(:,1)   self.paradigm{1}.stim.circle_rect(:,4)-self.paradigm{1}.stim.circle_rect(:,2)]);
            axis ij;
            for n = 1:8;
                rectangle('position',[C(n,:)]);
                text(C(n,1),C(n,2),sprintf('%s\nfile:%d\ndelta:%d',mat2str(n),self.paradigm{1}.stim.circle_file_id(n),self.paradigm{1}.stim.circle_order(n)),'fontsize',20);
            end
            hold on;
            rectangle('position',self.paradigm{1}.ptb.rect)
            hold on;
            axis tight;            
            xlim([0 1024]);
            ylim([0 768]);
            set(gca,'Color','none','xticklabel',[],'yticklabel',[],'xtick',linspace(1,1024,3),'ytick',linspace(1,768,3),'xgrid','on','ygrid','on')            
            circle(1024/2,768/2,190,'r');
            hold on
            circle(1024/2,768/2,380,'r');
            hold off;
            supertitle(mat2str(self.id));
            %it seems the PTB was drawing the face circle from 11:00 oclock
            %clockwise.
        end
    end
    methods %(spm_fmri analysis)
        function [stim_scanunit,stim_ids,stim_mbi] = analysis_StimTime2ScanUnit(self,L)
            %will return stim onsets in units of scan. Will check the Log
            %for stim onsets, it will discard those trials occuring outside
            %the first and last scans based on their time-stamps. This
            %methods accepts now Log data as input, instead of loading the
            %original from disk.
            %
            %In FearAmy, we have the reference measurements, i.e. scanner
            %stops for a while during the experiment. Those onsets are also
            %excluded (if the next scan unit is more than the TR far away).
            %
            
            scan_times          = L(L(:,2) == 0,1);%find all scan events and get their times
            scan_id             = 1:length(scan_times);%label pulses with increasing numbers, assumes no pulses is missing.
            stim_events         = find(L(:,2)==3);
            stim_times          = L(stim_events,1)';
            microblocks         = L(stim_events,end);
            stim_ids_all        = L(stim_events,3);%all stim ids including the to-be discarded ones.
            total_stims         = length(stim_times);
            trial               = 0;
            stim_ids            = [];
            stim_scanunit       = [];
            stim_mbi            = [];
            discarded           = 0;
            for stim_time = stim_times;%run stim by stim
                trial           = trial + 1;
                d               = scan_times - stim_time;
                first_positive  = find(d > 0,1);%find the first positive value
                decimal         = d(first_positive)./self.TR;
                if decimal < 1%if stimuli are shown but the scanner is not running
                    stim_scanunit        = [stim_scanunit; (scan_id(first_positive)-1) + decimal];
                    stim_ids             = [stim_ids     ; stim_ids_all(trial)];
                    stim_mbi             = [stim_mbi     ; microblocks(trial)];
                    fprintf('Stim %03i (%04i) onset converted %5.7g ==> %5.7g\n',trial,stim_ids(end), stim_time, stim_scanunit(end));
                else %probably between the phases.
                    discarded = discarded + 1;
                    fprintf('Found a stim onset happening between phases at %5.7g seconds.\n',stim_time);
                end
            end
            unusual_jumps = sum(diff(scan_times) > 15)-2;
            total_stims2  = discarded+length(stim_ids);
            cprintf([1 0 1],'%i stimulus onsets..\n',length(stim_ids));
            cprintf([1 0 1],'%i scanner  pulses..\n',length(scan_times));
            cprintf([1 0 1],'%i volumes present in the disk..\n',self.get_total_volumes(1));
            cprintf([1 0 1],'%i unusal jumps on the logged scan times..\n',unusual_jumps);
            cprintf([1 0 1],'%i stimuli discarded..\n',discarded);
            cprintf([1 0 1],'%i TOTAL stimuli..\n',discarded+length(stim_ids));
            %
            if unusual_jumps > 0 || total_stims ~= total_stims2
                keyboard
            end
        end
        function                                     analysis_spm_fir(self,nrun,model_num)
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            
            self.default_model_name                                 = 'fir_13_13';
            spm_dir                                                 = self.dir_spmmat(nrun,model_num);
            spm_path                                                = self.path_spmmat(nrun,model_num);
            
            if ~exist(spm_path);mkdir(spm_dir);end
            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            for session = nrun
                %load files using ...,1, ....,2 format
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans     = cellstr(spm_select('expand',self.path_epi(session,'r')));%use always the realigned data.
                %load the onsets
                dummy                                                      = load(self.path_model(session,model_num));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond      = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                %adjust onsets so that the FIR model contains also
                %prestimulus period.
                for nc = 1:length(dummy.cond)
                    dummy.cond(nc).onset = dummy.cond(nc).onset - 5;%push the onset 5s to the past
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond      = dummy.cond;
                %load nuissance parameters
                nuis                                                       = self.get_param_motion(session);
                nuis                                                       = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
                for nNuis = 1:size(nuis,2)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).val   = nuis(:,nNuis);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).name  = mat2str(nNuis);
                end
                %
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi     = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi_reg = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).hpf       = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                        = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length            = self.TR*20;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order             = 20;
            matlabbatch{1}.spm.stats.fmri_spec.volt                        = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                      = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                     = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask                        = {''};%{self.get_NativeMaskPath([48 58])};%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                         = 'none';
            spm_jobman('run', matlabbatch);%create SPM file first
            % now adapt for session effects.
            spm_fmri_concatenate(spm_path, [910 895 self.get_total_volumes(nrun)-910-895]);
            %estimation
            matlabbatch                                                    = [];
            matlabbatch{1}.spm.stats.fmri_est.spmmat                       = {spm_path};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical             = 1;
            spm_jobman('run', matlabbatch);
            %normalize the beta images right away
            beta_images                                                    = self.path_beta(nrun(1),model_num,'');%'' => with no prefix
            self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)
            self.VolumeSmooth(beta_images);%smooth the native images ('s_' will be added, resulting in 's_')
            beta_images                                                    = self.path_beta(nrun(1),model_num,'w_');%smooth the normalized images too.
            self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_ww_')
        end
        function                                     analysis_spm_fourier(self,nrun,model_num,order)
            %run the model MODEL_NUM for data in NRUN.
            %NRUN can be a vector, but then care has to be taken that
            %model_num is correctly set for different runs.
            
            self.default_model_name                                 = 'fourier_20_20';
            spm_dir                                                 = self.dir_spmmat(nrun,model_num);
            spm_path                                                = self.path_spmmat(nrun,model_num);
            
            if ~exist(spm_path);mkdir(spm_dir);end
            
            matlabbatch{1}.spm.stats.fmri_spec.dir                  = {spm_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units         = 'scans';%more robust
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT            = self.TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t        = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0       = 1;
            
            for session = nrun
                %load files using ...,1, ....,2 format
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).scans     = cellstr(spm_select('expand',self.path_epi(session,'r')));%use always the realigned data.
                %load the onsets
                dummy                                                      = load(self.path_model(session,model_num));
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond      = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                %adjust onsets so that the FIR model contains also
                %prestimulus period.
                for nc = 1:length(dummy.cond)
                    dummy.cond(nc).onset = dummy.cond(nc).onset - 5;%push the onset 5s to the past
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).cond      = dummy.cond;
                %load nuissance parameters
                nuis                                                       = self.get_param_motion(session);
                nuis                                                       = zscore([nuis [zeros(1,size(nuis,2));diff(nuis)] nuis.^2 [zeros(1,size(nuis,2));diff(nuis)].^2 ]);
                for nNuis = 1:size(nuis,2)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).val   = nuis(:,nNuis);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session).regress(nNuis).name  = mat2str(nNuis);
                end
                %
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi     = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).multi_reg = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session).hpf       = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact                        = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.length        = self.TR*20;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.order         = 20;
            matlabbatch{1}.spm.stats.fmri_spec.volt                        = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global                      = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                     = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask                        = {''};%{self.get_NativeMaskPath([48 58])};%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                         = 'none';
            spm_jobman('run', matlabbatch);%create SPM file first
            % now adapt for session effects.
            spm_fmri_concatenate(spm_path, [910 895 self.get_total_volumes(nrun)-910-895]);
            %estimation
            matlabbatch                                                    = [];
            matlabbatch{1}.spm.stats.fmri_est.spmmat                       = {spm_path};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical             = 1;
            spm_jobman('run', matlabbatch);
            %normalize the beta images right away
            beta_images                                                    = self.path_beta(nrun(1),model_num,'');%'' => with no prefix
            self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)
            self.VolumeSmooth(beta_images);%smooth the native images ('s_' will be added, resulting in 's_')
            beta_images                                                    = self.path_beta(nrun(1),model_num,'w_');%smooth the normalized images too.
            self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_ww_')
        end
        function                                     analysis_spm_firstlevel(self,nrun,model_num)            
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
            matlabbatch{1}.spm.stats.fmri_spec.mthresh                           = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask                              = {''};%add a proper mask here.;%add a proper mask here.
            matlabbatch{1}.spm.stats.fmri_spec.cvi                               = 'none';
            spm_jobman('run', matlabbatch);%create SPM file first
            % now adapt for session effects.
            spm_fmri_concatenate(path_spm, [910 895 self.get_total_volumes(nrun)-910-895]);
            
            matlabbatch = [];
            %estimation
            matlabbatch{1}.spm.stats.fmri_est.spmmat            = {path_spm};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;
            spm_jobman('run', matlabbatch);
            %                        
            matlabbatch = [];
            beta_images = self.path_beta(nrun(1),model_num,'');%'' => with no prefix
            %normalize the beta images right away            
            self.VolumeNormalize(beta_images);%normalize them ('w_' will be added)
%             self.VolumeSmooth(beta_images);%smooth the native images ('s_' will be added, resulting in 's_')
            beta_images = self.path_beta(nrun(1),model_num,'w_');%smooth the normalized images too.
            self.VolumeSmooth(beta_images);%('s_' will be added, resulting in 's_w_')
            %%                        
        end
        function                                     analysis_spm_fcontrast(self,nrun,model_num,beta_indices)
            %will compute spmf images based on beta images..            
            beta_images = self.path_beta(nrun,model_num,'');
            tbeta       = length(beta_indices);
            con         = 0;
            for n = beta_indices(:)'
                con                                                        = con + 1;%1
                matlabbatch{1}.spm.stats.con.spmmat                        = {self.path_spmmat(nrun,model_num)};
                matlabbatch{1}.spm.stats.con.consess{con}.fcon.name        = sprintf('%03d',n);
                matlabbatch{1}.spm.stats.con.consess{con}.fcon.weights     = {[circshift([1 zeros(1,tbeta-1)],[0 n-1])]}';
                matlabbatch{1}.spm.stats.con.consess{con}.fcon.sessrep     = 'none';
                matlabbatch{1}.spm.stats.con.delete                        = 1;
            end                                    
            spm_jobman('run', matlabbatch);%this will create a lot of spmF_ images.            
            spmf_images = self.path_contrast(nrun,model_num,'','F');
            self.VolumeNormalize(spmf_images);%normalize them ('w_' will be added)
            self.VolumeSmooth(spmf_images);%smooth the native images ('s_' will be added, resulting in 's_')
            spmf_images = self.path_contrast(nrun,model_num,'w_','F');
            self.VolumeSmooth(spmf_images);%('s_' will be added, resulting in 's_w_')
        end
        function                                     analysis_spm_tcontrast(self,nrun,model_num,beta_indices)
            %will compute spmf images based on beta images..            
            load(self.path_spmmat(nrun,model_num));
            tbeta       = length(beta_indices);
            con         = 0;
            matlabbatch = [];
            for n = beta_indices(:)'
                con                                                        = con + 1;%1
                matlabbatch{1}.spm.stats.con.spmmat                        = {self.path_spmmat(nrun,model_num)};
                matlabbatch{1}.spm.stats.con.consess{con}.tcon.name        = sprintf('%03d',n);
                matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights     = [circshift([1 zeros(1,tbeta-1)],[0 n-1])]';
                matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep     = 'none';
                matlabbatch{1}.spm.stats.con.delete                        = 1;
            end
            spm_jobman('run', matlabbatch);%this will create a lot of spmF_ images.            
            spmt_images = self.path_contrast(nrun,model_num,'','T');
            self.VolumeNormalize(spmt_images);%normalize them ('w_' will be added)
            self.VolumeSmooth(spmt_images);%smooth the native images ('s_' will be added, resulting in 's_')
            spmt_images = self.path_contrast(nrun,model_num,'w_','T');
            self.VolumeSmooth(spmt_images);%('s_' will be added, resulting in 's_w_')
        end
    end
    methods %(clearly fearamy specific)
        function analysis_CreateModel01(self)
            %create a model that separates conditions. It pools odd and ucs
            %condition into the same label. discards transition
            %microblocks. This will be used for multicondition mumfordian
            %analysis.
            model_num  = 1;
            for run = 1
                model_path               = self.path_model(run,model_num);%path to the model
                if ~exist(fileparts(model_path));
                    mkdir(fileparts(model_path));
                end
                L                 = self.get_log(run);
                stim_onsets       = L(:,2) == 3;
                mbi               = L(:,end);                
                i                 = ismember(mbi,[22 23 45 46])&stim_onsets;%transition mbi are 22, 23, 45, 46
                L(i,3)            = L(i,3)+1500; 
                L(L(:,3) == 1000,3) = 200

                L(L(:,3) == 1001,3) = 201;
                L(L(:,3) == 1002,3) = 201;
                
                %
                [scan,stim_id,stim_mbi]           = self.analysis_StimTime2ScanUnit(L);
                % sanity plot
                figure(1);i=L(:,2)==3;plot(L(i,3),L(i,end),'k.','markersize',20);title(mat2str(self.id));hold on;
                figure(1);;plot(stim_id,stim_mbi,'ro','markersize',10);title(mat2str(self.id));hold off
                %%
                counter                  = 0;
                for current_condition = unique(stim_id(stim_id(:) <= 500))'
                    counter                = counter + 1;
                    cond(counter).name     = mat2str(current_condition);
                    cond(counter).onset    = scan(stim_id == current_condition);
                    cond(counter).duration = zeros(1,length(cond(counter).onset));
                    cond(counter).tmod     = 0;
                    cond(counter).pmod     = struct('name',{},'param',{},'poly',{});
                end
                save(model_path,'cond');
            end
        end
        function analysis_CreateModel02(self)
            %will generate stim onsets ignoring all microblocks where UCS
            %is delivered as well as transition microblocks..
            model_num      = 2;
            L              = self.get_log(1);%get the log file, will discards all events before/after first/last scan.
            stim_onsets    = L(:,2) == 3;%all stim events.
            scan_onsets    = L((L(:,2) == 0),1);%all scan events.
            ucs_events     = L(:,2) == 5;%all ucs events
            odd_events     = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi            = L(:,end);%microblock identity for all events, this is inserted in the get_log
            
            %now we have to discard mbi where stimuli occured during a
            %transition. This can be achieved by search mbi indices that
            %occured between scans 910-911 and 1805-1806.
            % micro blocks during transition:
            transition_events = (L(:,1)>scan_onsets(910) & L(:,1)<scan_onsets(911))|(L(:,1)>scan_onsets(1805) & L(:,1)<scan_onsets(1806));
            %
            %find all indices of microblocks that are either a transition
            %event, oddball or ucs microblock
            i               = ismember(mbi,mbi(ucs_events|odd_events|transition_events))&stim_onsets;
            %%
            %add 500 degrees to "bad stimulus events".
            L(i,3)          = L(i,3)+500;
            %some more finetuning here basically for future convenience on labelling of conditions.
            i               = L(:,3) ==1500;%shifted null trials: put them back to 1000, we want to model null trials together irrespective of their microblock belongance.
            L(i,3)          = 1000;
            i               = L(:,3) ==1501;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)          = 500;
            i               = L(:,3) ==1502;%shifted oddball trials: put them back to 1002.
            L(i,3)          = 1002;
            i               = L(:,3) ==1000;%shifted oddball trials: put them back to 1002.
            L(i,3)          = 200;            
            %from this point on it is the same as
            %self.analysis_CreateModels01
            model_path               = self.path_model(1,model_num);
            [scan,stim_id,stim_mbi]  = self.analysis_StimTime2ScanUnit(L);%will discard all transition stimuli
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);plot(stim_id,stim_mbi,'ro','markersize',10);title(mat2str(self.id));hold off
            counter         = 0;
            for current_condition = unique(stim_id(:)')
                counter                = counter + 1;
                cond(counter).name     = mat2str(current_condition);
                cond(counter).onset    = scan(stim_id == current_condition);
                cond(counter).duration = zeros(1,length(cond(counter).onset));
                cond(counter).tmod     = 0;
                cond(counter).pmod     = struct('name',{},'param',{},'poly',{});
            end
            if ~exist(fileparts(model_path));
                mkdir(fileparts(model_path));
            end            
            save(model_path,'cond');
        end
        function analysis_CreateModel03(self)
            %same labelling as in MODEL02, but all stimuli from the clean
            %microblocks are added to a single regressor with pmods.
            %Corrupted microblock are treated the same but separately.
            %Nulls and odds are there too but not pmoded.
            %
            
            %this initial part is exactly the same as model_02
            model_num      = 3;            
            L              = self.get_log(1);%get the log file.
            stim_onsets    = L(:,2) == 3;%all stim events.
            scan_onsets    = L((L(:,2) == 0),1);%all scan events.
            ucs_events     = L(:,2) == 5;%all ucs events
            odd_events     = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi            = L(:,end);%microblock identity for all events.
            %now we have to discard mbi where stimuli occured during a
            %transition. This can be achieved by search mbi indices that
            %occured between scans 910-911 and 1805-1806
            % micro blocks during transition:
            transition_events = (L(:,1)>scan_onsets(910) & L(:,1)<scan_onsets(911))|(L(:,1)>scan_onsets(1805) & L(:,1)<scan_onsets(1806));
            %
            % find stim events which are appearing in a microblock where
            % there is an UCS (number 5);
            i                 = ismember(mbi,mbi(ucs_events|odd_events|transition_events))&stim_onsets;
            %
            %add 500 degrees to "bad stimulus events".
            L(i,3)            = L(i,3)+360;
            %some more finetuning here basically for future convenience on labelling of conditions.
            i                 = L(:,3) ==1360;%shifted null trials: put them back to 1000, we want to model null trials together irrespective of their microblock belongance.
            L(i,3)            = 1000;
            i                 = L(:,3) ==1361;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)            = 360;
            i                 = L(:,3) ==1362;%shifted oddball trials: put them back to 1002.
            L(i,3)            = 1002;
            i                 = L(:,3) ==1000;%shifted oddball trials: put them back to 1002.                        
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);            
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);;plot(stim_id,mbi_id,'ro','markersize',10);title(mat2str(self.id));hold off
            %%
            kappa   = .5;
            %create a weight vector for the derivative.
            res     = 8;
            x2      = [0:(res-1)]*(360/res)-135;
            x2      = [x2 - (360/res/2) x2(end)+(360/res/2)];            
            pmod    = NaN(length(stim_id),4);
            for ntrial = 1:length(stim_id)
                if stim_id(ntrial) < 1000
                    pmod(ntrial,1) = 1;%constant term
                    pmod(ntrial,2) = mbi_id(ntrial);%time
                    pmod(ntrial,3) = Tuning.VonMises( stim_id(ntrial),1,kappa,0,0);%amp                    
                end
            end
            %
            pmod(:,2:3)      = nandemean(pmod(:,2:3));                        
            pmod(:,4)        = pmod(:,2).*pmod(:,3);%time x amp                        
            pmod(:,2:end)    = nanzscore(pmod(:,2:end));
            %%                        
            cond             = [];
            %all valid trials
            i                = stim_id < 200;
            cond(1).name     = 'onsets';
            cond(1).onset    = scan(i);
            cond(1).duration = zeros(1,sum(i));
            cond(1).tmod     = 0;
            cond(1).pmod     = struct('name',{'time' 'amp' 'ampxtime' },'param',{ pmod(i,2) pmod(i,3) pmod(i,4) },'poly',{1 1 1 });
%             bar([cond(1).onset]',[cond(1).pmod(4).param]');                        
            %all the rest invalid trials            
            i                = stim_id > 200 & stim_id < 800;
            cond(2).name     = 'onsets_invalid';
            cond(2).onset    = scan(i);
            cond(2).duration = zeros(1,sum(i));
            cond(2).tmod     = 0;
            cond(2).pmod     = struct('name',{'time' 'amp' 'ampxtime' },'param',{pmod(i,2) pmod(i,3) pmod(i,4) },'poly',{1 1 1 });
            %null trials            
            i                = stim_id == 1000;
            cond(3).name     = 'null';
            cond(3).onset    = scan(i);
            cond(3).duration = zeros(1,sum(i));
            cond(3).tmod     = 0;
            cond(3).pmod     = struct('name',{},'param',{},'poly',{});            
            %add also the oddtrials
            i                = stim_id == 1002;
            cond(4).name     = 'odd';
            cond(4).onset    = scan(i);
            cond(4).duration = zeros(1,sum(i));
            cond(4).tmod     = 0;
            cond(4).pmod     = struct('name',{},'param',{},'poly',{});
            
            model_path       = self.path_model(1,model_num);
            model_dir        = fileparts(model_path);
            if ~exist(model_dir);mkdir(model_dir);end
            save(model_path,'cond');
        end        
        function analysis_CreateModel04(self)
            %same as 03, but with more pmod, namely damp/dkappa and its
            %time interactions.
            %
            
            %this initial part is exactly the same as model_02
            model_num      = 4;
            L              = self.get_log(1);%get the log file.
            stim_onsets    = L(:,2) == 3;%all stim events.
            scan_onsets    = L((L(:,2) == 0),1);%all scan events.
            ucs_events     = L(:,2) == 5;%all ucs events
            odd_events     = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi            = L(:,end);%microblock identity for all events.
            %now we have to discard mbi where stimuli occured during a
            %transition. This can be achieved by search mbi indices that
            %occured between scans 910-911 and 1805-1806
            % micro blocks during transition:
            transition_events = (L(:,1)>scan_onsets(910) & L(:,1)<scan_onsets(911))|(L(:,1)>scan_onsets(1805) & L(:,1)<scan_onsets(1806));
            %
            % find stim events which are appearing in a microblock where
            % there is an UCS (number 5);
            i                 = ismember(mbi,mbi(ucs_events|odd_events|transition_events))&stim_onsets;
            %
            %add 500 degrees to "bad stimulus events".
            L(i,3)            = L(i,3)+360;
            %some more finetuning here basically for future convenience on labelling of conditions.
            i                 = L(:,3) ==1360;%shifted null trials: put them back to 1000, we want to model null trials together irrespective of their microblock belongance.
            L(i,3)            = 1000;
            i                 = L(:,3) ==1361;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)            = 360;
            i                 = L(:,3) ==1362;%shifted oddball trials: put them back to 1002.
            L(i,3)            = 1002;
            i                 = L(:,3) ==1000;%shifted oddball trials: put them back to 1002.                        
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);            
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);;plot(stim_id,mbi_id,'ro','markersize',10);title(mat2str(self.id));hold off
            %%
            kappa   = .1;
            %create a weight vector for the derivative.
            res     = 8;
            x2      = [0:(res-1)]*(360/res)-135;
            x2      = [x2 - (360/res/2) x2(end)+(360/res/2)];
            deriv   = -abs(diff(Tuning.VonMises(x2,1,kappa,0,0)));
            pmod    = NaN(length(stim_id),6);
            for ntrial = 1:length(stim_id)
                if stim_id(ntrial) < 1000
                    pmod(ntrial,1) = 1;%constant term
                    pmod(ntrial,2) = mbi_id(ntrial);%time
                    pmod(ntrial,3) = Tuning.VonMises( stim_id(ntrial),1,kappa,0,0);%amp
                    pmod(ntrial,5) = deriv(mod(stim_id(ntrial)./45+4-1,8)+1);%sigma
                end
            end
            pmod(:,2:3)      = nandemean(pmod(:,2:3));%time and amp            
            pmod(:,5)        = nandemean(pmod(:,5));%dsigma mean corrected
            pmod(:,4)        = pmod(:,2).*pmod(:,3);%time x amp            
            pmod(:,6)        = pmod(:,2).*pmod(:,5);%time x dsigma
            pmod(:,2:end)    = nanzscore(pmod(:,2:end));
            %%                        
            cond             = [];
            %all valid trials
            i                = stim_id < 200;
            cond(1).name     = 'onsets';
            cond(1).onset    = scan(i);
            cond(1).duration = zeros(1,sum(i));
            cond(1).tmod     = 0;
            cond(1).pmod     = struct('name',{'time' 'amp' 'ampxtime' 'dyds' 'dydsxtime'},'param',{ pmod(i,2) pmod(i,3) pmod(i,4) pmod(i,5) pmod(i,6)},'poly',{1 1 1 1 1});
%             bar([cond(1).onset]',[cond(1).pmod(4).param]');                        
            %all the rest invalid trials            
            i                = stim_id > 200 & stim_id < 800;
            cond(2).name     = 'onsets_invalid';
            cond(2).onset    = scan(i);
            cond(2).duration = zeros(1,sum(i));
            cond(2).tmod     = 0;
            cond(2).pmod     = struct('name',{'time' 'amp' 'ampxtime' 'dyds' 'dydsxtime'},'param',{pmod(i,2) pmod(i,3) pmod(i,4) pmod(i,5) pmod(i,6)},'poly',{1 1 1 1 1});
            %null trials            
            i                = stim_id == 1000;
            cond(3).name     = 'null';
            cond(3).onset    = scan(i);
            cond(3).duration = zeros(1,sum(i));
            cond(3).tmod     = 0;
            cond(3).pmod     = struct('name',{},'param',{},'poly',{});            
            %add also the oddtrials
            i                = stim_id == 1002;
            cond(4).name     = 'odd';
            cond(4).onset    = scan(i);
            cond(4).duration = zeros(1,sum(i));
            cond(4).tmod     = 0;
            cond(4).pmod     = struct('name',{},'param',{},'poly',{});
            
            model_path       = self.path_model(1,model_num);
            model_dir        = fileparts(model_path);
            if ~exist(model_dir);mkdir(model_dir);end
%             save(model_path,'cond');
        end        
        function analysis_CreateModel05(self)
            %same as 04, but with quadratic terms on time pmods.
            %
            
            %this initial part is exactly the same as model_02
            model_num      = 5;
            L              = self.get_log(1);%get the log file.
            stim_onsets    = L(:,2) == 3;%all stim events.
            scan_onsets    = L((L(:,2) == 0),1);%all scan events.
            ucs_events     = L(:,2) == 5;%all ucs events
            odd_events     = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi            = L(:,end);%microblock identity for all events.
            %now we have to discard mbi where stimuli occured during a
            %transition. This can be achieved by search mbi indices that
            %occured between scans 910-911 and 1805-1806
            % micro blocks during transition:
            transition_events = (L(:,1)>scan_onsets(910) & L(:,1)<scan_onsets(911))|(L(:,1)>scan_onsets(1805) & L(:,1)<scan_onsets(1806));
            %
            % find stim events which are appearing in a microblock where
            % there is an UCS (number 5);
            i                 = ismember(mbi,mbi(ucs_events|odd_events|transition_events))&stim_onsets;
            %
            %add 500 degrees to "bad stimulus events".
            L(i,3)            = L(i,3)+360;
            %some more finetuning here basically for future convenience on labelling of conditions.
            i                 = L(:,3) ==1360;%shifted null trials: put them back to 1000, we want to model null trials together irrespective of their microblock belongance.
            L(i,3)            = 1000;
            i                 = L(:,3) ==1361;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)            = 360;
            i                 = L(:,3) ==1362;%shifted oddball trials: put them back to 1002.
            L(i,3)            = 1002;
            i                 = L(:,3) ==1000;%shifted oddball trials: put them back to 1002.                        
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);                        
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);;plot(stim_id,mbi_id,'ro','markersize',10);title(mat2str(self.id));hold off
            %%
            kappa   = .1;
            %create a weight vector for the derivative.
            res     = 8;
            x2      = [0:(res-1)]*(360/res)-135;
            x2      = [x2 - (360/res/2) x2(end)+(360/res/2)];
            deriv   = -abs(diff(Tuning.VonMises(x2,1,kappa,0,0)));
            pmod    = NaN(length(stim_id),6);
            for ntrial = 1:length(stim_id)
                if stim_id(ntrial) < 1000
                    pmod(ntrial,1) = 1;%constant term
                    pmod(ntrial,2) = mbi_id(ntrial);%onsets
                    pmod(ntrial,3) = Tuning.VonMises( stim_id(ntrial),1,kappa,0,0);%amp
                    pmod(ntrial,5) = deriv(mod(stim_id(ntrial)./45+4-1,8)+1);%dsigma
                end
            end
            %                        
            pmod(:,2:3)      = nandemean(pmod(:,2:3));%mean correct
            pmod(:,7)        = nandemean(pmod(:,2).^2);%poly expansion time;
            
            pmod(:,5)        = nandemean(pmod(:,5));%dsigma
            pmod(:,4)        = pmod(:,2).*pmod(:,3);%time x amp            
            pmod(:,6)        = pmod(:,2).*pmod(:,5);%time x dsigma
            pmod(:,8)        = pmod(:,7).*pmod(:,3);%time2 x amp
            pmod(:,9)        = pmod(:,7).*pmod(:,5);%time2 x dsigma
            pmod(:,2:end)    = nanzscore(pmod(:,2:end));
            %%                        
            cond             = [];
            %all valid trials
            i                = stim_id < 200;
            cond(1).name     = 'onsets';
            cond(1).onset    = scan(i);
            cond(1).duration = zeros(1,sum(i));
            cond(1).tmod     = 0;
            cond(1).pmod     = struct('name',{'time' 'amp' 'ampxtime' 'time2' ' ampxtime2' 'dyds' 'dydsxtime' 'dydsxtime2'},'param',{ pmod(i,2) pmod(i,3) pmod(i,4) pmod(i,7) pmod(i,8) pmod(i,5) pmod(i,6) pmod(i,9)},'poly',{1 1 1 1 1 1 1 1});
%             bar([cond(1).onset]',[cond(1).pmod(4).param]');                        
            %all the rest invalid trials            
            i                = stim_id > 200 & stim_id < 800;
            cond(2).name     = 'onsets_invalid';
            cond(2).onset    = scan(i);
            cond(2).duration = zeros(1,sum(i));
            cond(2).tmod     = 0;
            cond(2).pmod     = struct('name',{'time' 'amp' 'ampxtime' 'time2' 'ampxtime2' 'dyds' 'dydsxtime' 'dydsxtime2'},'param',{pmod(i,2) pmod(i,3) pmod(i,4) pmod(i,7) pmod(i,8) pmod(i,5) pmod(i,6) pmod(i,9)},'poly',{1 1 1 1 1 1 1 1});
            %null trials            
            i                = stim_id == 1000;
            cond(3).name     = 'null';
            cond(3).onset    = scan(i);
            cond(3).duration = zeros(1,sum(i));
            cond(3).tmod     = 0;
            cond(3).pmod     = struct('name',{},'param',{},'poly',{});            
            %add also the oddtrials
            i                = stim_id == 1002;
            cond(4).name     = 'odd';
            cond(4).onset    = scan(i);
            cond(4).duration = zeros(1,sum(i));
            cond(4).tmod     = 0;
            cond(4).pmod     = struct('name',{},'param',{},'poly',{});
            
            model_path       = self.path_model(1,model_num);
            model_dir        = fileparts(model_path);
            if ~exist(model_dir);mkdir(model_dir);end
            save(model_path,'cond');
        end          
        function analysis_CreateModel08(self)
            %chebyshev polynomial expansion
            %
            
            %this initial part is exactly the same as model_02
            model_num      = 8;
            L              = self.get_log(1);%get the log file.
            stim_onsets    = L(:,2) == 3;%all stim events.
            scan_onsets    = L((L(:,2) == 0),1);%all scan events.
            ucs_events     = L(:,2) == 5;%all ucs events
            odd_events     = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi            = L(:,end);%microblock identity for all events.
            %now we have to discard mbi where stimuli occured during a
            %transition. This can be achieved by search mbi indices that
            %occured between scans 910-911 and 1805-1806
            % micro blocks during transition:
            transition_events = (L(:,1)>scan_onsets(910) & L(:,1)<scan_onsets(911))|(L(:,1)>scan_onsets(1805) & L(:,1)<scan_onsets(1806));
            %
            % find stim events which are appearing in a microblock where
            % there is an UCS (number 5);
            i                 = ismember(mbi,mbi(ucs_events|odd_events|transition_events))&stim_onsets;
            %
            %add 500 degrees to "bad stimulus events".
            L(i,3)            = L(i,3)+360;
            %some more finetuning here basically for future convenience on labelling of conditions.
            i                 = L(:,3) ==1360;%shifted null trials: put them back to 1000, we want to model null trials together irrespective of their microblock belongance.
            L(i,3)            = 1000;
            i                 = L(:,3) ==1361;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)            = 360;
            i                 = L(:,3) ==1362;%shifted oddball trials: put them back to 1002.
            L(i,3)            = 1002;
            i                 = L(:,3) ==1000;%shifted oddball trials: put them back to 1002.                        
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);                        
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);;plot(stim_id,mbi_id,'ro','markersize',10);title(mat2str(self.id));hold off
            %%
            [pmodmat names] = self.get_chebyshev(5,6);            
            %%
            pmod    = NaN(length(stim_id),size(pmodmat,3));
            for ntrial = 1:length(stim_id)
                if stim_id(ntrial) < 1000
                    pmod(ntrial,:) = squeeze(pmodmat(mbi_id(ntrial),mod(stim_id(ntrial)./45+4-1,8)+1,:));
                end
            end
            %                        
            pmod    = nanzscore(pmod(:,2:end));
            names   = names(:,2:end);
            %%                        
            cond             = [];
            %all valid trials
            i                = stim_id < 200;
            cond(1).name     = 'onsets';
            cond(1).onset    = scan(i);
            cond(1).duration = zeros(1,sum(i));
            cond(1).tmod     = 0;
            cond(1).pmod     = struct('name',names,'param',num2cell(pmod(i,:),[1 size(pmod,2)]),'poly',num2cell(ones(1,size(pmod,2)),[1 size(pmod,2)]));
%             bar([cond(1).onset]',[cond(1).pmod(4).param]');                        
            %all the rest invalid trials            
            i                = stim_id > 200 & stim_id < 800;
            cond(2).name     = 'onsets_invalid';
            cond(2).onset    = scan(i);
            cond(2).duration = zeros(1,sum(i));
            cond(2).tmod     = 0;
            cond(2).pmod     = struct('name',names,'param',num2cell(pmod(i,:),[1 size(pmod,2)]),'poly',num2cell(ones(1,size(pmod,2)),[1 size(pmod,2)]));
            %null trials            
            i                = stim_id == 1000;
            cond(3).name     = 'null';
            cond(3).onset    = scan(i);
            cond(3).duration = zeros(1,sum(i));
            cond(3).tmod     = 0;
            cond(3).pmod     = struct('name',{},'param',{},'poly',{});            
            %add also the oddtrials
            i                = stim_id == 1002;
            cond(4).name     = 'odd';
            cond(4).onset    = scan(i);
            cond(4).duration = zeros(1,sum(i));
            cond(4).tmod     = 0;
            cond(4).pmod     = struct('name',{},'param',{},'poly',{});
            
            model_path       = self.path_model(1,model_num);
            model_dir        = fileparts(model_path);
            if ~exist(model_dir);mkdir(model_dir);end
            save(model_path,'cond');
        end          
        
        function analysis_CreateModel09(self)
            %zernike polynomial expansion
            %
            
            %this initial part is exactly the same as model_02
            model_num      = 9;
            L              = self.get_log(1);%get the log file.
            stim_onsets    = L(:,2) == 3;%all stim events.
            scan_onsets    = L((L(:,2) == 0),1);%all scan events.
            ucs_events     = L(:,2) == 5;%all ucs events
            odd_events     = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi            = L(:,end);%microblock identity for all events.
            %now we have to discard mbi where stimuli occured during a
            %transition. This can be achieved by search mbi indices that
            %occured between scans 910-911 and 1805-1806
            % micro blocks during transition:
            transition_events = (L(:,1)>scan_onsets(910) & L(:,1)<scan_onsets(911))|(L(:,1)>scan_onsets(1805) & L(:,1)<scan_onsets(1806));
            %
            % find stim events which are appearing in a microblock where
            % there is an UCS (number 5);
            i                 = ismember(mbi,mbi(ucs_events|odd_events|transition_events))&stim_onsets;
            %
            %add 500 degrees to "bad stimulus events".
            L(i,3)            = L(i,3)+360;
            %some more finetuning here basically for future convenience on labelling of conditions.
            i                 = L(:,3) ==1360;%shifted null trials: put them back to 1000, we want to model null trials together irrespective of their microblock belongance.
            L(i,3)            = 1000;
            i                 = L(:,3) ==1361;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)            = 360;
            i                 = L(:,3) ==1362;%shifted oddball trials: put them back to 1002.
            L(i,3)            = 1002;
            i                 = L(:,3) ==1000;%shifted oddball trials: put them back to 1002.                        
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);                        
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);;plot(stim_id,mbi_id,'ro','markersize',10);title(mat2str(self.id));hold off
            %%
            [pmodmat, names] = self.get_zernike;
            %%
            pmod    = NaN(length(stim_id),size(pmodmat,3));
            for ntrial = 1:length(stim_id)
                if stim_id(ntrial) < 1000
                    pmod(ntrial,:) = squeeze(pmodmat(mbi_id(ntrial),mod(stim_id(ntrial)./45+4-1,8)+1,:));
                end
            end
            %                        
            pmod    = nanzscore(pmod(:,2:end));
            names   = names(:,2:end);
            %%                        
            cond             = [];
            %all valid trials
            i                = stim_id < 200;
            cond(1).name     = 'onsets';
            cond(1).onset    = scan(i);
            cond(1).duration = zeros(1,sum(i));
            cond(1).tmod     = 0;
            cond(1).pmod     = struct('name',names,'param',num2cell(pmod(i,:),[1 size(pmod,2)]),'poly',num2cell(ones(1,size(pmod,2)),[1 size(pmod,2)]));
%             bar([cond(1).onset]',[cond(1).pmod(4).param]');                        
            %all the rest invalid trials            
            i                = stim_id > 200 & stim_id < 800;
            cond(2).name     = 'onsets_invalid';
            cond(2).onset    = scan(i);
            cond(2).duration = zeros(1,sum(i));
            cond(2).tmod     = 0;
            cond(2).pmod     = struct('name',names,'param',num2cell(pmod(i,:),[1 size(pmod,2)]),'poly',num2cell(ones(1,size(pmod,2)),[1 size(pmod,2)]));
            %null trials            
            i                = stim_id == 1000;
            cond(3).name     = 'null';
            cond(3).onset    = scan(i);
            cond(3).duration = zeros(1,sum(i));
            cond(3).tmod     = 0;
            cond(3).pmod     = struct('name',{},'param',{},'poly',{});            
            %add also the oddtrials
            i                = stim_id == 1002;
            cond(4).name     = 'odd';
            cond(4).onset    = scan(i);
            cond(4).duration = zeros(1,sum(i));
            cond(4).tmod     = 0;
            cond(4).pmod     = struct('name',{},'param',{},'poly',{});
            
            model_path       = self.path_model(1,model_num);
            model_dir        = fileparts(model_path);
            if ~exist(model_dir);mkdir(model_dir);end
            save(model_path,'cond');
        end          
        
        function analysis_CreateModel06(self)
            %will generate single trial based conditions sorted by
            %condition, will not discard mbi with UCS, so that once the
            %model is estimated I can simply plot the whole timexcond
            %microblock. However it discards microblocks during a transition period. 
            %%
            model_num       = 6;
            model_path      = self.path_model(1,model_num);
            L               = self.get_log(1);%get the log file.
            stim_onsets     = L(:,2) == 3;%all stim events.
            scan_onsets     = L((L(:,2) == 0),1);%all scan events.
            ucs_events      = L(:,2) == 5;%all ucs events
            odd_events      = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi             = L(:,end);%microblock identity for all events.
            %now we have to DELETE MBI's where stimuli occured during a transition.              
                        
            %
            % find stim events which are appearing in a microblock where
            % there is an UCS (number 5);
            i                 = ismember(mbi,[22 23 45 46])&stim_onsets;%transition mbi are 22, 23, 45, 46
            L(i,3)            = L(i,3)+1500;
%             i               = ismember(mbi,[10]);%&stim_onsets;%transition blocks;
%             L(i,:)          = [];
            
            i               = L(:,3) ==1001;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)          = 0;
            i               = L(:,3) ==1002;%shifted oddball trials: put them back to 1002.
            L(i,3)          = 0;                        
            %%
            %from this point on it is the same as            
            [scan,stim_id,stim_mbi]  = self.analysis_StimTime2ScanUnit(L);
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);plot(stim_id,stim_mbi,'ro','markersize',10);title(mat2str(self.id));hold off
            %we dont want to add additional collinearity on the DM by
            %modelling the null trials, so we just remove it
            invalid          = stim_id > 500;
            stim_id(invalid) = [];
            scan(invalid)    = [];            
            cond = [];       
            c = 0;
            for ncond = unique(stim_id)'
                for trial = find(stim_id == ncond)';%make a regressors for each single stimulus
                    c                 = c +1;
                    cond(c).name     = sprintf('cond: %i,mbi: %i',ncond,ceil(c./9));
                    cond(c).onset    = scan(trial);
                    cond(c).duration = 0;
                    cond(c).tmod     = 0;
                    cond(c).pmod     = struct('name',{},'param',{},'poly',{});
                end
            end
            %%
            if ~exist(fileparts(model_path));
                mkdir(fileparts(model_path));
            end            
            save(model_path,'cond');
        end    
        function analysis_CreateModel07(self)
            %Based on model3, However instead of a constant Gaussian, it
            %uses the behavioral ratings as pmod.
            %
            
            %this initial part is exactly the same as model_02
            model_num      = 7;            
            L              = self.get_log(1);%get the log file.
            stim_onsets    = L(:,2) == 3;%all stim events.
            scan_onsets    = L((L(:,2) == 0),1);%all scan events.
            ucs_events     = L(:,2) == 5;%all ucs events
            odd_events     = L(:,3) == 1002;%recover oddball events from the stim id as they are not logged in a specific channel.
            %
            mbi            = L(:,end);%microblock identity for all events.
            %now we have to discard mbi where stimuli occured during a
            %transition. This can be achieved by search mbi indices that
            %occured between scans 910-911 and 1805-1806
            % micro blocks during transition:
            transition_events = (L(:,1)>scan_onsets(910) & L(:,1)<scan_onsets(911))|(L(:,1)>scan_onsets(1805) & L(:,1)<scan_onsets(1806));
            %
            % find stim events which are appearing in a microblock where
            % there is an UCS (number 5);
            i                 = ismember(mbi,mbi(ucs_events|odd_events|transition_events))&stim_onsets;
            %
            %add 500 degrees to "bad stimulus events".
            L(i,3)            = L(i,3)+360;
            %some more finetuning here basically for future convenience on labelling of conditions.
            i                 = L(:,3) ==1360;%shifted null trials: put them back to 1000, we want to model null trials together irrespective of their microblock belongance.
            L(i,3)            = 1000;
            i                 = L(:,3) ==1361;%shifted ucs trials: put them to 500, which is the cond_id for CS+ with shock.
            L(i,3)            = 360;
            i                 = L(:,3) ==1362;%shifted oddball trials: put them back to 1002.
            L(i,3)            = 1002;
            i                 = L(:,3) ==1000;%shifted oddball trials: put them back to 1002.                        
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);            
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));hold on;
            figure(1);;plot(stim_id,mbi_id,'ro','markersize',10);title(mat2str(self.id));hold off
            %%            
            pmod    = NaN(length(stim_id),4);
            %collect ratings to be used as pmod.
            Ratings = zscore(self.rating.y_mean);%zscored ratings.            
            for ntrial = 1:length(stim_id)
                if stim_id(ntrial) < 1000
                    pmod(ntrial,1) = 1;%constant term
                    pmod(ntrial,2) = mbi_id(ntrial);%time
                    %pmod amplitude is obtained from the ratings.
                    current        = mod(unique(sort(stim_id(ntrial)/45+4))-1,8)+1;%current face index
                    pmod(ntrial,3) = Ratings( current);
                end
            end
            %
            pmod(:,2:3)      = nandemean(pmod(:,2:3));                        
            pmod(:,4)        = pmod(:,2).*pmod(:,3);%time x amp                        
            pmod(:,2:end)    = nanzscore(pmod(:,2:end));
            %%                        
            cond             = [];
            %all valid trials
            i                = stim_id < 200;
            cond(1).name     = 'onsets';
            cond(1).onset    = scan(i);
            cond(1).duration = zeros(1,sum(i));
            cond(1).tmod     = 0;
            cond(1).pmod     = struct('name',{'time' 'amp' 'ampxtime' },'param',{ pmod(i,2) pmod(i,3) pmod(i,4) },'poly',{1 1 1 });
%             bar([cond(1).onset]',[cond(1).pmod(4).param]');                        
            %all the rest invalid trials            
            i                = stim_id > 200 & stim_id < 800;
            cond(2).name     = 'onsets_invalid';
            cond(2).onset    = scan(i);
            cond(2).duration = zeros(1,sum(i));
            cond(2).tmod     = 0;
            cond(2).pmod     = struct('name',{'time' 'amp' 'ampxtime' },'param',{pmod(i,2) pmod(i,3) pmod(i,4) },'poly',{1 1 1 });
            %null trials            
            i                = stim_id == 1000;
            cond(3).name     = 'null';
            cond(3).onset    = scan(i);
            cond(3).duration = zeros(1,sum(i));
            cond(3).tmod     = 0;
            cond(3).pmod     = struct('name',{},'param',{},'poly',{});            
            %add also the oddtrials
            i                = stim_id == 1002;
            cond(4).name     = 'odd';
            cond(4).onset    = scan(i);
            cond(4).duration = zeros(1,sum(i));
            cond(4).tmod     = 0;
            cond(4).pmod     = struct('name',{},'param',{},'poly',{});
            
            model_path       = self.path_model(1,model_num);
            model_dir        = fileparts(model_path);
            if ~exist(model_dir);mkdir(model_dir);end
            save(model_path,'cond');
        end        
        
        function [out] = fit_pmf(self,varargin)
            %will load the pmf fit (saved in runXXX/pmf) if computed other
            %wise will read the raw pmf data (saved in runXXX/stimulation)
            %and compute a fit.
            
            if isempty(self.pmf)
                fprintf('PMF data doesn''t exist for subject %i...\nNot proceeding with the fit...\n',self.id);                
                out = [];
                return           
            elseif exist(self.path_data(2,'pmf')) && isempty(varargin)%Fit has been found, and force not required: just load it.
                %load directly or
                load(self.path_data(2,'pmf'));
%                 fprintf('PMF Fit found and loaded successfully for subject %i...\n',self.id);                
            elseif ~isempty(varargin) || ~exist(self.path_data(2,'pmf'))%force required or not yet computed: just compute it.
                
                %compute and save it.
                fprintf('Fitting PMF...\n')
                % define a search grid
                searchGrid.alpha  = linspace(0,100,10);    %structure defining grid to
                searchGrid.beta   = 10.^[-1:0.1:1];         %search for initial values
                searchGrid.gamma  = linspace(0,0.5,10);
                searchGrid.lambda = linspace(0,0.1,10);
                paramsFree        = [1 1 1 1];
                PF                = @PAL_Weibull;                
                options                      = PAL_minimize('options');
                options.MaxIter              = 10.^3;
                options.MaxFunEvals          = 10.^3;
                options.Display              = 'On';
                options.ToX                  = -10.^3;
                options.TolFun               = -10.^3;
                xlevels_HD                   = linspace(min(self.pmf.x),max(self.pmf.x),100);
                tchain                       = size(self.pmf.NumPos,2);
                for chain = 1:tchain                    
                    %fit the function using PAL
                    [params, LL, exitflag]       = PAL_PFML_Fit(self.pmf.x, self.pmf.NumPos(:,chain), self.pmf.OutOfNum(:,chain), searchGrid, paramsFree, PF,'lapseLimits',[0 .5],'guessLimits',[0 .5]);
                    out.params(chain,:)          = params;
                    out.LL(chain,:)              = LL;
                    out.exitflag(chain,:)        = exitflag;
                    out.y(chain,:)               = PF(params,xlevels_HD);
                    out.x                        = xlevels_HD;                                        
                end
                save(self.path_data(2,'pmf'),'out')
            end
        end
        function [out] = fit_rating(self)
            %will load the rating fit (saved in runXXX/rating) if computed other
            %wise will read the raw ratingdata (saved in runXXX/stimulation)
            %and compute a fit.
            
            fun        = self.selected_fitfun;;%vM function
            force      = 0;%repeat the analysis or load from cache            
            write_path = sprintf('%s/midlevel/rating_fun_%i.mat',self.pathfinder(self.id,1),fun);            
            
            if exist(write_path) && force ==0
                %load directly or
                load(write_path);
%                 fprintf('Rating Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif force == 1 || ~exist(write_path)
                %compute and save it.
                fprintf('Fitting Ratings...\n')                
                R                            = self.rating;
                %adapt to what Tuning.m wants to have.
                R.y                          = R.y(:)';
                R.x                          = R.x(:)';
                R.ids                        = R.ids(1);                
                %create a tuning object and make a single subject fit
                T                            = Tuning(R);                
                T.SingleSubjectFit(fun);                               
                %prepare data for outputting.
                out.params                   = T.fit_results{fun}.params(1,:);                
                out.LL                       = T.fit_results{fun}.pval(1,:);
                out.exitflag                 = T.fit_results{fun}.ExitFlag(1,:);
                out.y_hd                     = T.fit_results{fun}.y_fitted_HD(1,:);
                out.x_hd                     = T.fit_results{fun}.x_HD(1,:);
                out.y                        = T.fit_results{fun}.y_fitted(1,:);
                out.x                        = T.fit_results{fun}.x(1,:);
                out.fitfun                   = T.fit_results{fun}.fitfun;                
                %
                if fun == 8%if vM, then transform kappa to FWHM.
                    fprintf('Kappa to FWHM transformation + absolute(peakshift).\n');
                    out.params(:,2)            = vM2FWHM(out.params(:,2));
                    %out.params(:,3)            = abs(out.params(:,3));
                end
                save(write_path,'out')
            end
            
        end        
        function [out] = fit_facecircle(self,partition)
            %will fit function FUN to facecircle data with PARTITION
            %partitions.
            
            force = 0;%repeat the analysis or load from cache            
            fun   = self.selected_fitfun;
            write_path = sprintf('%s/midlevel/facecircle_fun_%i_partition_%i.mat',self.pathfinder(self.id,1),fun,partition);
            if exist(write_path) && force == 0
                %load directly or
                load(write_path);
%                 fprintf('Rating Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif force == 1 || ~exist(write_path)
                %compute and save it.
                fprintf('Fitting Ratings...\n')                
                R                            = self.get_facecircle(partition);
                %adapt to what Tuning.m wants to have.                
                R.y                          = self.circconv2(R.countw,[1 1]/2);;                                
                %create a tuning object and make a single subject fit
                T                            = Tuning(R);
                
                T.SingleSubjectFit(fun);
                %prepare data for outputting.
                out.params                   = T.fit_results{fun}.params; 
                out.LL                       = T.fit_results{fun}.pval;
                out.exitflag                 = T.fit_results{fun}.ExitFlag;
                out.y                        = T.fit_results{fun}.y_fitted_HD;
                out.x                        = T.fit_results{fun}.x_HD;
                out.fitfun                   = T.fit_results{fun}.fitfun;
                if fun == 8%if vM, then transform kappa to FWHM.
                    fprintf('Kappa to FWHM transformation + absolute(peakshift).\n');
                    out.params(:,2)            = vM2FWHM(out.params(:,2));
                    %out.params(:,3)            = abs(out.params(:,3));
                end
                save(write_path,'out');
            end
            
        end        
        function [out] = fit_pupil(self,fun)
            %will load the rating fit (saved in runXXX/rating) if computed
            %otherwise will read the raw ratingdata (saved in
            %runXXX/stimulation) and compute a fit.
            
            force      = 1;%repeat the analysis or load from cache            
            write_path = sprintf('%s/midlevel/pupil_fun_%i.mat',self.pathfinder(self.id,1),fun);            
            
            if exist(write_path) && force ==0
                %load directly or
                load(write_path);
%                 fprintf('Rating Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif force == 1 || ~exist(write_path)
                %compute and save it.
                fprintf('Fitting Pupil...\n')                
                dummy                        = self.pupil.get_singletrials(self.id,'clean');
                dummy                        = dummy(end);                
                R.y                          = Vectorize(dummy.mean(:,1:8))';
                R.x                          = Vectorize(dummy.x(:,1:8))';
                figure(200);imagesc(dummy.mean(:,1:8));drawnow;
                figure(201);plot(R.x(:),R.y(:),'o');hold on;plot(unique(R.x),accumarray(R.x(~isnan(R.y(:)))'./45+4,R.y(~isnan(R.y(:))),[8 1],@mean),'r');hold off;drawnow;
                %adapt to what Tuning.m wants to have.
                invalid                      = ~isnan(R.y);
                R.y                          = R.y(invalid);
                R.x                          = R.x(invalid);
                R.ids                        = self.id;
                %create a tuning object and make a single subject fit
                T                            = Tuning(R);                
                T.SingleSubjectFit(fun);                               
                %prepare data for outputting.
                out.params                   = T.fit_results{fun}.params(1,:);                
                out.LL                       = T.fit_results{fun}.pval(1,:);
                out.exitflag                 = T.fit_results{fun}.ExitFlag(1,:);
                out.y                        = T.fit_results{fun}.y_fitted_HD(1,:);
                out.x                        = T.fit_results{fun}.x_HD(1,:);
                out.fitfun                   = T.fit_results{fun}.fitfun;                
                %
                if fun == 8%if vM, then transform kappa to FWHM.
                    fprintf('Kappa to FWHM transformation + absolute(peakshift).\n');
                    out.params(:,2)            = vM2FWHM(out.params(:,2));
                    out.params(:,3)            = abs(out.params(:,3));
                end
                save(write_path,'out')
            end
            
        end        
        
        function [out] = fit_scr(self)
            %will load the rating fit (saved in runXXX/rating) if computed other
            %wise will read the raw ratingdata (saved in runXXX/stimulation)
            %and compute a fit.
            
            fun        = self.selected_fitfun;%
            force      = 0;%repeat the analysis or load from cache            
            write_path = sprintf('%s/midlevel/scr_fun_%i.mat',self.pathfinder(self.id,1),fun);            
            
            if exist(write_path) && force ==0
                %load directly or
                load(write_path);
%                 fprintf('SCR Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif force == 1 || ~exist(write_path)
                %compute and save it.
                fprintf('Fitting scr...\n')                
                R                            = self.scr;
                %adapt to what Tuning.m wants to have.
                R.y                          = R.y(:)';
                R.x                          = R.x(:)';
                R.ids                        = R.ids(1);                
                %create a tuning object and make a single subject fit
                T                            = Tuning(R);                
                T.SingleSubjectFit(fun);                               
                %prepare data for outputting.
                out.params                   = T.fit_results{fun}.params(1,:);                
                out.LL                       = T.fit_results{fun}.pval(1,:);
                out.exitflag                 = T.fit_results{fun}.ExitFlag(1,:);
                out.y_hd                     = T.fit_results{fun}.y_fitted_HD(1,:);
                out.x_hd                     = T.fit_results{fun}.x_HD(1,:);
                out.y                        = T.fit_results{fun}.y_fitted(1,:);
                out.x                        = T.fit_results{fun}.x(1,:);
                out.fitfun                   = T.fit_results{fun}.fitfun;                
                %
                if fun == 8%if vM, then transform kappa to FWHM.
                    fprintf('Kappa to FWHM transformation + absolute(peakshift).\n');
                    out.params(:,2)            = vM2FWHM(out.params(:,2));
                    %out.params(:,3)            = abs(out.params(:,3));
                end
                save(write_path,'out')
            end
            
        end        
 
        function writename = brainbehavior_analysis(self,sk)
            
            
            writename = sprintf('%sbrainbehavior_%02dmm.nii',self.path_midlevel(1),sk);
            if exist(writename) == 0;
                %load mumfordian analysis
                self.default_model_name = 'chrf_0_0_mumfordian';
                vol  = spm_vol(self.path_beta(1,0,''));
                tobewritten = vol(1);
                vol  = spm_smoothto16bit(vol,sk);
                
                data = spm_read_vols(vol);
                % reshape all to transform it to Condition x Voxel
                s    =  size(data);
                data = reshape(data,s(1)*s(2)*s(3),s(4));
                data = data';
                data = reshape(data,[65 11 size(data,2)]);;
                data = data(self.mbi_valid,1:8,:);
                data = squeeze(mean(data));
                %             data = permute(data,[2 1 3]);
                %             data = reshape(data,[8 size(data,2)*size(data,3)]);
                %%
                var               = self.rating.y_mean;
                r                 = corr(data,var(:));
                r                 = reshape(r,[s(1) s(2) s(3)]);                
                tobewritten.fname = writename;
                spm_write_vol(tobewritten,r);
                self.VolumeNormalize(writename);%normalize them ('w_' will be added);
            else
                cprintf([1 0 0],'Already computed...\n');
            end
        end
        function [volume]=get_brainbehavior_analysis(self)
            %[volume]=get_brainbehavior_analysis(self)
            %
            writename = sprintf('%sbrainbehavior.nii',self.path_midlevel(1));
            volume    = spm_read_vol(spm_vol(writename));
        end
    end
end
