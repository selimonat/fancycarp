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
        csp
        csn
        scr        
        pmf_param     = [];
        rating_param  = [];
        trio_session  = [];
        ratings       = [];
        total_run     = [];
        pmf
    end
    %%
    methods
        function s = Subject(id)%constructor
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
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn = s.paradigm{s.default_run}.stim.cs_neg;
                end
                s.scr = SCR(s);
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
        function rating     = get.ratings(self)
            %returns the CS+-aligned ratings for all the runs.
            for run = unique(self.dicom2run(:))';%don't count the first run
                if isfield(self.paradigm{run}.out,'rating')
                    if ~isempty(self.paradigm{run});
                        rating(run).y      = self.paradigm{run}.out.rating';
                        rating(run).y      = circshift(rating.y,[1 4-self.csp ]);
                        rating(run).x      = repmat([-135:45:180],size(self.paradigm{run}.out.rating,2),1);
                        rating(run).ids    = repmat(self.id,size(self.paradigm{run}.out.rating,2),size(self.paradigm{run}.out.rating,1));
                        rating(run).y_mean = mean(rating.y);
                        rating(run).y_std  = std(rating.y)./sqrt(2);                                                
                    else
                        fprintf('No rating present for this subject and run (%d) \n',nr);
                    end
                end
            end
        end
        function out        = get_scr(self,run,cond)
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
        function plot_facecircle(self)
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
        function [fg]=get_facecircle(self,partition)            
            %% get the masks to count the fixations as a f(cs+ distance)
            res      = self.screen_resolution;
            center   = res./2;%[y x]
            %% plot the fixation map with a mask            
            fix      = Fixmat(self.id,3);
            limits   = prctile(double(fix.start),linspace(0,100,partition+1));
            counter  = 0;
            for limit = [limits(1:end-1);limits(2:end)]
                %
                counter  = counter + 1;
                fix      = Fixmat(self.id,3);
                S        = fix.start;
                %
                valid    = (S >= limit(1))&(S < limit(2));%take only fixations which are within the limits.
                fix.replaceselection(valid);
                fix.ApplySelection;
                %
                X        = fix.x;
                Y        = fix.y;
                X        = (X-center(2));%center X and Y
                Y        = (Y-center(1));
                W        = double(abs(fix.stop - fix.start));
                Rank     = double(fix.fix);
%                 W        = W./sum(W);%weights which are proportional to fixation duration.
                %
                wedge_angle = 135:-45:-180;%this is the midpoint of the wedge.
                if length(unique(fix.trialid)) == 1;%%sanity check: there should be only one trial here.
                    [theta, rho]                = cart2pol(X,Y);%transform X and Y to
                    [theta]                     = rad2deg(theta);%angle of the fixation point in degrees
                    shift_angle                 = -mod(theta+22.5-45,360);%shift the angles so that we can get the wedge index.
                    wedge_index                 = mod(ceil(shift_angle./45)+3-1,8)+1;%moves together the same as PTB drawing
                end
                face_angle          = self.paradigm{1}.stim.circle_order(wedge_index)';%drawn condition in the wedge
                face_index          = face_angle./45+4;                
                A                   = [X;Y;double(W);theta;rho;wedge_index;wedge_angle(wedge_index);face_angle;face_index;Rank];
                %A is still in register with fixmat
                valid               = (A(5,:) > 190)&(A(5,:) < 380);%take only those fixations which are on the faces.
                A(:,~valid)         = [];%remove also fixations from the fixmat
                fix.replaceselection(valid);
                fix.ApplySelection;
                %
                [~,i]               = sort(A(6,:));%sort according to wedge_index i.e. drawing order.
                A                   = A(:,i);                
                
                [A(11,:),A(12,:)]   = pol2cart(deg2rad(A(8,:)+(A(4,:)-A(7,:))),A(5,:));%(new theta, old amp)
%                 figure(333);
%                 imagesc(A);
%                 set(gca,'ytick',1:size(A,1),'yticklabel',{'X' 'Y' 'W' 'angle' 'amp' 'wedge index' 'wedge angle' 'face angle' 'face index' 'Rank' 'NewX' 'NewY' });
%                 %                                
%                 %
%                 figure(3333);
%                 fix.kernel_fwhm = 30;
%                 subplot(2,1,1)
%                 fix.getcondmaps;
%                 M = fix.maps;
%                 imagesc(M);
%                 subplot(2,1,2);
%                 fix.duration_weight = 0;                
%                 fix.replaceXY(min(round(A(11,:)+center(2)),1023),min(round(A(12,:)+center(1)),767));
%                 fix.getcondmaps;
%                 M=fix.maps;
%                 imagesc(M)
                %
%                 figure(8);
                for nface = 1:8
                    fg.duration(counter,nface)  = sum(A(3,A(9,:) == nface));                    
                    fg.count(counter,nface)     = sum(A(9,:) == nface);
                end                
                ranks = unique(A(10,:));
                fg.Y = zeros(100,8);
                for rank = ranks(ranks<=100)
                    i                 = A(10,:) == rank;
                    fg.Y(rank,A(9,i)) = fg.Y(rank,A(9,i)) + 1;
                end
%                 figure;
%                 imagesc(fg.Y);drawnow
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
        function out        = get.pmf(self)
            %will load the raw pmf data.
            
            %first a double check for hiwi-fakups
            if self.csp == (self.paradigm{2}.stim.cs_plus)./45+1                
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
                out      = self.fit_pmf.params;
            else
                out = [];
                cprintf([1 0 0],'No pmf data for this subject.\n')
            end
        end         
        function out        = get.rating_param(self)
            %returns the parameters of the pmf fit (1 x parameter);            
            out      = [self.fit_rating.params(1,:)];
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
        function XYZvox     = get_mm2vox(self,XYZmm,vh)
            %brings points in the world space XYZmm to voxel space XYZvox
            %of the image in VH.
            XYZvox  = vh.mat\XYZmm;
            XYZvox  = unique(XYZvox','rows')';%remove repetitions
            XYZvox  = round(XYZvox);%remove decimals as these are going to be matlab indices.
        end
        function D          = get_data(self,file,mask_id)
            %will read the data specified in FILE 
            %FILE is the absolute path to a 3/4D .nii file.            
            %
            %MASK_ID is used to select voxels.
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
                78 76 85];
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
    methods %(analysis)        
        function [X,N,K]    = analysis_designmatrix(self,nrun,model_num)
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
            k                       = self.get_total_volumes(nrun);
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
            % get the filtering strcture a la spm. 
            run_borders              = [[0 910 910+895]+1;[910 910+895  self.get_total_volumes(nrun)]];
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
        function beta       = analysis_firstlevel(self,nrun,model_num,mask_id)
            %will compute beta weights "manually" without calling SPM.
            
            
            [X N K]   = self.spm_DesignMatrix(nrun,model_num);%returns the Design Matrix, Nuissance Matrix, and High-pass Filtering Matrix                        
            Y         = self.get_data(self.path_epi(nrun,'r'),mask_id);%return the realigned data.
            
%             GM        = 100;
%             g         = spm_global(spm_vol(self.path_epi(nrun)));    
%             factor    = GM./mean(g);
%             Y         = Y*factor;
            Y         = (Y-mean(Y(:)))./std(Y(:));                        
            Y         = spm_filter(K,Y);%high-pass filtering.
            %            
            DM        = [X N ones(size(X,1),1)];%append together Onsets, Nuissances and a constant
            DM        = spm_filter(K,DM);%filter also the design matrix
            DM        = spm_sp('Set',DM);
            DM        = spm_sp('x-',DM);% projector;
            beta      = DM*Y;
            beta      = beta';%(voxels x betas)
        end                
        function beta       = analysis_mumfordian(self,nrun,mask_id)
            
            %% get the nuissance and highpass filtering matrices.
            N                       = self.get_param_nuissance(nrun);            
            run_borders             = [[0 910 910+895]+1;[910 910+895  self.get_total_volumes(nrun)]];
            %
            K(1:size(run_borders,2)) = struct('HParam', self.HParam, 'row',    [] , 'RT',     self.TR ,'X0',[]);
            c = 0;            
            for b = run_borders
                c        = c + 1;
                K(c).row = b(1):b(2);
                K(c)     = spm_filter(K(c));
                K(c).X0  = [ones(length(K(c).row),1)*std(K(c).X0(:)) K(c).X0];                
            end   
            %% get the data and highpass filter right away.
            Y                        = self.get_data(self.path_epi(nrun,'r'),mask_id);%return the realigned data.
            Y                        = (Y-mean(Y(:)))./std(Y(:));                        
            Y                        = spm_filter(K,Y);%high-pass filtering.
            %% now we run across all trials and get a specific design matrix for each trial.                        
            L                        = self.get_log(1);
            [onsets,stim_id,mbi]     = self.analysis_StimTime2ScanUnit(L);            
            stim_id(stim_id<500)     = stim_id(stim_id<500)./45+4;
            stim_id(stim_id == 1000) = 9;
            stim_id(stim_id == 1001) = 10;
            stim_id(stim_id == 1002) = 11;                       
            %%
            new_cond                 = [];
            ttrial                   = length(stim_id);
            for ntrial = 1:ttrial
                new_cond{ntrial}(1).name     = num2str(stim_id(ntrial));
                new_cond{ntrial}(1).onset    = onsets(ntrial);
                new_cond{ntrial}(1).duration = 0;
                new_cond{ntrial}(1).tmod     = [];
                new_cond{ntrial}(1).pmod     = struct('name',{},'param',{},'poly',{});
                %
                new_cond{ntrial}(2).name     = 'all_other_stims';%exclude Null trials, UCS and oddball
                new_cond{ntrial}(2).onset    = setdiff(onsets(stim_id < 9),onsets(ntrial));%everything else
                new_cond{ntrial}(2).duration = 0;
                new_cond{ntrial}(2).tmod     = [];
                new_cond{ntrial}(2).pmod     = struct('name',{},'param',{},'poly',{});
                %
                new_cond{ntrial}(3).name     = 'ucs_oddballs';
                new_cond{ntrial}(3).onset    = setdiff(onsets(stim_id > 9),onsets(ntrial));%everything else;%everything else
                new_cond{ntrial}(3).duration = 0;
                new_cond{ntrial}(3).tmod     = [];
                new_cond{ntrial}(3).pmod     = struct('name',{},'param',{},'poly',{});
            end            
            %%
            
            for ntrial = 1:ttrial
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
                for i = 1:length(new_cond{ntrial});%one regressor for each condition
                    Sess.U(i).dt        = xBF.dt;%- time bin (seconds)
                    Sess.U(i).ons       = new_cond{ntrial}(i).onset;%- onsets    (in SPM.xBF.UNITS)
                    Sess.U(i).name      = {sprintf('%02d',i)};%- cell of names for each input or cause
                    %no parametric modulation here
                    Sess.U(i).dur    =  repmat(0,length(Sess.U(i).ons),1);%- durations (in SPM.xBF.UNITS)
                    Sess.U(i).P.name =  'none';
                    Sess.U(i).P.P    =  'none';
                    Sess.U(i).P.h    =  0;%- order of polynomial expansion
                    Sess.U(i).P.i    =  1;%- sub-indices of u pertaining to P
                end
                %
                k                       = self.get_total_volumes(nrun);
                SPM.xBF                 = xBF;
                SPM.nscan               = k;
                SPM.Sess                = Sess;
                SPM.Sess.U              = spm_get_ons(SPM,1);
                %
                % Convolve stimulus functions with basis functions
                [X,Xn,Fc]               = spm_Volterra(SPM.Sess.U,SPM.xBF.bf,SPM.xBF.Volterra);
                % Resample regressors at acquisition times (32 bin offset)
                X                       = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);                
                %%
                DM                     = [X N ones(size(X,1),1)];%append together Onsets, Nuissances and a constant
                DM                     = spm_filter(K,DM);%filter also the design matrix
                DM                     = spm_sp('Set',DM);
                DM                     = spm_sp('x-',DM);% projector;
                fprintf('Fitting Subject %03d''s %03dth trial (stim_id:%02d, onset:%3.5g) of %i (%2.3g percent))\n',self.id,ntrial,stim_id(ntrial),onsets(ntrial),ttrial,ntrial./ttrial*100)
                dummy                  = DM*Y;
                if ntrial == 1
                    beta = nan(max(mbi),11,length(dummy));
                end
                beta(mbi(ntrial),stim_id(ntrial),:)       = dummy(1);%(microblock,stim,voxel)                    
            end
        end
        function [result]   = analysis_roi_average(self,beta_files,mask_id)
            %returns interesting statistics with the roi MASK_ID on
            %activity maps specified in BETA_FILES.            
            D               = self.get_data(beta_files,mask_id)';
            result.mean     = mean(D);
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
           
            out = self.dir_spmmat(nrun,model_num);
            fprintf('Searching for beta images in:\n%s\n',out)
            out = spm_select('FPList',out,sprintf('^%sbeta_*',prefix'));
            if isempty(out)
                cprintf([1 0 0],'No beta images found, probably wrong prefix/run/etc is entered...\n');
                fprintf('%s\n',out)
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
            h = bar(mean(self.ratings.x),self.ratings.y_mean,1);%plot the data
            set(gca,'xtick',[0 180],'xticklabel',{'cs+' 'cs-'},'xgrid','on');%add labels
            SetFearGenBarColors(h);%colorize it.
            hold on;
            errorbar(mean(self.ratings.x),self.ratings.y_mean,self.ratings.y_std,'ko');%add error bars
            %if the fit is better than flat line, paint accordingly.
            if  (self.fit_rating.LL < -log10(.05))%plot simply a blue line if the fit is not significant.
                PlotTransparentLine(self.fit_rating.x(:),repmat(mean(self.ratings.y(:)),100,1),.35,'b','linewidth',2.5);                
            elseif self.fit_rating.fitfun(0,self.fit_rating.params) < self.fit_rating.fitfun(180,self.fit_rating.params);%if CS- > CS+, paint it blue too.
                PlotTransparentLine(self.fit_rating.x(:),self.fit_rating.y(:),.35,'b','linewidth',2.5);                                
            else
                PlotTransparentLine(self.fit_rating.x(:),self.fit_rating.y(:),.35,'r','linewidth',2.5);%this is the fit.
            end
            hold off;
            title(sprintf('id:%02d (+:%d)',self.id,self.csp),'fontsize',12);%subject and face id
            axis tight;box off;axis square;drawnow;ylim([0 10]);alpha(.5);
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
    end
    methods %(fmri analysis)
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
            
            self.default_model_name                                 = 'fir_20_20';
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
        function                                     analysis_spm_fourier(self,nrun,model_num)
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
            self.VolumeSmooth(beta_images);%smooth the native images ('s_' will be added, resulting in 's_')
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
        function analysis_CreateModel01(self,runs)
            %creates a model based on the stimulus onsets and logged pulses
            %in the log file. This model doesn't remove the microblocks
            %where an UCS event has occured but removes the transition
            %blocks. It is best suited for mumfordian single trial
            %analysis.
            model_num  = 1;
            for run = runs
                model_path               = self.path_model(run,model_num);%path to the model
                if ~exist(fileparts(model_path));
                    mkdir(fileparts(model_path));
                end
                L                        = self.get_log(run);
                %%
                [scan,stim_id]           = self.analysis_StimTime2ScanUnit(L);
                counter                  = 0;
                for current_condition = unique(stim_id(:))'
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
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));
            %from this point on it is the same as
            %self.analysis_CreateModels01
            model_path      = self.path_model(1,model_num);
            [scan,stim_id]  = self.analysis_StimTime2ScanUnit(L);
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
            %this model discards condition information and treat all
            %stimuli as if one single type, however add different pmod
            %values based on their condition labels. It relies on Model02,
            %which consists of "clean" microblocks. 
            %cond(1) = 'Valid stim onsets" based on conditions
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
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);            
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
            %this model discards condition information and treat all
            %stimuli as if one single type, however add different pmod
            %values based on their condition labels. It relies on Model02,
            %which consists of "clean" microblocks. 
            %cond(1) = 'Valid stim onsets" based on conditions
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
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);            
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
            save(model_path,'cond');
        end        
        function analysis_CreateModel05(self)
            %same as 04, but with quadratic terms on pmods.
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
            %sanity check plot
            figure(1);i=L(:,2)==3;plot(L(i,3),mbi(i),'k.','markersize',20);title(mat2str(self.id));
            %from this point on it is the same as
            %self.analysis_CreateModels01            
            [scan,stim_id,mbi_id]  = self.analysis_StimTime2ScanUnit(L);            
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
            cond(1).pmod     = struct('name',{'time' 'amp' 'ampxtime' 'time2' ' ampxtime2' 'dyds' 'dydsxtime'},'param',{ pmod(i,2) pmod(i,3) pmod(i,4) pmod(i,7) pmod(i,8) pmod(i,5) pmod(i,6)},'poly',{1 1 1 1 1 1 1});
%             bar([cond(1).onset]',[cond(1).pmod(4).param]');                        
            %all the rest invalid trials            
            i                = stim_id > 200 & stim_id < 800;
            cond(2).name     = 'onsets_invalid';
            cond(2).onset    = scan(i);
            cond(2).duration = zeros(1,sum(i));
            cond(2).tmod     = 0;
            cond(2).pmod     = struct('name',{'time' 'amp' 'ampxtime' 'time2' 'ampxtime2' 'dyds' 'dydsxtime'},'param',{pmod(i,2) pmod(i,3) pmod(i,4) pmod(i,7) pmod(i,8) pmod(i,5) pmod(i,6)},'poly',{1 1 1 1 1 1 1});
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
        function [out] = fit_rating(self,varargin)
            %will load the rating fit (saved in runXXX/rating) if computed other
            %wise will read the raw ratingsdata (saved in runXXX/stimulation)
            %and compute a fit.
            
            if exist(self.path_data(1,'rating')) && isempty(varargin)
                %load directly or
                load(self.path_data(1,'rating'));
%                 fprintf('Rating Fit found and loaded successfully for subject %i...\n',self.id);
                
            elseif ~isempty(varargin) || ~exist(self.path_data(2,'rating'))
                %compute and save it.
                fprintf('Fitting Ratings...\n')                
                R                            = self.ratings;
                %adapt to what Tuning.m wants to have.
                R.y                          = R.y(:)';
                R.x                          = R.x(:)';
                R.ids                        = R.ids(1);                
                %create a tuning object and make a single subject fit
                T                            = Tuning(R);
                fun = 8;
                T.SingleSubjectFit(fun);                               
                %prepare data for outputting.
                out.params                   = T.fit_results.params(1,:);                
                out.LL                       = T.fit_results.pval(1,:);
                out.exitflag                 = T.fit_results.ExitFlag(1,:);
                out.y                        = T.fit_results.y_fitted_HD(1,:);
                out.x                        = T.fit_results.x_HD(1,:);
                out.fitfun                   = T.fit_results.fitfun;
                if fun == 8%if vM, then transform kappa to FWHM.
                    fprintf('Kappa to FWHM transformation + absolute(peakshift).\n');
                    out.params(2)            = vM2FWHM(out.params(2));
                    out.params(3)            = abs(out.params(3));
                end
            end
            save(self.path_data(1,'rating'),'out')
        end        
    end
end
