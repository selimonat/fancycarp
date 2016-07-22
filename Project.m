classdef Project < handle
    % This is the PROJECT object that other objects will
    % be a child of. Here enters all the project specific data (e.g. the
    % subject ids) and methods (for e.g. getting paths from the dicom
    % server).
    %
    % The first set of properties has to be entered by hand. For example
    % TRIO_SESSIONS should be entered manually for your experiment. The
    % other properties drive from these. 
    %
    % To create the standard project structure you will need to call the
    % CreateFolderHierarchy method. Based on dicom2run, trio_sessions,
    % data_folders properties, the folder hierarchy will be created. These 
	% folders will later be populated by dump_epi and dump_hr methods. 
    %
    % Once you data is nicely loaded into the subject/run/ structure,
    % SanityCheck method comes very handy to ensure that your project
    % folders have all the size they should be.
    %
    % Don't forget to add spm, smr to your path, they are defined below,
    % however automatic changing of path names is not recommended.
	% 
	% Current methods:
	%
    % CreateFolderHierarchy: Will create a folder hierarchy to store data
	% SanityCheck: Will bar-plot the content of your folders
    % DicomDownload: connects to dicom server and dumps the epi and hr data
    % DicomTo4D: 3D 2 4D conversion
    % MergeTo4D: 3D 2 4D conversion
    % ConvertDicom: convert dicoms
    % dicomserver_request: unix-specific commands to talk with dicomserser.
    % dicomserver_paths: unix-specific commands to talk with dicomserser.
    % dartel_templates: location of dartel templates
	% SecondLevel_ANOVA: conducts secondlevel analysis.
    % VolumeGroupAverage: can average lots of images.
    % VolumeSmooth: smooths volumes
    % RunSPMJob: runs spm jobs.
    %
    % Feel free to improve this help section.
    %
    % 
    
    properties (Hidden, Constant)%adapt these properties for your project
        %All these properties MUST BE CORRECT and adapted to one owns
        %project

        path_project          = '/projects/fearamy/data/';        
        path_spm              = '/common/apps/spm12-6685/';        
        trio_sessions         = {  '' '' '' '' 'TRIO_17468' 'TRIO_17476' 'TRIO_17477' 'TRIO_17478' 'TRIO_17479' 'TRIO_17480' 'TRIO_17481' 'TRIO_17482' 'TRIO_17483' 'TRIO_17484' 'TRIO_17485' 'TRIO_17486' 'TRIO_17487' 'TRIO_17488' 'TRIO_17514' 'TRIO_17515' 'TRIO_17516' 'TRIO_17517'  'TRIO_17520' 'TRIO_17521' 'TRIO_17522' 'TRIO_17523' 'TRIO_17524' 'TRIO_17525' 'TRIO_17526' 'TRIO_17527' 'TRIO_17557' 'TRIO_17558' 'TRIO_17559' 'TRIO_17560'  'TRIO_17563' 'TRIO_17564' 'TRIO_17565' 'TRIO_17566' 'TRIO_17567' 'TRIO_17568' 'TRIO_17569' 'TRIO_17570' 'TRIO_17571' 'TRIO_17572'};
        dicom_serie_selector  = {  [] [] []   []      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [5 6 7]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]       [3 4 5]       [3 4 5]      [3 4 5]      [3 4 5]    [3 4 5]       [3 4 5]       [3 4 5]     [3 4 5]     [4 5 6]       [3 4 5]      [3 4 5]     [3 4 5]       [3 4 5]      [3 4 5]        [3 4 5]     [3 4 5]       [3 4 5]      [3 4 5]       [3 4 5]     [3 4 5]     [4 5 6]      [3 4 5]    };
        %this is necessary to tell matlab which series corresponds to which
        %run (i.e. it doesn't always corresponds to different runs)
        dicom2run             = [1 1 1];%how to distribute TRIO sessiosn to folders.
        data_folders          = {'eye' 'midlevel' 'mrt' 'scr' 'stimulation' 'pmf' 'rating'};%if you need another folder, do it here.
        TR                    = 0.99;                
        HParam                = 128;%parameter for high-pass filtering
        surface_wanted        = 0;%do you want CAT12 toolbox to generate surfaces during segmentation (0/1)                
        smoothing_factor      = 4;%how many mm images should be smoothened when calling the SmoothVolume method        
        path_smr              = sprintf('%s%ssmrReader%s',fileparts(which('Project')),filesep,filesep);%path to .SMR importing files in the fancycarp toolbox.
        gs                    = [5 6 7 8 9 12 14 16 17 18 19 20 21 24 25 27 28 30 32 34 35 37 39 40 41 42 43 44];
        bs                    = [10 11 13 15 22 23 26 29 31 33 36 38];
    end
    properties (Constant,Hidden) %These properties drive from the above, do not directly change them.
        tpm_dir               = sprintf('%stpm/',Project.path_spm); %path to the TPM images, needed by segment.         
        path_second_level     = sprintf('%sspm/',Project.path_project);%where the second level results has to be stored        
		current_time          = datestr(now,'hh:mm:ss');
        subject_indices       = find(cellfun(@(x) ~isempty(x),Project.trio_sessions));% will return the index for valid subjects (i.e. where TRIO_SESSIONS is not empty). Useful to setup loop to run across subjects.
        PixelPerDegree        = 29;
        screen_resolution     = [768 1024];
    end    
    properties (Hidden)
        atlas2mask_threshold  = 20;%where ROI masks are computed, this threshold is used.        
    end
	methods
        function DU         = SanityCheck(self,runs,measure,varargin)
            %DU = SanityCheck(self,runs,measure,varargin)
			%will run through subject folders and will plot their disk
            %space. Use a string in VARARGIN to focus only on a subfolder.
            %MEASURE has to be 'size' or 'amount', for disk usage and
            %number of files, respectively. This only works on Unix systems.
            cd(self.path_project);%
            total_subjects = length(Project.trio_sessions);
            DU = nan(total_subjects,length(runs));
            warning('off','all');
            for ns = 1:total_subjects
                for nr = runs
                    fprintf('Requesting folder size for subject %03d and run %03d\n',ns,nr);
                    subject_run = self.pathfinder(ns,nr);
                    if ~isempty(subject_run)
                        fullpath = fullfile(subject_run,varargin{:});                    
                        if strcmp(measure,'size')
                            [a b]       = system(sprintf('/usr/bin/du -cd0 %s | grep -e total | grep -oh -e [0-9]*',fullpath));
                             DU(ns,nr+1)= str2double(b)./1024;
                        elseif strcmp(measure,'amount')
                            [a b]       = system(sprintf('/usr/bin/find %s | wc -l',fullpath));
                            DU(ns,nr+1) = str2double(b);
                        else
                            fprint('Please enter a valid MEASURE...\n');
                            return
                        end
                    else
                        fprintf('Folder doesn''t exist: %s\n',subject_run);
                    end
                end                
            end
            %%
            if isempty(varargin)
                varargin{1} = '*';
            end
            n = 0;
            for nr = runs
                n = n + 1;
                subplot(length(runs),1,n)
                bar(DU(:,n));ylabel('MB or #','fontsize',10);xlabel('Subjects','fontsize',10);box off                
                title(sprintf('Subfolder: %s Run: %i',varargin{1},nr),'fontsize',10);
            end
            warning('on','all');
        end
        function              DicomDownload(self,source,destination)
            % Will download all the dicoms, convert them and merge them to
            % 4D.
            
            if ~ismac & ~ispc
                %proceeds with downloading data if we are NOT a mac not PC
                
                fprintf('%s:\n','DicomDownload');
                %downloads data from the dicom dicom server
                if exist(destination) == 0
                    fprintf('The folder %s doesn''t exist yet, will create it...\n',destination)
                    mkdir(destination)
                end
                if ~isempty(source) & ~isempty(destination)%both paths should not be empty.
                    fprintf('Calling system''s COPY function to dump the data...%s\nsource:%s\ndestination:%s\n',self.current_time,source,destination)
                    a            = system(sprintf('cp -r %s/* %s',source,destination));
                    if a ~= 0
                        fprintf('There was a problem while dumping...try to debug it now\n');
                        keyboard
                    else
                        fprintf('COPY finished successully %s\n',self.current_time)
                    end
                else
                    fprintf('Either source or destination is empty\n');
                    return
                end
            else
                fprintf('DicomDownload: To use DicomDownload method you have to use one of the institute''s linux boxes\n');
            end
        end
        function              DicomTo4D(self,destination)
            %A wrapper over conversion and merging functions.
            
            %start with conversion
            self.ConvertDicom(destination);
            %finally merge 3D stuff to 4D and rename it data.nii.
            self.MergeTo4D(destination);
        end
        function              MergeTo4D(self,destination)
            %will create data.nii consisting of all the [f,s]TRIO images
            %merged to 4D. the final name will be called data.nii.
            % merge to 4D
            files       = spm_select('FPListRec',destination,'^fTRIO');
            fprintf('MergeTo4D:\nMerging (%s):\n',self.current_time);
            matlabbatch{1}.spm.util.cat.vols  = cellstr(files);
            matlabbatch{1}.spm.util.cat.name  = 'data.nii';
            matlabbatch{1}.spm.util.cat.dtype = 0;
            self.RunSPMJob(matlabbatch);            
            fprintf('Finished... (%s)\n',self.current_time);
            fprintf('Deleting 3D images in (%s)\n%s\n',self.current_time,destination);
            files = cellstr(files);
            delete(files{:});
        end
        function              ConvertDicom(self,destination)
            % dicom conversion. ATTENTION: dicoms will be converted and
            % then deleted. Assumes all files that start with MR are
            % dicoms.
            %
            
            matlabbatch = [];
            files       = spm_select('FPListRec',destination,'^MR');            
            if ~isempty(files)
                fprintf('ConvertDicom:\nFound %i files...\n',size(files,1));
                
                matlabbatch{1}.spm.util.import.dicom.data             = cellstr(files);
                matlabbatch{1}.spm.util.import.dicom.root             = 'flat';
                matlabbatch{1}.spm.util.import.dicom.outdir           = {destination};
                matlabbatch{1}.spm.util.import.dicom.protfilter       = '.*';
                matlabbatch{1}.spm.util.import.dicom.convopts.format  = 'nii';
                matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
                               
                fprintf('Dicom conversion s#%i... (%s)\n',self.id,self.current_time);
                self.RunSPMJob(matlabbatch);
                fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
                % delete dicom files
                fprintf('Deleting DICOM images in (%s)\n%s\n',self.current_time,destination);
                files = cellstr(files);
                delete(files{:});
                fprintf('Finished... (%s)\n',self.current_time);                
            else
                fprintf('No dicom files found for %i\n',self.id);
            end
        end        
        function [result]   = dicomserver_request(self)
            %will make a normal dicom request. use this to see the state of
            %folders
            results = [];
            if ~ismac & ~ispc
                fprintf('Making a dicom query, sometimes this might take long (so be patient)...(%s)\n',self.current_time);
                [status result]  = system(['/common/apps/bin/dicq --series --exam=' self.trio_session]);
                fprintf('This is what I found for you:\n');
            else
                fprintf('To use dicom query you have to use one of the institute''s linux boxes\n');
            end
            
        end
        function [paths]    = dicomserver_paths(self)
            paths = [];
            if ~ismac & ~ispc
                fprintf('Making a dicom query, sometimes this might take long (so be patient)...(%s)\n',self.current_time)
                [status paths]   = system(['/common/apps/bin/dicq -t --series --exam=' self.trio_session]);
                paths            = strsplit(paths,'\n');%split paths
            else
                fprintf('To use dicom query you have to use one of the institute''s linux boxes\n');
            end
        end
        function [data_path]= pathfinder(self,subject,run)
            %gets the path
            %Example: s.pathfinder(s.id,[]) => will return the path to
            %subject
            % empty [] above can be replaced with any phase number.
            data_path = self.path_project;
            for no = [subject run]                
                file_list        = dir(data_path);%next round of the forloop will update this
                i                = regexp({file_list(:).name},sprintf('[0,a-Z]%d$',no));%find all folders which starts with
                i                = find(cellfun(@(bla) ~isempty(bla), i  ));
                if ~isempty(i)
                    folder       = getfield(file_list(i),'name');
                    data_path    = sprintf('%s%s%s',data_path,filesep,folder);
                else
                    data_path    = [];
                    return
                end                
            end
            data_path(end+1)         = filesep;
        end
        function [out]      = dartel_templates(self,n)
            %returns the path to Nth Dartel template                                    
            out = fullfile(self.path_spm,'toolbox','cat12','templates_1.50mm',sprintf('Template_%i_IXI555_MNI152.nii',n) );
        end
        function              CreateFolderHierarchy(self)
            %Creates a folder hiearchy for a project. You must run this
            %first to create an hiearchy and fill this with data.
            for ns = 1:length(self.trio_sessions)
                for nr = 0:length(self.dicom2run)
                    for nf = 1:length(self.data_folders)                        
                        path2subject = sprintf('%s%ssub%03d%srun%03d%s%s',self.path_project,filesep,ns,filesep,nr,filesep,self.data_folders{nf});
                        if ~isempty(self.trio_sessions{ns})
                            a = fullfile(path2subject);
                            mkdir(a);
                        end
                    end
                end
            end
        end
        function t          = get.current_time(self)
            t = datestr(now,'hh:mm:ss');
        end
        function files      = collect_files(self,selector)
            %will collect all the files from all subjects that matches the
            %SELECTOR. Example selectors could be:
            %'run001/spm/model_03_chrf_0_0/s_w_beta_0001.nii'; 
            %
            % The selector string is simply appended to subjects' path. And
            % its existence is checked, in case of inexistenz method
            % stops.
            
            files = [];
            for ns = self.subject_indices                
                current = sprintf('%s/sub%03d/%s',self.path_project,ns,selector);
                if exist(current) ~= 0
                    files = [files ; current];
                else
                    cprintf([0 0 1],'Invalid Selector\n');
                    keyboard;%for sanity checks
                end
            end
        end
        function out        = path_atlas(self,varargin)
            %path to subjects native atlas, use VARARGIN to slice out a
            %given 3D volume.
            out = sprintf('%satlas/data.nii',self.path_project);
            if nargin > 1
                out = sprintf('%s,%d',out,varargin{1});
            end
        end
    end
    methods(Static) %SPM analysis related methods.       
        function RunSPMJob(matlabbatch)
            %will run the spm matlabbatch using the parallel toolbox.
            fprintf('Will call spm_jobman...\n');
           
            for n = 1:length(matlabbatch)
                fprintf('Running SPM jobman %i...\n',n);
                spm_jobman('run', matlabbatch(n));
            end            
        end                        
    end
    methods %getters
        function out        = getgroup_pmf_param(self)
            %collects the pmf parameters for all subjects in a table.
            out = [];
            for ns = self.subject_indices
                s       = Subject(ns);
                if ~isempty(s.pmf)
                    out = [out;[s.id s.pmf_param(:)']];
                end
            end
            out = array2table(out,'variablenames',{'subject_id' 'csp_alpha' 'csn_alpha' 'pooled_alpha' 'csp_beta' 'csn_beta' 'pooled_beta' 'csp_gamma' 'csn_gamma' 'pooled_gamma' 'csp_lamda'    'csn_lamda'    'pooled_lamda' });
        end
        function out        = getgroup_rating_param(self)
            %collects the rating parameters for all subjects in a table.
            out = [];
            for ns = self.subject_indices
                s   = Subject(ns);
                if ~isempty(s.pmf)%if there is no pmf data
                    out = [out;[s.id s.rating_param]];
                end
            end            
            out = array2table(out,'variablenames',{'subject_id' 'amp' 'fwhm' 'mu' 'offset' 'std' });
        end        
        function out        = path_beta_group(self,nrun,model_num,prefix,varargin)
            %returns the path for beta images for the second level
            %analysis. it simply loads single subjects' beta images and
            %replaces the string sub0000 with second_level.           
            dummy     = self.path_beta(nrun,model_num,prefix,varargin{:});
            for n = 1:size(dummy,1)
                out(n,:) = strrep(dummy(n,:),sprintf('sub%03d',self.id),'second_level');
            end
        end        
        function out        = path_contrast_group(self,nrun,model_num,prefix,type)
            %returns the path for F,T images for the second level
            %analysis. it simply loads single subjects' images and
            %replaces the string sub0000 with second_level.
            
            dummy     = self.path_contrast(nrun,model_num,prefix,type);
            for n = 1:size(dummy,1)
                out(n,:) = strrep(dummy(n,:),sprintf('sub%03d',self.id),'second_level');
            end
        end        
    end
    methods %methods that does something on all subjects one by one
        function VolumeGroupAverage(self,selector)
            %Averages all IMAGES across all subjects. This can only be done
            %on normalized images. The result will be saved to the
            %project/midlevel/. The string SELECTOR is appended to the
            %folder RUN.             
            %
            %For example to take avarage skullstripped image use: RUN = 0,
            %SELECTOR = 'mrt/w_ss_data.nii' (which is the normalized skullstripped
            %image). This method will go through all subjects and compute an average skull stripped image.
                        
            files       = self.collect_files(selector);            
            v           = spm_vol(files);%volume handles
            V           = mean(spm_read_vols(v),4);%data
            target_path = regexprep(files(1,:),'sub[0-9][0-9][0-9]','midlevel');%target filename
            dummy       = v(1);
            dummy.fname = target_path;
            mkdir(fileparts(target_path));
            spm_write_vol(dummy,V);
        end   
        function VolumeSmooth(self,files)
            %will smooth the files by the factor defined as Project
            %property.
            matlabbatch{1}.spm.spatial.smooth.data   = cellstr(files);
            matlabbatch{1}.spm.spatial.smooth.fwhm   = repmat(self.smoothing_factor,[1 3]);
            matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
            matlabbatch{1}.spm.spatial.smooth.im     = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's_';            
            spm_jobman('run', matlabbatch);
        end
        function SecondLevel_ANOVA(self,run,model,beta_image_index)
            % This method runs a second level analysis for a model defined in MODEL using beta images indexed in BETA_IMAGE_INDEX.
            
            %store all the beta_images in a 3D array            
            beta_files = [];
            for ns = self.subject_indices
                s          = Subject(ns);
%                 beta_files = cat(3,beta_files,s.path_contrast(run,model,'s_w_','T',beta_image_index)');%2nd level only makes sense with smoothened and normalized images, thus prefix s_w_
                  beta_files = cat(3,beta_files,s.path_beta(run,model,'s_w_',beta_image_index)');%2nd level only makes sense with smoothened and normalized images, thus prefix s_w_
            end
            %            
            c = 0;
            for ind = 1:length(beta_image_index)
                %take single beta_images across all subjects and store them
                %in a cell
                c = c +1;                
                files                                                              = squeeze(beta_files(:,ind,:))';
                matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(c).scans = cellstr(files);
            end
            %
            spm_dir                                                          = regexprep(self.dir_spmmat(run,model),'sub...','second_level');%convert to second-level path
            matlabbatch{1}.spm.stats.factorial_design.dir                    = cellstr(spm_dir);
            matlabbatch{1}.spm.stats.factorial_design.des.anova.dept         = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.variance     = 1;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca        = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova       = 0;
            matlabbatch{1}.spm.stats.factorial_design.cov                    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov              = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;
            %
            matlabbatch{2}.spm.stats.fmri_est.spmmat                         = {[spm_dir '/SPM.mat']};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical               =  1;
            %
            spm_jobman('run', matlabbatch);
        end
        function B = SecondLevel_Mumfordian(self,nrun,mask_id)
            %
            B = [];
            for ns = self.subject_indices
                s        = Subject(ns);
                beta     = s.analysis_mumfordian(nrun,mask_id);
                B        = cat(3,B,mean(beta,3));
            end
        end
    end
    methods %plotters
        function plot_ss_ratings(self)
            %will plot all single subject ratings in a big subplot
            tsubject = length(self.subject_indices);
            [y x]    = GetSubplotNumber(tsubject);
            c =0;
            for ns = self.subject_indices                
                c        = c+1;
                s        = Subject(ns);
                subplot(y,x,c)                
                s.plot_rating;
            end
            supertitle(self.path_project,1)
        end
        function plot_ss_pmf(self,chains)
            %will plot all single subject ratings in a big subplot
            tsubject = length(self.subject_indices);
            [y x]    = GetSubplotNumber(tsubject);
            c =0;
            for ns = self.subject_indices                
                c        = c+1;
                s        = Subject(ns);
                subplot(y,x,c)                                
                s.plot_pmf(chains);                
            end
            supertitle(self.path_project,1)
        end
        function plot_activitymap_normalized(self,files)
            %plots the data in a volume as average intensity projection, as
            %a 3xN panel, where N is the number of volumes. Masking
            %operates in the normalized space.
            
            %%                        
            %%                    
            files      = vertcat(files{:});
            vh         = spm_vol(files);%volume handle
            data_mat   = spm_read_vols(vh);%read the activity data, bonus: it checks also for orientation.
            [XYZmm]    = self.get_nativeatlas2mask(120);%world space;
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
            %%
            fields_alpha = {'x_alpha' 'y_alpha' 'z_alpha' };

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
        function plot_volume(self,filename)
            %will plot the volume using spm_image;
            global st                        
            spm_image('init',filename)
            spm_clf;
            st.callback = @self.test
            spm_orthviews('Image',filename)
%             spm_image('display',filename)
            spm_orthviews('AddColourBar',h(1),1);
            spm_orthviews('AddContext',h(1));
%             keyboard
        end
        function test(self)
            global st;
            coor = (st.centre);            
            self.default_model_name                                 = 'fir_15_15';
            %betas = self.path_beta(1,2,'s_',1:135);
            betas  = FilterF(sprintf('%s/midlevel/run001/spm/model_02_fir_15_15/s_w_beta_*',self.path_project));
            betas  = vertcat(betas{:});
            XYZvox = self.get_mm2vox([coor ;1],spm_vol(betas(1,:)))
            d = spm_get_data(betas,XYZvox);
            d = reshape(d,15,9);            
            figure(10001);   
            set(gcf,'position',[1370 1201 416  1104])
            set(gca, 'ColorOrder', GetFearGenColors, 'NextPlot', 'replacechildren');
            subplot(3,1,1)
            plot(d);box off;
            subplot(3,1,2)            
            imagesc(d);colormap jet
            subplot(3,1,3)                        
            set(gca, 'ColorOrder', [linspace(0,1,15);zeros(1,15);1-linspace(0,1,15)]', 'NextPlot', 'replacechildren');            
            data = [];
            data.x   = repmat(linspace(-135,180,8),15,1);
            data.y   = d(:,1:8);
            data.ids = 1:15;
            Y = repmat([1:15]',1,8);            
            h=bar(mean(d));            
            grid on;
            self.default_model_name                                 = 'chrf_0_0';
            %
            betas  = FilterF(sprintf('%s/midlevel/run001/spm/model_03_chrf_0_0/s_w_beta_*',self.path_project));
            betas  = vertcat(betas{:});
            XYZvox = self.get_mm2vox([coor(:); 1],spm_vol(betas(1,:)));
            d      = spm_get_data(betas,XYZvox);
            Weights2Dynamics(d(:)');          
        end
        
        
        function plot_overlay(self,file1,file2)
            spm_image('init',file1)
            spm_clf;
            pause(0.5);                                    
            v1 = spm_vol(file1);
            v2 = spm_vol(file2);        
            
            global st                                                
            st.callback = @self.test
            
            
            h  = spm_orthviews('Image' ,v1 );            
            spm_orthviews('Addtruecolourimage', h(1), file2, jet, 0.8 );
            spm_orthviews('Redraw');
            spm_orthviews('AddColourBar',h(1),1);
            spm_orthviews('AddContext',h(1));
            spm_orthviews('Caption', h(1),file2);
        end
        function hist_pmf_param(self,fields)
            %will generate an histogram of the pmf parameters
            out = self.getgroup_pmf_param;
            counter = 0;
            for fields = fields
                counter = counter+1;
                x   = linspace(0,120,11);
                c   = histc(out.(fields{1}),x);
                subplot(1,length(fields),counter)
                bar(x,c,1,'k');
                axis tight;
                box off;
                xlabel('alpha');                
                M = mean(out.(fields{1}));
                S = std(out.(fields{1}))./sqrt(size(M,1));
                hold on;
                plot([M M],ylim,'r','linewidth',4);
                hold off;
                title(sprintf('%s: %3.4g ± %3.4g (M±SEM)',fields{1},M,S),'interpreter','none');
            end
        end
        function hist_rating_param(self)
            %%will generate an histogram of the rating parameters.
            out     = self.getgroup_rating_param;
            counter = 0;            
            for fields = {'amp' 'kappa' 'mu'}
                counter = counter+1;                
                Y   = out.(fields{1});
                x   = linspace(min(Y),max(Y),20);
                c   = histc(Y,x);                
                subplot(1,3,counter)
                bar(x,c,1,'k');
                axis tight;
                box off;
                xlabel('alpha');                
                M = mean(out.(fields{1}));
                S = std(out.(fields{1}))./sqrt(size(M,1));
                hold on;
                plot([M M],ylim,'r');
                hold off;
                title(sprintf('%s: %3.4g ± %3.4g (M±SEM)',fields{1},M,S),'interpreter','none');
            end
        end
    end
end
