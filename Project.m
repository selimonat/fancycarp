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
        smoothing_factor      = 1:10;%how many mm images should be smoothened when calling the SmoothVolume method        
        path_smr              = sprintf('%s%ssmrReader%s',fileparts(which('Project')),filesep,filesep);%path to .SMR importing files in the fancycarp toolbox.
        gs                    = [5 6 7 8 9 12 14 16 17 18 19 20 21 24 25 27 28 30 32 34 35 37 39 40 41 42 43 44];
        subjects              = [];
        bs                    = [10 11 13 15 22 23 26 29 31 33 36 38];        
        %
        mbi_ucs               = [5 13 16 20 22 28 33 36 40 45 47 53 54 61];
        mbi_oddball           = [6 65];
        mbi_transition        = [22 23 45 46];
        mbi_invalid           = sort([Project.mbi_ucs,Project.mbi_oddball,Project.mbi_transition]);
        mbi_valid             = setdiff(1:65,Project.mbi_invalid);
        roi                   = struct('xyz',{[30 0 -24] [30 22 -8] [33 22 6] [39 22 -6] [44 18 -13] [-30 18 -6] [-9 4 4] [40 -62 -12]},'label',{'rAmy' 'rAntInsula' 'rAntInsula2' 'rAntInsula3' 'rFronOper' 'lAntInsula' 'BNST' 'IOC'} );%in mm
    end
    properties (Constant,Hidden) %These properties drive from the above, do not directly change them.
        tpm_dir               = sprintf('%stpm/',Project.path_spm); %path to the TPM images, needed by segment.         
        path_second_level     = sprintf('%sspm/',Project.path_project);%where the second level results has to be stored        
		current_time          = datestr(now,'hh:mm:ss');
        subject_indices       = find(cellfun(@(x) ~isempty(x),Project.trio_sessions));% will return the index for valid subjects (i.e. where TRIO_SESSIONS is not empty). Useful to setup loop to run across subjects.
        PixelPerDegree        = 29;
        screen_resolution     = [768 1024];
        path_stim             = sprintf('%sstim/data.png',Project.path_project); %path to the TPM images, needed by segment.         
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
        function XYZmm      = get_XYZmmNormalized(self,mask_id)
            %Will return XYZ coordinates from ROI specified by MASK_INDEX
            %thresholded by the default value. XYZ values are in world 
            %space, so they will need to be brought to the voxel space of
            %the EPIs.            
            mask_handle = spm_vol(self.path_atlas(mask_id));%read the mask
            mask_ind    = spm_read_vols(mask_handle) > self.atlas2mask_threshold;%threshold it            
            [X Y Z]     = ind2sub(mask_handle.dim,find(mask_ind));%get voxel indices
            XYZ         = [X Y Z ones(sum(mask_ind(:)),1)]';%this is in hr's voxels space.
            XYZmm       = mask_handle.mat*XYZ;%this is world space.            
            XYZmm       = unique(XYZmm','rows')';%remove repetititons.
        end
        function D          = getgroup_data(self,file,mask_id)
            %will read the data specified in FILE 
            %FILE is the absolute path to a 3/4D .nii file. Important point
            %is that the file should be in the normalized space, as MASK_ID
            %refers to a normalized atlas.
            %
            %MASK_ID is used to select voxels.
            %                                    
            vh      = spm_vol(file);
            if spm_check_orientations(vh)
                XYZmm   = self.get_XYZmmNormalized(mask_id);
                XYZvox  = self.get_mm2vox(XYZmm,vh(1));%in EPI voxel space.
                D       = spm_get_data(vh,XYZvox);
            else
                fprintf('The data in\n %s\n doesn''t have same orientations...',file);
            end            
        end
        function sub        = get_selected_subjects(self,criteria,inversion)
            %will select subjects based on different CRITERIA. Use
            %INVERSION to invert the selection, i.e. use the unselected
            %subjects. LIST is the indices of subjects, NAME is the name of
            %the group, used to write folder names, identifying a given
            %group of subjects. So dont get confused, when INVERSION is 0,
            %you get the selected subjects (generaly performing better).
            
            if nargin == 2
                inversion = 1;%will select good subjects            
            end
            
            if criteria == 0
                sub.name = '00_all';
                sub.list = self.subject_indices;
            
            elseif     criteria == 1
                %
                sub.name   = '01_rating';
                %
                out        = self.getgroup_rating_param;%load the ratings of all subjects
                select     = (out.rating_LL > -log10(.05))&(out.rating_mu > -90)&(out.rating_mu < 90);%select based on tuning
                select     = select == inversion;
                sub.list   = out.subject_id(select);%subject indices.
                %               
            elseif criteria == 2
                %
                sub.name = '02_detection'; 
                out      = self.getgroup_detected_face;
                select   = abs(out.detection_face) <= 45;
                select   = select==inversion;
                sub.list = out.subject_id(select);
                
            elseif criteria == 3
                
                %same as rating tuning, but with saliency data.
                sub.name = '03_saliency';
                out      = self.getgroup_facecircle;                
                select   = (out.facecircle_LL > -log10(.05))&(out.facecircle_mu > -90)&(out.facecircle_mu < 90);%select based on tuning                
                select   = select==inversion;
                sub.list = out.subject_id(select);
                
            elseif criteria == 4
                
                sub.name    = '04_perception';                
                out     = self.getgroup_pmf_param;%all threshoild values
                m       = nanmedian(out.pmf_pooled_alpha);
                select  = out.pmf_pooled_alpha <= m;%subjects with sharp pmf
                select  = select==inversion;
                sub.list= out.subject_id(select);
                
            end
            
            %invert the group name necessary.
            if ~inversion
                sub.name = [sub.name '_minus'];
            end
        end
% % %         function out        = path_spmmat_SecondLevel_ANOVA(self,run,model,sk)
% % %             %returns the path to second level analysis for a given
% % %             %smoothening kernel.
% % %             
% % %             out = regexprep(self.path_spmmat(run,model),'sub...',sprintf('second_level_%02dmm',sk));
% % %             
% % %         end
    end    
    methods %getters
        function out        = getgroup_all_param(self)
            %
            out = self.getgroup_pmf_param;
            out = join(out,self.getgroup_rating_param);
            out = join(out,self.getgroup_facecircle);
            out = join(out,self.getgroup_detected_oddballs);
            out = join(out,self.getgroup_detected_face);                        
            
        end
        function out        = getgroup_pmf_param(self)
            %collects the pmf parameters for all subjects in a table.
            out = [];
            for ns = self.subject_indices
                s       = Subject(ns);
                if ~isempty(s.pmf)
                    out = [out;[s.id s.pmf_param(:)' s.fit_pmf.LL']];
                else
                    out = [out;[s.id NaN(1,size(out,2)-1)]];
                end
            end
            out = array2table(out,'variablenames',{'subject_id' 'pmf_csp_alpha' 'pmf_csn_alpha' 'pmf_pooled_alpha' 'pmf_csp_beta' 'pmf_csn_beta' 'pmf_pooled_beta' 'pmf_csp_gamma' 'pmf_csn_gamma' 'pmf_pooled_gamma' 'pmf_csp_lamda'    'pmf_csn_lamda'    'pmf_pooled_lamda' 'pmf_csp_LL' 'pmf_csn_LL' 'pmf_pooled_LL' });
        end
        function out        = getgroup_rating_param(self)
            %collects the rating parameters for all subjects in a table.
            out = [];
            for ns = self.subject_indices%run across all the subjects
                s   = Subject(ns);
                out = [out;[s.id s.rating_param s.fit_rating.LL]];                
            end            
            out = array2table(out,'variablenames',{'subject_id' 'rating_amp' 'rating_fwhm' 'rating_mu' 'rating_offset' 'rating_std' 'rating_LL'});
        end       
        function out        = getgroup_facecircle(self)
            %collects the rating parameters for all subjects in a table.
            out = [];
            for ns = self.subject_indices
                s   = Subject(ns);                
                out = [out;[s.id s.fit_facecircle(1,8).params s.fit_facecircle(1,8).LL]];                
            end            
            out = array2table(out,'variablenames',{'subject_id' 'facecircle_amp' 'facecircle_fwhm' 'facecircle_mu' 'facecircle_offset' 'facecircle_std' 'facecircle_LL'});
        end       
        function out        = getgroup_detected_oddballs(self)            
            odd = [];
            for ns = self.subject_indices
                s       = Subject(ns);
                odd     = [odd ;[s.id sum(s.detected_oddballs)]];
            end
            out = array2table(odd,'variablenames',{'subject_id' 'detection_oddball'});
        end
        function out        = getgroup_detected_face(self)            
            face = [];
            for ns = self.subject_indices
                s      = Subject(ns);
                face   = [face ;[s.id s.detected_face]];                
            end
            out = array2table(face,'variablenames',{'subject_id' 'detection_face'});
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
        function out        = getgroup_all_spacetime(self,who,sk)            
            
            out = [];
            for ns = self.get_selected_subjects(who)'
                if ns ~= 28
                    fprintf('Collecting subject %03d\n',ns);
                    s     = Subject(ns,'pupil');
                    dummy = s.get_all_spacetime(sk);
                    out   = cat(5,out,dummy);
                else
                    fprintf('28 excluded~~~\n')
                end
            end            
        end
        function out        = getgroup_pupil_spacetime(self,who)
            %set WHO to select subject groups, 1:good, 0:bad;
            out = [];
            for ns = self.get_selected_subjects(who)'
                fprintf('Collecting subject %03d\n',ns);
                s     = Subject(ns,'pupil');
                dummy = s.get_pupil_spacetime;
                out   = cat(4,out,dummy);
            end
        end
        function  name      = get_atlasROIname(self,roi_index)
            %get_atlasROIname(self,roi_index)
            %
            %returns the name of the probabilistic atlas at ROI_INDEX along
            %the 4th dimension of the atlas data.            
            cmd           = sprintf('sed "%dq;d" %s/data.txt',roi_index,fileparts(self.path_atlas));
            [~,name]      = system(cmd);
            name          = deblank(name);
        end
        function roi        = get_atlaslabel(self,XYZmm)
            %returns the label of the ROI at location XYZmm using the
            %Project's atlas.
                        
            atlas_labels = regexprep(self.path_atlas,'nii','txt');
            roi.name     = {};
            roi.prob     = {};
            %            
            V            = spm_vol(self.path_atlas);
            V            = V(63:end);%take only the bilateral ones.
            %get the probabilities at that point
            [p i]    = sort(-spm_get_data( V, V(1).mat\[XYZmm(:);1]));
            %black list some useless areas
            [~,a]      = system(sprintf('cat %s | grep Cerebral',atlas_labels));
            black_list = strsplit(a);
            %list the first 6 ROIs
%             fprintf('Cursor located at (mm): %g %g %g\n', x,y,z);
            counter = 0;
            while counter <= 3
                counter               = counter + 1;
                %get the current ROI name
                current_name = self.get_atlasROIname(i(counter)+62);%the first 62 are bilateral rois, not useful.          
                %include if the probability is bigger than zero and also if the current roi is not
                %one of the stupid ROI names, e.g. cerebral cortex or so...             
                if -p(counter) > 0 & isempty(cell2mat( regexp(current_name,black_list)))
                    %isempty == 1 when there is no match, which is good
                    %                    
                    roi.name     = [roi.name current_name];
                    roi.prob     = [roi.prob -p(counter)];
                    fprintf('ROI: %d percent chance for %s.' , roi.prob{end}, roi.name{end}(1:end-1));
                    fprintf('\n');                    
                end
            end
            %if at the end of the loop we are still poor, we say it
            if isempty(roi.prob)
                roi.name{1} = 'NotFound';
                roi.prob(1) = NaN;
            end
            
            fprintf('\n\n');
        end
    end
    methods %methods that does something on all subjects one by one
        function              VolumeGroupAverage(self,selector)
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
            V           = mean((spm_read_vols(v)).^2,4);%data
            target_path = regexprep(files(1,:),'sub[0-9][0-9][0-9]','midlevel');%target filename
            dummy       = v(1);
            dummy.fname = target_path;
            mkdir(fileparts(target_path));
            spm_write_vol(dummy,V);
        end   
        function              VolumeSmooth(self,files)
            %will smooth the files by the factor defined as Project
            %property.
            for sk = self.smoothing_factor
                matlabbatch{1}.spm.spatial.smooth.data   = cellstr(files);
                matlabbatch{1}.spm.spatial.smooth.fwhm   = repmat(sk,[1 3]);
                matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
                matlabbatch{1}.spm.spatial.smooth.im     = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('s%02d_',sk);
                spm_jobman('run', matlabbatch);
            end
        end        
        function [xSPM]     = SecondLevel_ANOVA(self,ngroup,run,model,beta_image_index,sk)
            % This method runs a second level analysis for a model defined
            % in MODEL using beta images indexed in BETA_IMAGE_INDEX. The
            % same model can be ran at different smoothening levels (SK) as well
            % as for different subject groups (NGROUP, see
            % self.get_selected_subjects).
            %
            % Results are saved at the project root. xSPM structure is also
            % saved in the analysis SPM dir. If this file is present, it is
            % directly loaded, so the second level analysis is not reran.
            %
            %
            
            
            
            subjects   = self.get_selected_subjects(ngroup);                        
            spm_dir    = regexprep(self.dir_spmmat(run,model),'sub...','second_level');%convert to second-level path, replace the sub... to second-level.
            spm_dir    = sprintf('%s%02dmm/group_%s/',spm_dir,sk,subjects.name);
            xspm_path  = sprintf('%sxSPM.mat',spm_dir);
            if exist(xspm_path) == 0;
                                
                beta_files = [];
                for ns = subjects.list(:)'
                    s          = Subject(ns);
                    %store all the beta_images in a 3D array
                    beta_files = cat(3,beta_files,s.path_beta(run,model,sprintf('s%02d_w_',sk),beta_image_index)');%2nd level only makes sense with smoothened and normalized images, thus prefix s_w_
                end
                %
                c = 0;
                for ind = 1:length(beta_image_index)
                    %take single beta_images across all subjects and store them
                    %in a cell
                    c                                                                  = c +1;
                    files                                                              = squeeze(beta_files(:,ind,:))';
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(c).scans = cellstr(files);
                end
                %depending on subject group and smoothening, generate
                %a directory name for this spm analysis.
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
                %UPDATE THE SPM.mat with the contrast
                %create the contrast that we are interested in, this could
                %be later extended into a for loop
                load(sprintf('%s/SPM.mat',spm_dir))
                SPM = rmfield(SPM,'xCon');%result the previous contrasts, there shouldnt by any
                tbetas = length(beta_image_index);
                if model == 7
                    SPM.xCon(1) = spm_FcUtil('set','eoi','F','c',[[0 0 0 ]' eye(3)]',SPM.xX.xKXs);
                elseif model == 2
                    SPM.xCon(1) = spm_FcUtil('set','eoi','F','c',[eye(8)]',SPM.xX.xKXs);
                    kappa       = .5;
                    SPM.xCon(2) = spm_FcUtil('set','tuning','F','c',[zscore(Tuning.VonMises( linspace(-135,180,8),1,kappa,0,0))]',SPM.xX.xKXs);
                else
                    keyboard
                end
                
                save(sprintf('%s/SPM.mat',spm_dir),'SPM');%save the SPM with the new xCon field
                %xSPM is used to threshold according to a contrast.
                xSPM = struct('swd', spm_dir,'title','eoi','Ic',1,'n',1,'Im',[],'pm',[],'Ex',[],'u',.00001,'k',0,'thresDesc','none');
                %replace 'none' to make FWE corrections.
                [SPM xSPM] = spm_getSPM(xSPM);%now get the tresholded data, this will fill in the xSPM struct with lots of information.
                save(sprintf('%s/SPM.mat',spm_dir),'SPM');%save the SPM with the new xCon field
                save(xspm_path,'xSPM');
            else
                load(xspm_path);                
%                 t                                 = spm_list('table',xSPM{2}.rating);
            end            

        end        
        function B          = SecondLevel_Mumfordian(self,nrun,mask_id)
            %
            B = [];
            for ns = self.subject_indices
                s        = Subject(ns);
                beta     = s.analysis_mumfordian(nrun,mask_id);
                B        = cat(3,B,mean(beta,3));
            end
        end
        function              LatencyMap(self,midlevel_path,xBF)
            %Will read all beta images in midlevel_path consisting of a FIR
            %model.
            beta_files      = strvcat(FilterF(sprintf('%smidlevel/',self.path_project),midlevel_path,'beta'));
            total_beta      = size(beta_files,1);
            target_path     = strvcat(regexp(fileparts(beta_files(1,:)),'run0.*','match'));
            target_path     = sprintf('%smidlevel/%s/',self.path_project,target_path);
            target_path     = regexprep(target_path,'spm','latencymaps');
            if exist(target_path) == 0
                mkdir(target_path);
            end                
            fprintf('Found %i beta files here...\n',total_beta);
            if total_beta == 0
                fprintf('Found 0 beta files here...\n');
                keyboard
            end
            %% read the data [x,y,z,time,cond], and put it into a [x*y*z,time,cond];
            time_points     = size(xBF.bf,2);
            total_condition = total_beta/time_points;
            D               = spm_read_vols(spm_vol(beta_files));%[x,y,z,time,cond]
            ori_size        = size(D);
            D               = reshape(D,size(D,1)*size(D,2)*size(D,3),total_beta/total_condition,total_condition);%[x*y*z,time,cond]
            D               = mean(D,3)';%average across conditions for more reliable estimation, this is the first step in this line.
            cprintf([1 0 0],'ATTENTION ONLY FIRST 13 values are taken...\n');
            D               = xBF.bf(:,1:13)*D(1:13,:);
            %%
            Time            = linspace(0,xBF.order*self.TR,xBF.order*xBF.T+1)-5;
            pre_stim        = Time <= 0;
            pre_std         = std(D(pre_stim,:));
            post_std        = std(D(~pre_stim,:));
            post_amp        = mean(D(~pre_stim,:));
            for p = linspace(10,90,10)
                threshold  = prctile(post_amp,p);
                for factor = linspace(1,5,10)
                    valid           = post_std > (pre_std*factor);%valid voxels with high SNR
                    valid           = valid & (post_amp > threshold);                    
                    [~,Latency]     = max(D);%sample with peak value.
                    Latency         = Time(Latency);
                    valid           = valid & (Latency<7);
                    Latency(~valid) = 0;
                    %cancel all voxels with low SNR
                    Latency         = reshape(Latency,ori_size(1),ori_size(2),ori_size(3));
                    %%
                    dummy           = spm_vol(beta_files(1,:));
                    dummy.fname     = sprintf('%s/latency_map_factor_%1.2g_%1.2g.nii',target_path,factor,p);
                    spm_write_vol(dummy,Latency);
                end
            end
        end
        function nface      = analysis_selected_face(self,varargin)
            %will plot an histogram of detected faces 
            subjects = self.subject_indices
            if nargin > 1
                subjects = varargin{1};
            end
            nface = [];            
            for ns = subjects;
                s     = Subject(ns);
                nface = [nface s.detected_face];                
            end  
            self.plot_newfigure;
            self.plot_bar(histc(nface,linspace(-135,180,8)));
        end
        function [count_facewc,count_facew]=analysis_facecircle(self)
            % will conduct a group level face saliency analysis, will also
            % save the weights associated to each position. The analysis is
            % limited to the first 60 fixations, this seems to be the
            % number of fixations where everybody has equal numbers i.e.
            % number of participants with 70, 80, 90 fixations drops
            % linearly with number of fixations.
            %
            viz = 1;
            raw   = [];%will contain fix by fix all the info
            for ns = self.subject_indices;
                s         = Subject(ns);
                [dummy]   = s.get_facecircle(1);%with no partition.
                raw       = [raw dummy.raw];
            end
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
            %                            9/Face Index, this is the filename of the stimulus
            %                            10/Rank; this is the order of the fixation along the 30 seconds
            %                            11/START time
            %                            12/angle of fixation maps aligned to CS+.
            %                            13/old amp.
            %                            14/subject index.
            counts       = histc( raw(10,:) , 1:100);
            limit        = 60;
            fprintf('I will take the first %i fixations\n',limit);
            %delete all fixations above LIMIT
            valid        = raw(10,:) <= limit;
            raw          = raw(:,valid);
            PositionBias = accumarray(raw(6,:)',1);
            W            = 1./PositionBias;
            save(sprintf('%smidlevel/run003/facecircle_position_weights_limit_%i.mat',s.path_project,limit),'W');
            %%
            if viz
                figure(1);set(gcf,'position',[67   350   327   752]);
                bar(circshift(PositionBias,[2 0]),.9,'k')
                set(gca,'color','none','xticklabel',{'' '' '' '12 Uhr' '' '' '' '6 Uhr'});xlabel('positions');box off;ylabel('# fixations');xlim([.5 8.5]);ylim([0 700]);SetTickNumber(gca,3,'y');                            
                % sanity check
                figure(2);set(gcf,'position',[67   350   327   752]);
                bar(accumarray(raw(6,:)',W(raw(6,:)')),.9,'k');
                set(gca,'color','none','xticklabel',{'' '' '' '12 Uhr' '' '' '' '6 Uhr'},'ytick',1);xlabel('positions');box off;ylabel('# wfixations');xlim([.5 8.5]);ylim([0 1.5]);
                % overall saliency profile            
                figure(3);set(gcf,'position',[67   350   327   752]);
                h=bar(accumarray(raw(9,:)',W(raw(6,:)')),.9);
                SetFearGenBarColors(h);
                set(gca,'color','none','xticklabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'ytick',1,'ygrid','on');xlabel('positions');box off;ylabel('# wfixations');xlim([.5 8.5]);ylim([0 1.5]);
                set(get(h,'Children'),'FaceAlpha',1)
                % overall saliency without correction
                figure(4);set(gcf,'position',[67   350   327   752]);
                h=bar(accumarray(raw(9,:)',1),.9);
                SetFearGenBarColors(h)                
                set(gca,'color','none','xticklabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'ygrid','on');xlabel('positions');box off;ylabel('# wfixations');xlim([.5 8.5]);
                set(get(h,'Children'),'FaceAlpha',1)
            end
            %%  sort fixations rank by rank
            count_facew=[];
            count_posw =[];
            for nrank = 1:limit;
                valid                = raw(10,:) == nrank;                       
                count_facew(nrank,:)  = accumarray(raw(9,valid)',W(raw(6,valid)'),[8 1]);
                count_posw(nrank,:)   = accumarray(raw(6,valid)',W(raw(6,valid)'),[8 1]);                
            end
            %%
            if viz
                figure(4);set(gcf,'position',[680         681        1001         417]);
                subplot(1,4,1);
                imagesc(1:8,1:60,count_facew);thincolorbar;axis xy;ylabel('fixation order');xlabel('condition');
                set(gca,'xticklabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'xtick',1:8);                
                subplot(1,4,2);
                count_facewc = conv2([count_facew count_facew count_facew],ones(2,2),'same');
                count_facewc = count_facewc(:,9:16);
                imagesc(count_facewc);axis xy;set(gca,'xticklabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'xtick',1:8,'yticklabel',{});                
                subplot(1,4,3);
                imagesc(1:8,1:60,cumsum(count_facewc));axis xy;set(gca,'yticklabel',{''});set(gca,'xticklabel',{'CS+'  'CS-'},'xtick',[4 8],'xgrid','on');
                subplot(1,4,4);
                contour(1:8,1:60,cumsum(count_facewc),30);axis xy;set(gca,'yticklabel',{''});set(gca,'xticklabel',{'CS+' 'CS-'},'xtick',[4 8],'xgrid','on');
            end
            %% model comparison
            data.y           = cumsum(count_facewc);
            data.x           = repmat(linspace(-135,180,8),size(data.y,1),1);
            data.ids         = 1:size(data.y,1);
            t2               = Tuning(data);
            t2.visualization = 0;
            t2.gridsize      = 10;
            t2.SingleSubjectFit(8);
            t2.SingleSubjectFit(7);
            winfun = [8 7];
            if viz
                figure(5);clf;
                for nfix = 1:limit
                    hold on;
                    X = t2.fit_results{7}.x_HD(1,:);
                    
                    if (t2.fit_results{8}.pval(nfix) > t2.fit_results{7}.pval(nfix))%Gaussian wins
                        winfun = 8
                        Est    = t2.fit_results{winfun}.params(nfix,:);
                        Est(4) = 0;
                        c = 'r';
                    else% cosine wins
                        winfun = 7
                        Est    = t2.fit_results{winfun}.params(nfix,:);
                        Est(3) = Est(1);
                        c = 'b';
                    end
                    Y      = t2.fit_results{winfun}.fitfun(X,Est);
                    plot3(X,repmat(nfix,1,size(X,2)),Y,'color',c,'linewidth',1);
                end
                hold off;
                view(43,32)
                ylim([0 60])
                DrawPlane(gcf,.4);
                
                set(gca,'xtick',[0 180],'xticklabel',{'high' 'low'},'xgrid','on','ytick',[0 30 60],'ztick',[0 .35 .7],'ygrid','on','color','none')
                xlabel('similarity');
                ylabel('fixation');
                figure(6)
                plot(t2.fit_results{7}.pval,'linewidth',3);
                hold on;
                plot(t2.fit_results{8}.pval,'r','linewidth',3);
                box off;
                xlabel('fixations');
                ylabel('LRT (-log(p))');
                ylim([0 8])
                legend({'Cosine vs. Null','Gaussian vs. Null'});legend boxoff;
                set(gca,'color','none')
            end
        end        
        
    end
    
    methods %plotters
        function plot_normalized_surface(self);
            cat_surf_display(struct('data',{{'/Volumes/feargen2/project_fearamy/data/midlevel/run000/mrt/surf/lh.central.w_ss_data.gii' '/Volumes/feargen2/project_fearamy/data/midlevel/run000/mrt/surf/rh.central.w_ss_data.gii'}},'multisurf',0))
        end
        function plot_ss_ratings(self)
            figure;set(gcf,'position',[5           1        1352        1104]);
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
            figure;set(gcf,'position',[5           1        1352        1104]);
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
        function plot_ss_facecircle(self,partition,fun)
            figure;set(gcf,'position',[5           1        1352        1104]);
            %will plot all single subject ratings in a big subplot
            tsubject = length(self.subject_indices);
            [y x]    = GetSubplotNumber(tsubject);
            c =0;
            for ns = self.subject_indices                
                c        = c+1;
                s        = Subject(ns);
                subplot(y,x,c)                
                s.plot_facecircle(partition,fun);
            end
            supertitle(self.path_project,1)
        end
        function plot_ss_detected_oddballs(self)
            o = self.getgroup_detected_oddballs;
            o = o.detection_oddball;
            bar(1:40,o,.9,'k');
            xlabel('subjects');
            ylabel('# detected oddballs out of 2');
            box off;
            ylim([0 3]);
            xlim([0 41]);
            set(gca,'ytick',[0 1 2],'xtick',1:40,'xticklabel',5:44);
            title('oddball performace');
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
%             st.callback = @self.test
            spm_orthviews('Image',filename)
%             spm_image('display',filename)
            spm_orthviews('AddColourBar',h(1),1);
            spm_orthviews('AddContext',h(1));
%             keyboard
        end        
        function plot_newfigure(self)
            figure;
            set(gcf,'position',[671   462   506   477]);
        end
        function test3(self)
            
            global st;
            try               
                coor   = (st.centre);                                
                betas  = self.path_beta(1,6,'s_w_',1:488);                
                tbeta  = size(betas,1);
                XYZvox = self.get_mm2vox([coor ;1],spm_vol(betas(1,:)));
                d = spm_get_data(betas,XYZvox);
                d = reshape(d,tbeta/8,8);                
                d = conv2([d d d],ones(2,2),'same');
                d = d(:,9:16);
                figure(10005);
                set(gcf,'position',[919        1201         861        1104])                
                subplot(1,3,1)
                imagesc(d);colorbar;axis xy;
                subplot(1,3,2)
                imagesc(cumsum(d));colorbar
                axis xy;
                subplot(1,3,3);
                bar(mean(d));                
                
                
            catch
                fprintf('call back doesn''t run\n')
            end
        end        
        function cb_fir_vs_hrf(self)
            global st
            try               
                coor                   = (st.centre);
                figure(10002);
                set(gcf,'position',[950        1201         416        1104])                
                set(gca, 'ColorOrder', GetFearGenColors, 'NextPlot', 'replacechildren');
                %% get BF for FIR and HRF, this should be exactly the same
                %ones we used for the first-levels
                TR              = self.TR;
                xBF.T           = 16;
                xBF.T0          = 1;
                xBF.dt          = TR/xBF.T;
                xBF.UNITS       = 'scans';
                xBF.Volterra    = 1;
                xBF.name        = 'Finite Impulse Response';
                xBF.order       = 15;
                xBF.length      = 15*TR;
                fir_xBF         = spm_get_bf(xBF);
                %
                TR              = self.TR;
                xBF.T           = 16;
                xBF.T0          = 1;
                xBF.dt          = TR/xBF.T;
                xBF.UNITS       = 'scans';
                xBF.Volterra    = 1;
                xBF.name        = 'hrf';
                xBF.order       = 1;
                xBF.length      = 15*TR;
                hrf_xBF         = spm_get_bf(xBF);
                %
                %% now load the betas for both the FIR and HRF models.
                betas           = FilterF(sprintf('%s/midlevel/run001/spm/model_02_fir_15_15/s_w_beta_*',self.path_project));
                betas           = vertcat(betas{:});                
                XYZvox          = self.get_mm2vox([coor ;1],spm_vol(betas(1,:)));
                d               = spm_get_data(betas,XYZvox);
                d               = reshape(d,15,9);%this is the FIR model fit.
                subplot(2,1,1);
                set(gca, 'ColorOrder', GetFearGenColors, 'NextPlot', 'replacechildren');
                plot(d);
                %% now load the hrf model
                betas           = FilterF(sprintf('%s/midlevel/run001/spm/model_02_chrf_0_0/s_w_beta_*',self.path_project));
                betas           = vertcat(betas{:});
                XYZvox          = self.get_mm2vox([coor ;1],spm_vol(betas(1,:)));
                B               = spm_get_data(betas,XYZvox);%these are 9 betas images
                Fit             = hrf_xBF.bf(:,1:16:end,:)*B';
                subplot(2,1,2);
                set(gca, 'ColorOrder', GetFearGenColors, 'NextPlot', 'replacechildren');
                Time            = [1:518]*hrf_xBF.dt;                
                plot(Time,Fit);xlim([0 15]);

            end
        end        
        function test2(self)
            
            global st;
            try
                TR                     =.99;
                xBF.T                  = 16;
                xBF.T0                 = 1;
                xBF.dt                 = TR/xBF.T;
                xBF.UNITS              = 'scans';
                xBF.Volterra           = 1;
                xBF.name               = 'Fourier set';
                xBF.order              = 20;
                xBF.length             = 20*TR;
                xBF                    = spm_get_bf(xBF);
                %
                Time                   = linspace(0,20*TR,20*16+1)-5;
                %
                
                coor = (st.centre);
                self.default_model_name                                 = 'fourier_20_20';
                %betas = self.path_beta(1,2,'s_',1:135);
                betas  = FilterF(sprintf('%s/midlevel/run001/spm/model_02_fourier_20_20/s_w_beta_*',self.path_project));
                betas  = vertcat(betas{:});
                XYZvox = self.get_mm2vox([coor ;1],spm_vol(betas(1,:)))
                d = spm_get_data(betas,XYZvox);
                d = reshape(d,41,9);
                d(13:end,:) = 0;
                d = xBF.bf*d;                
                figure(10002);
                set(gcf,'position',[950        1201         416        1104])                
                subplot(3,1,1)
                set(gca, 'ColorOrder', GetFearGenColors, 'NextPlot', 'replacechildren');
                plot(Time,d);
                box off;
                hold on;
                plot(Time,mean(d(:,1:8),2),'k--','linewidth',3);                
                plot([0 0], ylim,'k');
                xlabel('time (s)');
                hold off;
                subplot(3,1,2)
                imagesc(d);colormap jet;
                subplot(3,1,3)                
                set(gca, 'ColorOrder', [linspace(0,1,15);zeros(1,15);1-linspace(0,1,15)]', 'NextPlot', 'replacechildren');
                data = [];
                data.x   = repmat(linspace(-135,180,8),15,1);
                data.y   = d(:,1:8);
                data.ids = 1:15;
                Y = repmat([1:15]',1,8);
                h=bar(mean(d));
%                 SetFearGenBarColors(h);
                grid on;
                self.default_model_name                                 = 'chrf_0_0';
                %
                betas  = FilterF(sprintf('%s/second_level/run001/spm/model_05_chrf_0_0/beta_*',self.path_project));
                betas  = vertcat(betas{:});
                XYZvox = self.get_mm2vox([coor(:); 1],spm_vol(betas(1,:)));
                d      = spm_get_data(betas,XYZvox);
                Weights2Dynamics(d(:)');
            catch
                fprintf('call back doesn''t run\n')
            end
        end
        function test(self)
            
            global st;
            try
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
                subplot(3,1,1)
                set(gca, 'ColorOrder', GetFearGenColors, 'NextPlot', 'replacechildren');
                plot(d);box off;
                hold on;plot(mean(d(:,1:8),2),'k--','linewidth',3);                
                hold off;
                subplot(3,1,2)
                imagesc(d);colormap jet;
                subplot(3,1,3)                
                set(gca, 'ColorOrder', [linspace(0,1,15);zeros(1,15);1-linspace(0,1,15)]', 'NextPlot', 'replacechildren');
                data = [];
                data.x   = repmat(linspace(-135,180,8),15,1);
                data.y   = d(:,1:8);
                data.ids = 1:15;
                Y = repmat([1:15]',1,8);
                h=bar(mean(d));
%                 SetFearGenBarColors(h);
                grid on;
                self.default_model_name                                 = 'chrf_0_0';
                %
                betas  = FilterF(sprintf('%s/second_level/run001/spm/model_05_chrf_0_0/beta_*',self.path_project));
                betas  = vertcat(betas{:});
                XYZvox = self.get_mm2vox([coor(:); 1],spm_vol(betas(1,:)));
                d      = spm_get_data(betas,XYZvox);
                Weights2Dynamics(d(:)');
            catch
                fprintf('call back doesn''t run\n')
            end
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
            spm_orthviews('Addtruecolourimage', h(1), file2, jet, 0.4 );
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
                title(sprintf('%s: %3.4g ? %3.4g (M?SEM)',fields{1},M,S),'interpreter','none');
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
                title(sprintf('%s: %3.4g ? %3.4g (M?SEM)',fields{1},M,S),'interpreter','none');
            end
        end
        function plot_spacetime(self,Y)
            
            %%
            Ys  = self.circconv2(Y);%smooth
            Ycs = cumsum(Ys)./repmat([1:size(Ys,1)]',1,8);            
            t   = self.fit_spacetime(Ycs);
            %%
            figure;set(gcf,'position',[25 7 831 918]);
            subplot(2,2,1);imagesc((Ys')');colorbar;axis xy;set(gca,'xtick',[4 8],'xticklabel',{'CS+' 'CS-'});
            subplot(2,2,2);imagesc((Ycs')');colorbar;axis xy;set(gca,'xtick',[4 8],'xticklabel',{'CS+' 'CS-'});
            subplot(2,2,3);bar(((t.fit_results{8}.params(:,1))));box off;title('amplitude')
            subplot(2,2,4);bar(((vM2FWHM(t.fit_results{8}.params(:,2)))));box off;title('fwhm')
        end
    end
    methods (Static)
        function plot_bar(Y)
            % Fear tuning are stored in columns.
            %%            
            tcol = size(Y,2);
            X    = linspace(-135,180,8)';
            X    = repmat(X,1,tcol) + repmat([0:360:360*(tcol-1)],8,1);
            h    = bar(X(:),Y(:),.9);
            cmap = GetFearGenColors;
            cmap = cmap(1:8,:);
            SetFearGenBarColors(h,repmat(cmap,tcol,1)');
            %
            hold on;
            set(gca,'xtick',X(:),'xticklabel',{'' '' '' 'CS+' '' '' '' 'CS-'});
            box off;
            set(gca,'color','none');
            xlim([-155.2500  200.2500+360*(tcol-1)]);            
            drawnow;            
            axis tight;box off;axis square;drawnow;alpha(.5);
            if tcol > 1
                mx = max(X)
                for ncol = 1:tcol
                    plot(repmat(mx(ncol)+45/2,1,2),ylim,'k-.')
                end
            end
        end
        function Yc = circconv2(Y,varargin)
            %circularly smoothes data using a 2x2 boxcar kernel;
            
            if nargin == 1
%                 kernel = make_gaussian2D(3,3,3,1.5,3,2);
%                   kernel = make_gaussian2D(7,3,4,1.5,7,2);
                  kernel = make_gaussian2D(7,3,10,1.5,7,2);
%                 kernel = make_gaussian2D(3,3,3,1.5,2,2);
%                 kernel = make_gaussian2D(5,3,4,1.5,3,2);
                kernel = kernel./sum(kernel(:));
                kernel = flipud(kernel);
            elseif nargin == 2
                kernel = varargin{1};
            end            
            [Yc]= conv2([Y Y Y],kernel);
            Yc  = Yc(size(kernel,1):end,10:17);

        end
        function [t2]=fit_spacetime(Y);
            data.y           = Y(:,1:8);
            data.x           = repmat(linspace(-135,180,8),size(data.y,1),1);
            data.ids         = 1:size(data.y,1);
            t2               = Tuning(data);
            t2.visualization = 1;
            t2.gridsize      = 10;
            t2.SingleSubjectFit(8);
%             t2.SingleSubjectFit(7);
%             t2.SingleSubjectFit(3);
        end
        function RunSPMJob(matlabbatch)
            %will run the spm matlabbatch using the parallel toolbox.
            fprintf('Will call spm_jobman...\n');
           
            for n = 1:length(matlabbatch)
                fprintf('Running SPM jobman %i...\n',n);
                spm_jobman('run', matlabbatch(n));
            end            
        end                        
    end
end
