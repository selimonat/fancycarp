classdef Project < handle
    % This is the PROJECT object that many of the other objects will
    % be a child of. Here enters all the project specific data (e.g. the
    % subject ids) and methods (for e.g. getting paths from the dicom
    % server).
    %
    % To start with, trio_sessions should be entered manually for your
    % experiment. Later these variables will be used by different methods
    % in the SUBJECT object to transfer data from the dicom server.
    %
    % To create the standard project structure you will need to call the
    % CreateFolderHierarchy method. Based on dicom2run, trio_sessions,
    % data_folders properties, the hierarchy will be created.
    %
    % Once you data is nicely loaded into the subject/run/ structure,
    % SanityCheck method comes very handy to ensure that your project
    % folders have all the size they should be.
    %
    % Feel free to improve this help section.
    %
    % please ensure that smrREADER, SPM12, Palamedes are in your path.
    
    properties (Hidden, Constant)
        %All these properties MUST BE CORRECT and adapted to one owns
        %project
        
        trio_sessions         = {  '' '' '' 'TRIO_17455' 'TRIO_17468' 'TRIO_17476' 'TRIO_17477' 'TRIO_17478' 'TRIO_17479' 'TRIO_17480' 'TRIO_17481' 'TRIO_17482' 'TRIO_17483' 'TRIO_17484' 'TRIO_17485' 'TRIO_17486' 'TRIO_17487' 'TRIO_17488' 'TRIO_17514' 'TRIO_17515' 'TRIO_17516' 'TRIO_17517'  'TRIO_17520' 'TRIO_17521' 'TRIO_17522' 'TRIO_17523' 'TRIO_17524' 'TRIO_17525' 'TRIO_17526' 'TRIO_17527' 'TRIO_17557' 'TRIO_17558' 'TRIO_17559' 'TRIO_17560'  'TRIO_17563' 'TRIO_17564' 'TRIO_17565' 'TRIO_17566' 'TRIO_17567' 'TRIO_17568' 'TRIO_17569' 'TRIO_17570' 'TRIO_17571' 'TRIO_17572'};
        dicom_serie_selector  = {  [] [] []   [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [5 6 7]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]       [3 4 5]       [3 4 5]      [3 4 5]      [3 4 5]    [3 4 5]       [3 4 5]       [3 4 5]     [3 4 5]     [4 5 6]       [3 4 5]      [3 4 5]     [3 4 5]       [3 4 5]      [3 4 5]        [3 4 5]     [3 4 5]       [3 4 5]      [3 4 5]       [3 4 5]     [3 4 5]     [4 5 6]      [3 4 5]    };
        %this is necessary to tell matlab which series corresponds to which
        %run (i.e. it doesn't always corresponds to different runs)
        dicom2run             = repmat({[1 1 1]},1,length(Project.dicom_serie_selector));
        data_folders          = {'eye' 'midlevel' 'mrt' 'scr' 'stimulation'};%        
        palamedes_path        = '/Users/onat/Documents/Code/Matlab/palamedes1_8_0/Palamedes/';
        spm_path              = '/Users/onat/Documents/Code/Matlab/spm12/';  
        TR                    = 0.99;
        path_project          = '/projects/fearamy/data/';
        path_stimuli          = '';%optional
    end
    properties (Constant) %project specific properties        
        condition_labels      = {'null' '1' '2' '3' '4' '5' '6' '7' '8' 'ucs' 'odd'};
        colors                = [ [0 0 0]; 0.0784 0.3284 1.0000;0.5784    0.0784    1.0000;1.0000    0.0784    0.8284;1.0000    0.0784    0.0784;1.0000    0.8284    0.0784;0.5784    1.0000    0.0784;0.0784    1.0000    0.3284;0.0784    1.0000    1.0000;0.0784    0.0784    0.0784;0.5784    0.5784    0.5784  ;[.8 0 0];[.8 0 0]];
        line                  = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '.' '.'};
        symbol                = {'.' '.' '.' '.' '.' '.' '.' '.' '.' 'p' 's'};
        total_subject         = [];
        font_style            = {'fontsize' 12};        
    end
    
    methods
        function DU = SanityCheck(self,runs,varargin)
            %will run through subject folders and will plot their disk
            %space. Use a string in VARARGIN to focus only on a subfolder
            total_subjects = length(Project.trio_sessions);
            DU = nan(total_subjects,length(runs));
            for ns = 1:total_subjects
                for nr = runs
                    fprintf('Requesting folder size for subject %03d and run %03d\n',ns,nr);
                    fullpath = fullfile(self.pathfinder(ns,nr),varargin{:});
                    if exist(fullpath)
                        [a b]       = system(sprintf('/usr/bin/du -cd0 %s | grep -e total | grep -oh -e [0-9]*',fullpath));
                        DU(ns,nr+1) = str2double(b)./1024;
                    else
                        fprintf('Folder doesn''t exist: %s\n',fullpath);
                    end
                end
            end
            %%
            n = 0;
            for nr = runs
                n = n + 1;
                subplot(3,1,n)
                bar(DU(:,n));ylabel('MB');xlabel('Subjects');box off
                title(sprintf('Subfolder %s\n Run: %03d',varargin{:},nr))
            end
            
        end
        function DicomDownload(self,source,destination)
            % Will download all the dicoms, convert them and merge them to
            % 4D.
            
            fprintf('%s:\n','DicomDownload');
            %downloads data from the dicom dicom server
            if exist(destination) == 0
                fprintf('The folder %s doesn''t exist yet, will create it...\n',destination)
                mkdir(destination)
            end
            if ~isempty(source) & ~isempty(destination)%both paths should not be empty.
                fprintf('Calling system''s COPY function to dump the data...%s\nsource:%s\ndestination:%s\n',self.gettime,source,destination)
                a            = system(sprintf('cp -r %s/* %s',source,destination));
                if a ~= 0
                    fprintf('There was a problem while dumping...try to debug it now\n');
                    keyboard
                else
                    fprintf('COPY finished successully %s\n',self.gettime)
                end
            else
                fprintf('Either source or destination is empty\n');
                return
            end            
        end
        function DicomTo4D(self,destination)
            %A wrapper over conversion and merging functions.
            
            %start with conversion
            if self.ConvertDicom(destination);
                %finally merge 3D stuff to 4D and rename it data.nii.
                self.MergeTo4D(destination);
            end
        end
        function MergeTo4D(self,destination)
            %will create data.nii consisting of all the [f,s]TRIO images
            %merged to 4D. the final name will be called data.nii.
            % merge to 4D
            
            matlabbatch = [];
            files       = spm_select('FPListRec',destination,'^[f,s]TRIO');
            fprintf('MergeTo4D:\nMerging (%s):\n%s\n',self.gettime,files);
            if ~isempty(files)
                matlabbatch{1}.spm.util.cat.vols  = cellstr(files);
                matlabbatch{1}.spm.util.cat.name  = 'data.nii';
                matlabbatch{1}.spm.util.cat.dtype = 0;
            end
            
            if ~isempty(matlabbatch)
                self.RunSPMJob(matlabbatch);
            end
            fprintf('Finished... (%s)\n',self.gettime);
            fprintf('Deleting 3D images in (%s)\n%s\n',self.gettime,destination);
            files = cellstr(files);
            delete(files{:});
        end
        function success    = ConvertDicom(self,destination)
            % dicom conversion. ATTENTION: dicoms will be deleted
            % and the converted. Assumes all files that start with MR are
            % dicoms.
            %
            success =0;
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
                               
                fprintf('Dicom conversion s#%i... (%s)\n',self.id,self.gettime);
                self.RunSPMJob(matlabbatch);
                fprintf('Finished... (%s)\n',datestr(now,'hh:mm:ss'));
                % delete dicom files
                fprintf('Deleting DICOM images in (%s)\n%s\n',self.gettime,destination);
                files = cellstr(files);
                delete(files{:});
                fprintf('Finished... (%s)\n',self.gettime);
                success = 1;
            else
                fprintf('No dicom files found for %i\n',self.id)                
                fprintf('No file found...\n')
            end
        end        
        function [result]   = dicomserver_request(self)
            %will make a normal dicom request. use this to see the state of
            %folders
            fprintf('Making a dicom query, sometimes this might take long (so be patient)...(%s)\n',self.gettime);
            [status result]  = system(['/common/apps/bin/dicq --series --exam=' self.trio_session]);
            fprintf('This is what I found for you:\n');
            result
        end
        function [paths]    = dicomserver_paths(self)
            fprintf('Making a dicom query, sometimes this might take long (so be patient)...(%s)\n',self.gettime)
            [status paths]   = system(['/common/apps/bin/dicq -t --series --exam=' self.trio_session]);
            paths            = strsplit(paths,'\n');%split paths
        end
        function [data_path]= pathfinder(self,subject,run)
            %gets the path
            %Example: s.pathfinder(s.id,[]) => will return the path to
            %subject
            % empty [] above can be replaced with any phase number.
            data_path = self.path_project;
            for no = [subject run]
                file_list        = dir(data_path);
                i                = regexp({file_list(:).name},sprintf('[0,a-Z]%d$',no));
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
    end
    methods(Static)
        function t          = gettime
            t = datestr(now,'hh:mm:ss');
        end
        function RunSPMJob_inParallel(matlabbatch)
            %will run the spm matlabbatch using the parallel toolbox.
            fprintf('Will call spm_jobman in parallel with 4 cores...\n');
            if isempty(gcp)
                parpool(4);
            end
            parfor n = 1:length(matlabbatch)
                fprintf('Running SPM jobman %i...\n',n);
                spm_jobman('run', matlabbatch(n));
            end
            delete(gcp);
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
    methods %project specific methods
        function degree    = stimulus2degree(self,stim_id)
            %will transform condition indices to distances in degrees from
            %the csp face. stim_id is a cell array.
            
            ind_valid     = find(cellfun(@(x) ~isempty(x),regexp(stim_id,'[0-9]')));
            degree        = stim_id;
            for i = ind_valid(:)'
                degree{i} = mat2str(MinimumAngle( 0 , (stim_id{i}-self.csp)*45 ));
            end
        end
        function set_feargen_colors(h,color_ind);
            %if H is a handle of a barplot, it will colorize it with
            %typical feargen colors.
            h     = get(h,'children');
            tbar  = length(get(h,'YData'));
            set(h,'CData', repmat(1:tbar,1,tbar/tbar),'edgecolor','none');
            colormap(Project.colors(color_ind,:));
        end
    end
end
