classdef Project < handle
    %% please ensure that smrREADER, SPM12, Palamedes are in your path.
    
    properties (Hidden, Constant)
        trio_sessions         = {  '' '' '' 'TRIO_17455' 'TRIO_17468' 'TRIO_17476' 'TRIO_17477' 'TRIO_17478' 'TRIO_17479' 'TRIO_17480' 'TRIO_17481' 'TRIO_17482' 'TRIO_17483' 'TRIO_17484' 'TRIO_17485' 'TRIO_17486' 'TRIO_17487' 'TRIO_17488' 'TRIO_17514' 'TRIO_17515' 'TRIO_17516' 'TRIO_17517'  'TRIO_17520' 'TRIO_17521' 'TRIO_17522' 'TRIO_17523' 'TRIO_17524' 'TRIO_17525' 'TRIO_17526' 'TRIO_17527' 'TRIO_17557' 'TRIO_17558' 'TRIO_17559' 'TRIO_17560'  'TRIO_17563' 'TRIO_17564' 'TRIO_17565' 'TRIO_17566' 'TRIO_17567' 'TRIO_17568' 'TRIO_17569' 'TRIO_17570' 'TRIO_17571' 'TRIO_17572'};
        dicom_serie_selector  = {  [] [] []   [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [5 6 7]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]      [3 4 5]       [3 4 5]       [3 4 5]      [3 4 5]      [3 4 5]    [3 4 5]       [3 4 5]       [3 4 5]     [3 4 5]     [4 5 6]       [3 4 5]      [3 4 5]     [3 4 5]       [3 4 5]      [3 4 5]        [3 4 5]     [3 4 5]       [3 4 5]      [3 4 5]       [3 4 5]     [3 4 5]     [4 5 6]      [3 4 5]    };
        %this is necessary to tell matlab which series corresponds to which
        %run (i.e. it doesn't always corresponds to different runs)
        dicom2run             = repmat({[1 1 1]},1,length(Project.dicom_serie_selector));
        colors                = [ [0 0 0];circshift(hsv(8),[3 0]);[.8 0 0];[.8 0 0]];
        line                  = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '.' '.'};
        symbol                = {'.' '.' '.' '.' '.' '.' '.' '.' '.' 'p' 's'};
        total_subject         = [];
        palamedes_path        = '/Users/onat/Documents/Code/Matlab/palamedes1_8_0/Palamedes/';
    end
    properties (Constant)
        path_project       = '/projects/fearamy/data/';
        path_stimuli       = '';
        condition_labels   = {'null' '1' '2' '3' '4' '5' '6' '7' '8' 'ucs' 'odd'};
        TR                 = 0.99;
    end
    
    methods
    end
    
    methods
         
        function DU = SanityCheck(self,runs,varargin)
            %will run through subject folders and will plot their disk
            %space. Use a string in VARARGIN to focus only on a subfolder            
            DU = nan(self.total_subject,length(runs));
            for ns = 1:self.total_subject
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
                    fprintf('There was a problem while dumping...\n');
                    keyboard
                else
                    fprintf('COPY finished successully %s\n',self.gettime)
                end
            else
                fprintf('Either source or destination is empty\n');
            end
        end
        function [result]=dicomserver_request(self)
            %will make a normal dicom request. use this to see the state of
            %folders
            fprintf('Making a dicom query, sometimes this might take long (so be patient)...(%s)\n',self.gettime);
            [status result]  = system(['/common/apps/bin/dicq --series --exam=' self.trio_session]);
            fprintf('This is what I found for you:\n');
            result
        end
        function [paths]=dicomserver_paths(self)
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
        
              
        function degree    = stimulus2degree(self,stim_id)
            %will transform condition indices to distances in degrees from
            %the csp face. stim_id is a cell array. 
            
            ind_valid     = find(cellfun(@(x) ~isempty(x),regexp(stim_id,'[0-9]')));
            degree        = stim_id;
            for i = ind_valid(:)'
                degree{i} = mat2str(MinimumAngle( 0 , (stim_id{i}-self.csp)*45 ));
            end
        end
        function out = get.total_subject()
           %returns the number of total subject in the project folder
            out = length(Project.trio_sessions);
        end
        
        
    end
    methods(Static)
        function t          = gettime
            t = datestr(now,'hh:mm:ss');
        end
        
        function RunSPMJob(matlabbatch)
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
        
        function color     = condition2color(cond_id)
            %returns the color associated to a condition
            color = cond_id/45+4;
        end
        
    end
    
end
