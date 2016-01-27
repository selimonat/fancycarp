classdef ProjectMR < handle
    properties (Hidden, Constant)
        trio_sessions       = { 'TRIO-17370' 'TRIO_17399' };
        dicom_folders       = {   [1:22]   [6 7 8] };
        dicom2run           = {   [1:22]   [1 1 1] };
        trio_ids            = { 'V13132' 'V9434'};
    end
    properties (Constant)
        path_project       = sprintf('%s%sDesktop%sfearamy/',homedir,filesep,filesep);        
        subjects           = [1 2];
        TR                 = 0.99;
    end
    
    methods        
    end
    
    methods        
        function DicomDownload(self,source,destination)
            
            %downloads data from the dicom dicom server
            if exist(destination) == 0
                fprintf('The folder %s doesn''t exist yet, will create it...\n',destination)
                mkdir(destination)
            end
            fprintf('Calling system''s COPY function to dump the data...\n')
            a            = system(sprintf('cp -vr %s/* %s',source,destination));
            a = 0;
            if a ~= 0
                fprintf('There was a problem while dumping...\n');
                keyboard
            end
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
        function path2data  = path2data(self,subject,run,varargin)
            % s.path2data(53,4) will return the path to the subject's phase 4
            % s.path2data(53,4,'eye') return the path to the eye data file at the
            % 4th phase.
            
            %will return the path to phase/data_type/
            path2data = self.pathfinder(subject , run);
            if length(varargin) >= 1 
                path2data = sprintf('%s%s%sdata.mat',path2data,varargin{1},filesep);
            end
            if length(varargin) == 2
                path2data = regexprep(path2data,'mat',varargin{2});
            end
        end                                       
        
    end    
end
