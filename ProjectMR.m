classdef ProjectMR < handle
    properties (Hidden, Constant)
        
    end
    properties (Constant)
        path_project       = sprintf('%s%sDesktop%sfearamy/',homedir,filesep,filesep);        
        trio_names         = { 'TRIO_12476' };
        trio_folders       = {   [1:22]    };
        subjects           = [1]
    end
    
    methods        
    end
    
    methods        
        function stim = find_stim(self,varargin)
            %will return the path to the Nth (varargin) stimulus. If not
            %specified the average stim will be returned.                       
            if ~isempty(varargin)
                stim = sprintf('%s%02d.bmp',Project.path_stimuli,n);            
            else
                stim = sprintf('%save.bmp',Project.path_stimuli);
            end
        end                
        function [data_path]=pathfinder(self,subject,run)
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
        function path2data = path2data(self,subject,run,varargin)
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