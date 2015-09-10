classdef Project < handle
    properties (Hidden, Constant)
        colors            = [ [0 0 0];circshift( hsv(8), [3 0] );[.4 .4 .4];[.8 .8 .8]];
        line              = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '.' '.'};
        symbol            = {'.' '.' '.' '.' '.' '.' '.' '.' '.' 'p' 's'};
        marker_size       = {20  20  20  20  20  20  20  20  20  10  10};
        line_width        = {2  2  2  2  2  2  2  2  2  1  1};
        condition_indices = {1000 45 90 135 180 225 270 315 360 1001 1002};
        screen_resolution = [1200 1600];
        PixelPerDegree    = 20;
    end
    properties (Hidden,Constant)
        path_project      = '/Users/onat/Google Drive_ninhamburg/EthnoMaster/BDNF/';
        condition_labels  = {'null' '1' '2' '3' '4' '5' '6' '7' '8' 'ucs' 'odd'};
        plot_style
        subjects_600      = [27,37:65];
        subjects_1500     = [6:26,28:36];        
        BDNF              = [2 1 1 1 2 1 2 1 2 2 1 1 1 2 1 1 2 2 2 2 1 1 1 2 1 1 1 1 1 1 1];
    end
    
    methods
        
        function ls = get.plot_style()
            %produce a cell array of string for each condition that specifies plotting
            %attributes.
            for ncond = [1:11];
                ls(Project.condition_indices{ncond}).line        = {{'line',  Project.line{ncond}}};
                ls(Project.condition_indices{ncond}).color       = {{'color', Project.colors(ncond,:)}};
                ls(Project.condition_indices{ncond}).symbol      = {{'marker',Project.symbol{ncond}}};
                ls(Project.condition_indices{ncond}).marker_size = {{'markersize',Project.marker_size{ncond}}};
                ls(Project.condition_indices{ncond}).line_width = {{'linewidth',Project.line_width{ncond}}};
            end
        end
        
        function [data_path]=pathfinder(self,subject,run)
            %gets the path
            %Example: s.pathfinder(s.id,[]) => will return the path to
            %subject
            % empty [] above can be replaced with any phase number.
            % 
            
            data_path = self.path_project;
            for no = [subject run]
                file_list        = dir(data_path);
                i                = regexp({file_list(:).name},sprintf('[0,a-Z]%d$',no));
                i                = find(cellfun(@(bla) ~isempty(bla), i  ));
                if ~isempty(i)
                    folder       = getfield(file_list(i),'name');
                    data_path    = sprintf('%s/%s',data_path,folder);
                else
                    data_path    = [];
                    return
                end
            end
            data_path(end+1)         = '/';
        end
        
        function path2data = path2data(self,subject,run,varargin)
            % s.path2data(53,4) will return the path to the subject's phase 4
            % s.path2data(53,4,'eye') return the path to the eye data file at the
            % 4th phase.
            
            %will return the path to phase/data_type/
            path2data = self.pathfinder(subject , run);
            if length(varargin) >= 1 
                path2data = sprintf('%s%s/data.mat',path2data,varargin{1});
            end
            if length(varargin) == 2
                path2data = regexprep(path2data,'mat',varargin{2});
            end
        end
    end
    
end