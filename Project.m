classdef Project < handle
    properties (Hidden, Constant)
        colors            = [ [0 0 0];circshift( hsv(8), [3 0] );[.8 0 0];[.8 0 0]];
        line              = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '.' '.'};
        symbol            = {'.' '.' '.' '.' '.' '.' '.' '.' '.' 'p' 's'};
        marker_size       = {20 20 20 20 20 20 20 20 20 20 20};
        line_width        = {2  2  2  2  2  2  2  2  2  1  1};
        condition_indices = {1000 45 90 135 180 225 270 315 360 1001 1002};
        screen_resolution = [1200 1600];
        PixelPerDegree    = 37;
    end
    properties (Hidden,Constant)

        path_project      = sprintf('%s%sDocuments%sBehavioralExperiments%sfearcloud%s',homedir,filesep,filesep,filesep,filesep)
%         path_project      = sprintf('%s%sGoogle Drive%sEthnoMaster%sBDNF%s',homedir,filesep,filesep,filesep,filesep);
%         scr_blocknames    = {'test_rating' 'test' 'cond_rating' 'cond' 'base_rating' 'base' };
        scr_blocknames    = {'test_rating' 'test' 'cond_rating' 'cond' 'base_rating' 'base' };

        path_stimuli      = sprintf('%sstimuli%s',Project.path_project,filesep);
        ETMaskpath        = sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep)
        condition_labels  = {'null' '1' '2' '3' '4' '5' '6' '7' '8' 'ucs' 'odd'};
        plot_style
        subjects_600      = [27,37:40,42:65];
        subjects_1500     = [6:26,28,30:36];        
        BDNF              = [2 1 1 1 2 1 2 1 2 2   1  1  1  2 1  1  2  2  2  2  1  1  1  2  1  1  1  1  1  1  1  2  2  2  2  2  2  2  2  2  2  1  1  1  1   1 2  1  1  2  1  1  2  1  1  2  1  1  1  2  2  1  1  1  1  1  2  2];
        subjects_bdnf     = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68];
        gender            = [ones(40,1); ones(40,1)*2];
        
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
        
        function out = getMask(self,varargin)
            %varargin is ET_feargen, ET_discr, SCR, PMF,RATE          
            a = load(sprintf('%smidlevel%ssubjmasks%s%smask.mat',self.path_project,filesep,filesep,varargin{1}));
            dummy = fieldnames(a);
            out   = a.(dummy{1});
        end
        
        function [o]=tRuns(self)
            %% returns the total number of runs in a folder
            [~, d] = spm_select('FPList',s.path,'^run');
            o      = length(d);
        end
                
    end
    
end
