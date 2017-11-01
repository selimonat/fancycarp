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
        path_project      = [];%refer to get function at the bottom

        scr_blocknames    = {'test_rating' 'test' 'cond_rating' 'cond' 'base_rating' 'base' };

        path_stimuli      = sprintf('%sstimuli%s',Project.path_project,filesep);
        ETMaskpath        = sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep)
        condition_labels  = {'null' '1' '2' '3' '4' '5' '6' '7' '8' 'ucs' 'odd'};
        plot_style        
        BDNF              = [2 1 1 1 2 1 2 1 2  2  1  1  1  2 1  1  2  2  2  2  1  1  1  2  1  1  1  1  1  1  1  2  2    2  2  2  2  2  2     2  1  1  1  1  1  2  1  1  2  1  1  2   1  1  2  1  1  1  2  2  1  1  1  1  1  2  2  2  2  1  2  1 NaN 1  1  2  1   2  2  2  2];
        subjects          = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33   35 36 37 38 39 40   41 42 43 44 45 46 47 48 49  50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78  79 80 81 82];
        subjects_ET       = [1 1 1 1 1 1 1 1 1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    1  1  1  0  1  1    1  1  1  1  1  1  1  1  1   1  1  0  1  1  0  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  1  1   1  1  1  1] == 1;
        subjects_scr      = [1 0 1 1 1 1 1 1 1  1  1  1  0	1  1  1  1	0  1  1  1  1  1  1  1  1  1  1  1  1  0  0  1    1  1  1  1  0  1    1  1  1  1  1  1  1  1  1   0  1	1  1  1  1  1  0  0  1  0  1  1  0  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1   1  1  1  0] == 1;
        cspface           = [1 1 2 3 2 4 3 5 4  5  6  7  8  6  1  2  7  8  2  8  1  2  3  1  4  5  6  7  8  5  4  5  2    3  4  6  5  6  1    3  3  7  1  5  2  7  6  4   8  8  3  1  4  7  1  1  5  2  5  2  4  6  6  8  8  6  4  3  7  2  5  3  2  5  6  2  7   8  6  4  7];   
        %the following questionnaire data is from the Z project, use those for now.
        CTQ_SFB           = [1 0 0 0 1 1 0 0 0  0  1  0  1  0  0  0  1  0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  1    0  0  0  0  1  0    1  0  1  0  0  1  1  0  0   0  0  0  0  0  0 NaN 0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  1  0  0   0  0  0  0];
        LTE_SFB           = [1 0 1 1 0 0 0 1 0	1  1  1  1	1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  0  1    1  1	1  0  1  1    1  1  1  1  1  0  1  1  0   1  1  1  1  1  1 NaN 1  1  1  1  1  1  0  1  1  1  1  1  1  0  1  1  1  1  1  0  1  1   1  1  1  0];
        % those come from BDNF data collection
        CTQ_BDNF          = [1 0 0 0 1 1 0 0 0  0  1  0  1  0  0  0  0  0  0  0  0  1  1  0  0  1  0  0  0  0  0  0  0    0  1  0  0  1  0    1  0  1  0  1  1  1  1  0   0  0  1  0  0  0  1  0  0  1  1  1  0  0  0  0  0  0  1  0  1  0  1  0  0  1  1  0  0   0  0 NaN NaN];
        LTE_BDNF          = [1 0 1 1 1 1 0 1 0  0  1  0  0  1  1  1  1  0  1  0  1  1  1  1  1  1  1  0  1  1  0  0  1    1  1  1  1  1  1    1  0  1  1  1  0  1  0  1   1  1	1  0  1  1  1  1  0  1  1  1  1  0  1  0  1  1  1  1  0  1  1  1  1  1  0  1  1   1  1	1  0];
        
        age               = [20	28	25	27	26	26	30	23	23	32	32	26	27	26	21	22	33	26	27	23	26	31	29	25	31	36	33	34	30	25	25	25	25	20	31	22	29	32	25	25	26	27	35	22	22	26	28	28	23	23	21	26	25	25	30	28	25	31	28	30	34	24	21	24	24	25	28	28	22	21	28	29	25	25	31	23	34	28	32	36	31];

        gender            = [ones(39,1); ones(42,1)*2];
        
    end
    
    methods
        
        function stim = find_stim(self,varargin)
            %will return the path to the Nth (varargin) stimulus. If not
            %specified the average stim will be returned.                       
            if ~isempty(varargin)
                stim = sprintf('%s%02d.png',Project.path_stimuli,varargin{:});            
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
    
    
        function o = get.path_project
            %
            o = fullfile(FPSA_FearGen('get_path_project'),'data',filesep);
            %add a filesep if necessary
            if ~strcmp(o(end),filesep)
                cprintf([1 0 0],'Correcting project path...\n');
                o(end+1) = filesep;
            end
        end
    end
    methods (Static)
        function plot_bar(X,Y,SEM)
            % input vector of 8 faces in Y, angles in X, and SEM. All
            % vectors of 1x8;
            %%
            cmap  = GetFearGenColors;
            if length(Y)==11;
                cmap = [cmap; [0 0 0 ]];
            end
            tbar  = length(Y);
            for i = 1:tbar
                try
                    h(i)    = bar(X(i),Y(i),40,'facecolor',cmap(i,:),'edgecolor','none','facealpha',.8);
                catch
                    h(i)    = bar(X(i),Y(i),40,'facecolor',cmap(i,:),'edgecolor','none');
                end
                hold on;
            end
            %%
            hold on;
            if nargin == 3
                errorbar(X,Y,SEM,'k.');%add error bars
            end
            
            %%
            set(gca,'xtick',X,'xticklabel',{'' '' '' 'CS+' '' '' '' 'CS-'});
            box off;
            set(gca,'color','none');
            xlim([0 tbar+1])
            drawnow;
            axis tight;box off;axis square;drawnow;alpha(.5);
            
        end
        
        
    end
    
end
