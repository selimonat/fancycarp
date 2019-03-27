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
       
%         path_project        = '/Users/onat/Documents/BehavioralExperiments/fearcloud/';
%         path_project      = sprintf('%s%sGoogle Drive%sEthnoMaster%sBDNF%s',homedir,filesep,filesep,filesep,filesep);
%         path_project      = sprintf('%s%sDocuments%sExperiments%sBDNF%sdata%s',homedir,filesep,filesep,filesep,filesep,filesep)
        path_project      = sprintf('%s%sDocuments%sExperiments%sFearCloud_Eyelab%sdata%s',homedir,filesep,filesep,filesep,filesep,filesep)
%          path_project      = sprintf('%s%sDocuments%sExperiments%sPlaPil%sdata%s',homedir,filesep,filesep,filesep,filesep,filesep)
        scr_blocknames    = {'test_rating' 'test' 'cond_rating' 'cond' 'base_rating' 'base' };

        path_stimuli      = sprintf('%sstimuli%s',Project.path_project,filesep);
        ETMaskpath        = sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep)
        condition_labels  = {'null' '1' '2' '3' '4' '5' '6' '7' '8' 'ucs' 'odd'};
        plot_style
        subjects         = [ 6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 40 42 43 44 45 46 47 48 49 50 51 52 53 54 55  56  57 58 59  60 61 62 63 64 65];
        age              = [29 30 nan 29 27 27 35 26 26 26 21 35 33 19 29 19 23 21 34 28 27 35 26 20 22 27 26 23 31 26 19 20 23 24 23 19 29 30 25 24 25 27 21 34 28 18 34 nan 25 nan 24 26 nan 24 26 27 25 29];
        subjects_600     = [ 0	0	0  0  0	 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1  1  1   1  1  1  1  1  1];
        subjects_1500    = [ 1	1	1  1  1	 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0  0  0   0  0  0  0  0  0];
        ET_fg            = [ 1  0   1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1   1   1  1  1   1  1  1  1  1  1];
        ET_pmf           = [ 1  1   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1   1   1  1  0   1  1  1  1  1  1];
        subjects_scr     = [ 1  1   1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1  1  1   1  1  1  1  1  1];
        subjectID_600      = [27,37:40,42:65]; %absolute numbers if needed
        subjectID_1500     = [6:26,28,30:36];  %absolute numbers if needed
        gender            = ones(length(Project.subjects),1);%all male anyway
        
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
            o = fullfile(FPSA_FearGen_MSc('get_path_project'),'data',filesep);
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
                errorbar(X,Y,SEM,'k.','LineWidth',1.5);%add error bars
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
