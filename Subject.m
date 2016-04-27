classdef Subject < Project
    %   Getter Methods are used to get stuff e.g. ratings, epis, stimulation
    %   log, behavioral data etc. It also implements downloading data from the
    %   dicom server.
    %
    %   path_tools Methods are used to return the path to different data
    %   types, not really interesting, used internally.
    %
    %   plotter Methods plot subject specific data.
    %
    %   fmri analysis methods contain first-level modelling.
    %
    %
    properties (Hidden)
        paradigm
        default_run   = 3;
        align_ratings = 1;
    end
    properties (SetAccess = private)
        id
        path
        csp
        csn
        ratings
    end
    %%
    methods
        function s = Subject(id)%constructor
            s.id              = id;
            s.path            = s.pathfinder(s.id,[]);
            if exist(s.path)
                runs = dir(s.path);
                a    = regexp({runs.name},'run[0-9][0-9][0-9]','match');
                a    = a(cellfun(@(x) ~isempty(x),a));%take all non empty entries;                
                for runs = a;                    
                    nrun             = regexp(runs{1}{1},'[0-9][0-9][0-9]','match');
                    nrun             = str2num(nrun{1});
                    s.paradigm{nrun} = s.load_paradigm(nrun);
                end
                s.csp = s.paradigm{end}.stim.cs_plus;
                s.csn = s.paradigm{end}.stim.cs_neg;
            else
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path)
            end
        end
    end
    
    methods %(Getters)
        
        function p          = load_paradigm(self,nrun)            
            filename = self.path2data(nrun,'stimulation');
            p = [];            
            p = load(filename);
            p = p.p;               
            
        end
        
        function rating = get.ratings(self)
            %returns the CS+-aligned ratings for all the runs
            for run = 1:3
                if ~isempty(self.paradigm{run})                    
                    rating(run).y      = self.paradigm{run}.out.rating';
                    if self.align_ratings
                        rating(run).y      = circshift(rating(run).y,[1 4-self.csp ]);%center it to CS+
                    end
                    rating(run).x      = repmat([-135:45:180],size(self.paradigm{run}.out.rating,2),1);
                    rating(run).ids    = repmat(self.id,size(self.paradigm{run}.out.rating,2),size(self.paradigm{run}.out.rating,1));
                    rating(run).y_mean = mean(rating(run).y,1);                    
                end                
            end
        end
    end
    
    methods %fmri path_tools which are related to the subject
        function path2data  = path2data(self,run,varargin)
            % s.path2data(53,4) will return the path to the subject's phase 4
            % s.path2data(53,4,'eye') return the path to the eye data file at the
            % 4th phase.
            
            if nargin < 2
                fprintf('you have to have at least one input for me...\n');
                return
            end
            %will return the path to phase/data_type/            
            path2data = self.pathfinder(self.id , run);
            if length(varargin) >= 1
                path2data = sprintf('%s%s%sdata.mat',path2data,varargin{1},filesep);
            end
            if length(varargin) == 2
                path2data = regexprep(path2data,'mat',varargin{2});
            end
        end
    end
end
