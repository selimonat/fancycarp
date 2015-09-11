classdef Subject < Project
    properties (Hidden)
        paradigm
    end
    properties (SetAccess = private)
        id
        path
        csp
        csn
        scr        
    end
    methods
        function s = Subject(id)%constructor
            s.id              = id;
            s.path            = s.pathfinder(s.id,[]);
            if exist(s.path)
                for nrun = 1:5
                    s.paradigm{nrun} = s.load_paradigm(nrun);
                end                
                s.csp = s.paradigm{2}.stim.cs_plus;
                s.csn = s.paradigm{2}.stim.cs_neg;                                             
%                 s.scr = SCR(s);                
            else                
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path)
            end
        end
    end
    
    methods
        
        function p         = load_paradigm(self,nrun)
            %HAST TO GO TO THE PROJECT ACTUALLY TOGETHER WITH
            %CONDTION_>COLOR DESCRIPTION
            filename = self.path2data(self.id,nrun,'stimulation');
            p = [];
            if exist(filename)
                p = load(filename);
                p = p.p;
                %transform id to labels
                if isfield(p,'presentation')
                    %                 p.presentation.stim_label = self.condition_labels(p.presentation.cond_id+1);
                    p.presentation.dist(p.presentation.dist < 500)   = p.presentation.dist(p.presentation.dist < 500) + 180;
                    p.presentation.dist(p.presentation.dist == 500)  = 1001;%ucs
                    p.presentation.dist(p.presentation.dist == 1000) = 1002;%odd
                    p.presentation.dist(isnan(p.presentation.dist))  = 1000;%null
                end
            end
        end
        function degree    = stimulus2degree(self,stim_id)
            %will transform condition indices to distances in degrees from
            %the csp face. stim_id is a cell array. This is a subject
            %method as it depends on the subject specific CSP face.
            
            ind_valid     = find(cellfun(@(x) ~isempty(x),regexp(stim_id,'[0-9]')));
            degree        = stim_id;
            for i = ind_valid(:)'
                degree{i} = mat2str(MinimumAngle( 0 , (stim_id{i}-self.csp)*45 ));
            end
        end
        function color     = condition2color(self,cond_id)
            cond_id/45+4;
        end
        function rating    = GetRating(self,run,align)
            %align optional, default = 1;
            if nargin < 3
                align = 1;
            end
            % s is a subject instance
            rating = [];
            if ~isempty(self.paradigm{run})
                rating.y      = self.paradigm{run}.out.rating';
                if align
                    rating.y  = circshift(rating.y,[1 4-self.csp ]);
                end
                rating.x      = repmat([-135:45:180],size(self.paradigm{run}.out.rating,2),1);
                rating.i      = repmat(run          ,size(self.paradigm{run}.out.rating,2),size(self.paradigm{run}.out.rating,1));
                rating.y_mean = mean(rating.y);
            else
                fprintf('no rating present for this subject and run (%d) \n',run);
            end            
        end
            
    end
end