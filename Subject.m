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
        pmf
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
                s.scr = SCR(s);
                s.pmf = s.getPMF;
            else
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path)
            end
        end
    end
    
    methods
        function out = getPMF(self)
            load(sprintf('%smidlevel%sweibull%sdata.mat',Project.path_project,filesep,filesep));
            out = data(self.id);
            %lines in the output are 
            %1/ CS+ before
            %2/ CS- before
            %3/ CS+ after
            %4/ CS- after
            out.subject_alpha = mean(out.params1(1:2,1),1);
            out.subject_beta  = mean(out.params1(1:2,2),1);
        end
        
        function pmfplot(self)
            plotpath = sprintf('%s%sp05%sfigures%sfearcloud_FitPMFs_RE.fig',self.path,filesep,filesep,filesep);
            if exist(plotpath)
            openfig(plotpath);
            else
                fprintf('no figure found!')
            end
        end
                   
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
        function out    = GetSubSCR(self,run,cond)
            if nargin < 3
                cond=1:8;
            end
            conddummy=[-135:45:180 500 1000 3000];
            % s is a subject instance
            out = [];
            cutnum = self.scr.findphase(run);
            self.scr.cut(cutnum);
            self.scr.run_ledalab;
            self.scr.plot_tuning_ledalab(cond);
            out.y = self.scr.fear_tuning;
            out.x = conddummy(cond);
            out.ind = cutnum;
        end
    end
end
