classdef Subject < Project
    properties (Hidden)
        paradigm
        default_run     = 2;
        mean_correction = 0;%decides if mean correction should be applied
        align           = 1;%should ratings be aligned to CS+ face
        fit_method      = 8;%defines method for Tuning Fit on feargen profiles
    end
    properties (SetAccess = private)
        id
        path
        csp
        csn
        scr
        feargen_rating
        feargen_scr
        groupinfo
    end
    methods
        function s = Subject(id)%constructor
            
            s.id              = id;
            s.path            = s.pathfinder(s.id,[]);
            if exist(s.path)
                fprintf('Subject %02d (%s)\n',id,s.path);
                for nrun = 1:5
                    s.paradigm{nrun} = s.load_paradigm(nrun);
                end
                s.csp = s.paradigm{s.default_run}.stim.cs_plus;
                s.csn = s.paradigm{s.default_run}.stim.cs_neg;
                
                s.scr            = SCR(s);
                                
%                 s.feargen_rating = s.get_fit('rating');
                s.feargen_scr    = s.get_fit('scr');
                
                s.groupinfo = s.get_groupinfo;
                                
            else
                fprintf('Subject %02d doesn''t exist somehow :(\n %s\n',id,s.path)
            end
        end
    end
    
    methods
        function out    = get_groupinfo(self)
            load(sprintf('%smidlevel%sgroups.mat',self.path_project,filesep));
            out.groups = group(self.id,:);
            out.tags   = tags;
        end
        function rating    = get_rating(self,run)
            % returns raw ratings in a format that is compatible with
            % Tuning object
            rating.y  = [];
            rating.x  = [];
            rating.ids = self.id;
            if ~isempty(self.paradigm{run})
                rating.y      = self.paradigm{run}.out.rating';
                if self.align
                    rating.y  = circshift(rating.y,[1 4-self.csp ]);
                end
                rating.y      = rating.y(:)';
                rating.x      = sort(repmat([-135:45:180],1,size(self.paradigm{run}.out.rating,2)));
                if self.mean_correction
                    rating.y = rating.y - mean(rating.y);
                end
            else
                warning('no rating present for this subject and run (%d) \n',run);
            end
        end               
        function out    = get_fit(self,modality) 
            %returns fit results for ratings or scr, modality is a string.
            forcefit = 1;
            if strcmp(modality,'rating');
                fun = @(x) self.get_rating(x);
                phases = 2:4;
            elseif strcmp(modality,'scr');
                fun = @(x) self.get_scr(x);
                phases = [2 4];
            end
            %
            fit_results =[];%will be filled in below;
            for ph = phases
                %
                phpath          = sprintf('%s%sp0%g%smidlevel%s%sfit_method_%d.mat',self.path,filesep,ph,filesep,filesep,modality,self.fit_method);
                if exist(phpath)&&~forcefit
                    cprintf([0 1 0],'Rating data cached for phase %02d, modality %s, will load it...\n',ph,modality);
                    load(phpath);
                    out(ph)     = fit_results;
                elseif ~exist(phpath) || forcefit
                    if forcefit
                        cprintf([1 0 0],'Forced refitting data, starting procedure with method %g...\n',self.fit_method)
                    end
                    if ~exist(phpath)
                        fprintf('No parameters found for phase %d (modality: %s), starting fitting procedure with method %g...\n',ph,modality,self.fit_method)
                    end
                    data                = fun(ph);
                    if ~isempty(data.y)
                        t               = Tuning(data);
                        t.visualization = 1;
                        t.gridsize      = 10;
                        fprintf('Phase %g... ',ph);
                        t.SingleSubjectFit(self.fit_method);
                        fit_results     = t.fit_results;
                        save(phpath,'fit_results');
                    end
                else
                    cprintf([1 0 0],'No %s data found for phase %02d to load or fit, doing nothing...\n',modality,ph);
                end
                %store the results for outputting and prepare the table                
                if ~isempty(fit_results)
                    out(ph)  = fit_results;
                    %adapt the parameter name to phase and modality;
                    for nnn = out(ph).param_table.Properties.VariableNames;
                        out(ph).param_table.Properties.VariableNames{nnn{1}} = sprintf('%s_%s_%02d',modality,nnn{1},ph);
                    end                
                else
                    out(ph) = struct();
                end
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
        function scr       = get_scr(self,run)
            %similar function to get_rating it will return the scr data
            %exactly in the same format suitable for Tuning object
            
            scr.x      = [-135:45:180];
            scr.ids    = self.id;
            scr.y      = [];
            if self.scr.ok;
                indices                         = {'' [1:8] [9:16] [17:24]};
                [out_z, out_raw, single_trials] = self.scr.ledalab_summary;
                
                %scr.y            = out_raw(indices{run});
                scr.y    = single_trials(:,indices{run});                
                scr.x    = repmat(scr.x,[size(scr.y,1) 1]);
                scr.y    = scr.y(:)';
                i        = isnan(scr.y);
                scr.y(i) = [];
                scr.x(i) = [];
            else                
                scr.y            = [];                
                warning('no scr present for this subject and run (%d) \n',run);
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
        
        function out    = GetSubSCRgraphs(self,run,cond)
            if nargin < 3
                cond=1:8;
            end
            conddummy=[-135:45:180 500 1000 3000];
            % s is a subject instance
            out = [];
            cutnum = self.scr.findphase(run);
            if length(self.scr.BlockNames)> 1
                self.scr.cut(cutnum);
            else
                warning('SCR was already cut, please check if correct phase!')
            end
            self.scr.run_ledalab;
            out.y = self.scr.ledalab.mean(:,cond);
            out.conds = conddummy(cond);
            out.phase = run;
        end
        function out    = GetSubSCRbars(self,run,cond)
            if nargin < 3
                cond=1:8;
            end
            conddummy=[-135:45:180 500 1000 3000];
            % s is a subject instance
            out = [];
            cutnum = self.scr.findphase(run);
            if length(self.scr.BlockNames)> 1
                self.scr.cut(cutnum);
            else
                warning('SCR was already cut, please check if correct phase!')
            end
            self.scr.run_ledalab;
            self.scr.plot_tuning_ledalab(cond);
            out.y = self.scr.fear_tuning;
            out.x = conddummy(cond);
            out.ind = cutnum;
            close(gcf);
        end
        function [o]    = tRuns(self)
            %% returns the total number of runs in a folder
            a      = dir(self.path);
            %include all runs except run000
            o      = cellfun( @(x) ~isempty(x), regexp({a(:).name},'run[0-9][0-9][0-9]')).*cellfun( @(x) isempty(x), regexp({a(:).name},'run000'));
            o      = sum(o);
        end
        function fwhm   = vM2FWHM(kappa)
            %fwhm = vM2FWHM(amp,centerX,kappa,offset)
            %transforms a given vonMises function's kappa parameter to FWHM. Kappa
            %parameter has no intuition, however FWHM is easily understandable.
            amp     = 1;
            centerX = 0;
            offset  = 0;
            
            X           = linspace(-180,180,100000);%degrees
            Y           = Tuning.VonMises(X,amp,kappa,centerX,offset);%requires degrees, converts to rads inside.
            
            half_height = [max(Y)-min(Y)]./2+min(Y);%
            d           = abs(Y - half_height);
            
            [~,i]       =  min(d(X > centerX));
            i = i+sum(X < centerX);
            [~,i2]      =  min(d(X < centerX));
            
            % plot(abs(Y-half_height));
            
            fwhm        =  abs(diff(X([i i2])));
        end
    end
end
