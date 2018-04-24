classdef Group < Project
    properties (Hidden,Constant)
        mean_correction = 0;%decides if mean correction should be applied
        align_tunings   = 1;%should ratings be aligned to CS+ face
    end
    properties
        subject
        ids
        csps
        table_pmf
        ratings
        total_subjects
        tunings
        fit_results
    end
    
    methods
        
        %%
        function group = Group(subjects)
            c = 0;
            for s = subjects(:)'
                fprintf('subject: %03d\n',s)
                c                = c+1;
                dummy            = Subject(s);
                group.subject{c} = dummy;
                group.ids        = subjects;                
            end
        end
        %%
        function out = get.ratings(self)                        
            %will collect group ratings as a matrix in out(run).y,
            %out(run).x etc.
            for s = 1:self.total_subjects
                for fields = fieldnames(self.subject{s}.ratings)'
                    out.(fields{1})(s,:) = self.subject{s}.ratings.(fields{1})(:);
                end
            end
        end
        %%        
        function out = get.table_pmf(self)            
            %returns pmf parameters as a table
            out = [];
            for s = self.subject                
                out = [out ;s{1}.pmf_parameters];
            end            
        end
        %%        
        function csps = get.csps(self)   
            %returns the csp face for all the group... In the future one
            %could make a get_subject_field function
            for s = 1:length(self.subject)
                csps(s) = self.subject{s}.csp;
            end
        end   
        %% 
        function out = get.total_subjects(self)
            out = length(self.ids);
        end
        %%
        function plot_ratings(self)
            %%
            [y x]  =  GetSubplotNumber(self.total_subjects);
            for ns = 1:self.total_subjects
                subplot(y,x,ns);
                self.feargen_plot(self.ratings.y_mean(1,:,ns));
                box off;
                title(sprintf('s: %d, cs+: %d',self.ids(ns),self.subject{ns}.csp),self.font_style{:}); 
            end
            EqualizeSubPlotYlim(gcf);
            supertitle(self.path_project,1);
        end        
        %%
        function model_ratings(self,run,funtype)
            %will fit to ratings from RUN the FUNTYPE and cache the result
            %in the midlevel folder.            
            T = [];%future tuning object
            filename               = sprintf('%smidlevel/Tunings_Run_%03d_FunType_%03d_N%s.mat',self.path_project,run,funtype,sprintf('%s\b\b',sprintf('%ito',self.ids([1 end]))));
            if exist(filename) == 0
                %create a tuning object and fits FUNTYPE to it.
                T  = Tuning(self.ratings(run));%create a tuning object for the RUN for ratings.                
                T.SingleSubjectFit(funtype);%call fit method from the tuning object                
                save(filename,'T');
            else
                fprintf('Will load the tuning parameters from the cache:\n%s\n',filename);
                load(filename);
            end
            %get the relevant data from the tuning object            
            self.fit_results = T.fit_results;
        end               
        %%
        function feargen_plot(self,data)
            %elementary function to make feargen plots
            h = bar(data,1);            
            self.set_feargen_colors(h,2:9);
            set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'},'xgrid','on',self.font_style{:});
            axis tight;
        end         
        %%
        function ModelSCR(self,run,funtype)
            %create a tuning object and fits FUNTYPE to it.
            self.tunings.scr = Tuning(self.getSCRs(run));%create a tuning object for the RUN for SCRS.
%             self.tunings.scr.SingleSubjectFit(funtype);%call fit method from the tuning object
        end
        function getSCRtunings(self,run,funtype)
            self.ModelSCR(run,funtype);
        end
        
        function [out] = getSCRmeans(self,phase)
            for n = 1:length(self.ids)
                    ind = self.subject{n}.scr.findphase(phase);
                    self.subject{n}.scr.cut(ind);
                    self.subject{n}.scr.run_ledalab;
                    out(n,:) = mean(self.subject{n}.scr.ledalab.mean(1:800,:));
            end
        end
      
       

       
        function plotPMFbars(self)
            means     = reshape(mean(self.pmf.params1(:,1,:),3),2,2);%compute the mean
            stds      = reshape(std(self.pmf.params1(:,1,:),0,3),2,2);
            sem       = stds/sqrt(length(self.ids));
            
            fig=figure;
            [h,e] = barwitherr(sem,means);
            set(gca,'XTickLabel',{'CS+','CS-'})
            set(e,'LineWidth',1.5)
            set(h(1), 'FaceColor','r')
            set(h(2), 'FaceColor',[143/255 0 0 ])
            ylim([20 80])
            ylabel('threshold \alpha (degrees)')
            legend('before','after','orientation','horizontal','location','southoutside')
        end
        function [out labels] = parameterMat(self)
            labels = {'csp_before_alpha' 'csp_after_alpha' 'csn_before_alpha' 'csn_after_alpha' ...
                      'csp_before_beta' 'csp_after_beta' 'csn_before_beta' 'csn_after_beta' ...                     
                      'csp_improvmt' 'csn_improvmnt' ...
                      'csp_imprvmtn_cted' ...
                      'rating_cond' ... 
                      'rating_test' ... 
                      'SI'...
                      'SCR ampl'};
            out = [self.pmf.csp_before_alpha,...
                   self.pmf.csp_after_alpha,...              
                   self.pmf.csn_before_alpha,...
                   self.pmf.csn_after_alpha,...   
                   self.pmf.csp_before_beta,...
                   self.pmf.csp_after_beta,...              
                   self.pmf.csn_before_beta,...
                   self.pmf.csn_after_beta,... 
                   self.pmf.csp_before_alpha - self.pmf.csp_after_alpha,...
                   self.pmf.csn_before_alpha - self.pmf.csn_after_alpha,...                   
                   (self.pmf.csp_before_alpha-self.pmf.csp_after_alpha)-(self.pmf.csn_before_alpha-self.pmf.csn_after_alpha),...
                   self.sigma_cond,...
                   self.sigma_test,...
                   self.SI,...
                   self.SCR_ampl];
            
               try
               if strcmp(self.tunings.rate{3}.singlesubject{1}.funname,'vonmisses_mobile')
                   
                   for s = 1:size(out,1)
                       out(s,15) = self.tunings.rate{3}.singlesubject{s}.Est(3);
                       out(s,16) = self.tunings.rate{4}.singlesubject{s}.Est(3);
                   end
                   out(:,14) = out(:,13) - out(:,12);
                   labels = [labels(1:14) { 'mu_cond' 'mu_test'}];
               end
               end
        end        
        function PlotRatingFit(self,subject)
            
            
            if ~isempty(self.tunings.rate)
               
                i    =  find(self.ids == subject);
                ave  = mean(reshape(self.tunings.rate{3}.y(i,:),2,8));                                             
                x    = mean(reshape(self.tunings.rate{3}.x(i,:),2,8));
                x_HD = linspace(min(x),max(x),1000);
                h    = figure(100);clf
                subplot(1,2,1)
                title(sprintf('Sub: %i, Likelihood: %03g (p = %5.5g)',subject,self.tunings.rate{3}.singlesubject{i}.Likelihood,self.tunings.rate{3}.singlesubject{i}.pval));
                hold on;               
                plot(x_HD,self.tunings.rate{3}.singlesubject{i}.fitfun(x_HD,self.tunings.rate{3}.singlesubject{i}.Est),'ro','linewidth',3);
                plot(x,ave, 'b','linewidth', 3);
                ylabel('Cond')
                drawnow;
                grid on;
                
                subplot(1,2,2)
                ave  = mean(reshape(self.tunings.rate{4}.y(i,:),2,8));  
                title(sprintf('CSP: %i, Likelihood: %03g (p = %5.5g)',self.subject{i}.csp,self.tunings.rate{4}.singlesubject{i}.Likelihood,self.tunings.rate{4}.singlesubject{i}.pval));
                hold on;
                plot(x_HD,self.tunings.rate{4}.singlesubject{i}.fitfun(x_HD,self.tunings.rate{4}.singlesubject{i}.Est),'ro','linewidth',3);
               
                 plot(x,ave, 'b','linewidth', 3);
                ylabel('Test')
                EqualizeSubPlotYlim(h);               
                drawnow;
                grid on;
                pause
            else
                fprintf('No tuning object found here yet...\n');
            end
        end                                
        %%       
        function [rating] = PlotRatings(self,runs)            
            hvfigure;
            trun = length(runs);
            crun = 0;
            for run = runs(:)'%for each run make a subplot column
                crun    = crun + 1;
                %
                subplot(2,trun,crun);                
                rating  = self.Ratings(run);%collect group ratings                
                imagesc(rating.y,[0 10]);thincolorbar('vert');%single subject data
                set(gca,'xticklabel',{'CS+' 'CS-'},'xtick',[4 8],'fontsize',20,'yticklabel',{''});
                colormap hot
                %
                subplot(2,trun,crun+trun);
                [y x] = hist(rating.y);
                y     = y./repmat(sum(y),size(y,1),1)*100;%make it a percentage
                imagesc(rating.x(1,:),x,y,[0 75]);axis xy;
                thincolorbar('vert');
                hold on                
                h     = errorbar(mean(rating.x),mean(rating.y),std(rating.y),'g-');
                axis xy;
                set(gca,'xticklabel',{'CS+' 'CS-'},'xtick',[0 180]);                
                hold off;
            end
        end
        function PlotRatingResults(self)%plots conditioning and test, in the usual bar colors. With GroupFit Gauss/Mises Curve visible
            %%
            f=figure;
            subplot(1,2,1);
            h = bar(unique(self.tunings.rate{3}.x(1,:)),self.tunings.rate{3}.y_mean);SetFearGenBarColors(h);
            hold on;
            errorbar(unique(self.tunings.rate{3}.x(1,:)),self.tunings.rate{3}.y_mean,self.tunings.rate{3}.y_std./sqrt(length(self.ids)),'k.','LineWidth',2);
            xlim([-160 200]);
            box off
            set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
            x = linspace(self.tunings.rate{3}.x(1,1),self.tunings.rate{3}.x(1,end),100);
            plot(x ,  self.tunings.rate{3}.groupfit.fitfun( x,self.tunings.rate{3}.groupfit.Est(:,1:end-1)) ,'k--','linewidth',2);
%             plot(x ,  self.tunings.rate{3}.singlesubject{1}.fitfun( x,mean(self.tunings.rate{3}.params(:,1:end-1))) ,'k--','linewidth',1);
            hold off
%             set(gca,'fontsize',14);
            axis square
            t=title('Conditioning');set(t,'FontSize',14);
            %
            subplot(1,2,2);
            h = bar(unique(self.tunings.rate{4}.x(1,:)),self.tunings.rate{4}.y_mean);SetFearGenBarColors(h);hold on;
            errorbar(unique(self.tunings.rate{4}.x(1,:)),self.tunings.rate{4}.y_mean,self.tunings.rate{4}.y_std./sqrt(length(self.ids)),'k.','LineWidth',2);
            EqualizeSubPlotYlim(gcf);
            box off
            xlim([-160 200]);
            set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
            x = linspace(self.tunings.rate{4}.x(1,1),self.tunings.rate{4}.x(1,end),100);
%             plot(x ,  self.tunings.rate{4}.singlesubject{1}.fitfun( x,mean(self.tunings.rate{4}.params(:,1:end-1))) ,'k','linewidth',1);
            plot(x ,  self.tunings.rate{4}.groupfit.fitfun( x,self.tunings.rate{4}.groupfit.Est(:,1:end-1)) ,'k','linewidth',2);
            x = linspace(self.tunings.rate{3}.x(1,1),self.tunings.rate{3}.x(1,end),100);
% % %            plot(x ,  self.tunings.rate{3}.singlesubject{1}.fitfun( x,mean(self.tunings.rate{3}.params(:,1:end-1))) ,'k--','linewidth',1);
%             plot(x , self.tunings.rate{3}.groupfit.fitfun( x,self.tunings.rate{3}.groupfit.Est(:,1:end-1)) ,'k--','linewidth',2);
%             set(gca,'fontsize',14);
            axis square
            t=title('Test');set(t,'FontSize',14);
            annotation(f,'textbox',[0.78 0.65 0.1 0.1],'String',['SI = ' num2str(nanmean(self.SI))],'FitBoxToText','off','LineStyle','none');
            hold off
        end
        %%
        function [scr]    = getSCRs(self,run)
            %will collect the ratings from single subjects
            scr.y = [];
            scr.x = [];
            scr.ids = [];
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                    dummy = self.subject{s}.GetSubSCR(run);
                    if ~isempty(dummy)
                        scr.y   = [scr.y; dummy.y];
                        scr.x   = [scr.x; dummy.x];
                        scr.ids = [scr.ids; self.ids(s)];
                    end
                end
            end
        end
        
        
    end
    methods %(mri))
        function Fit2ndlevel(self,nrun,modelnum,namestring,varargin)
            %% 2ndlevel 8conds ANOVA
            clear2ndlevel = 1;
            versiontag = 0;
            prefix = 's6_wCAT_';
            foldersuffix = '';
            
            if nargin == 5 %one varargin is given
                versiontag = varargin{1};
            elseif nargin == 6
                versiontag = varargin{1};
                foldersuffix = varargin{2};
            end
            
            if strcmp(namestring,'CSdiff')
                cons2collect   = 1;
            elseif strcmp(namestring,'8conds')
                %go through all conditions, i.e. 8 (B/T) or 2 (C)
                if ismember(nrun,[1 3])
                    cons2collect = 2:9;
                elseif nrun == 2
                    cons2collect = 2:3;
                end
                
            elseif strcmp(namestring,'8conds_rate')
                switch nrun
                    case 1
                        cons2collect = 1:8;%WIP
                    case 2
                        cons2collect = 2:4;
                    case 3
                        cons2collect = 2:4;
                end
            elseif strcmp(namestring,'VMdVM')
                cons2collect = 2:3;
            end
            
            dependencies = 1;
            unequalvar   = 1;
                        start = tic;
            fprintf('Starting 2nd Level for model %02d, run %02d, named %s, with foldersuffix ''%s''...\n',modelnum,nrun,namestring,foldersuffix);
            
            path2ndlevel = fullfile(self.path_second_level,sprintf('model_%02d_chrf_%01d_%s_%s%s',modelnum,versiontag,namestring,self.nrun2phase{nrun},foldersuffix));
            if exist(path2ndlevel) && (clear2ndlevel==1);
                system(sprintf('rm -fr %s*',strrep(path2ndlevel,'//','/')));
            end%this is AG style.
            if ~exist(path2ndlevel)
                mkdir(path2ndlevel);
            end
              clear matlabbatch
              
                load(self.subject{1}.path_spmmat(nrun,modelnum));
            % collect all subjects' con images for every cond
            c = 0;
            for ncon = cons2collect(:)'
                c = c+1;
                clear files
                fprintf('Getting con_%04d, %s from sub ',ncon, SPM.xCon(ncon).name)
                for ns = 1:numel(self.ids)
                    files(ns,:) = strrep(self.subject{ns}.path_con(nrun,modelnum,prefix,ncon),'sub004',sprintf('sub%03d',self.ids(ns)));
                    fprintf('%d..',self.ids(ns))
                end
                fprintf('.done. (N=%02d).\n',ns)
                matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(c).scans = cellstr(files); %one cond at a time, but all subs
            end
            
            % specify rest of the model
            matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(path2ndlevel);
            matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = dependencies;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = unequalvar;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = -Inf;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {[self.path_groupmeans '/ave_wCAT_s3_ss_data.nii']};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
          
            
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {[path2ndlevel filesep 'SPM.mat']};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            spm_jobman('run',matlabbatch);
            done = toc(start);
            fprintf('\n\nDone estimating 2nd Level for HRF model %d, called %s %s, version %d, N = %02d subs in %05.2f mins. (Simple one-way anova)\n',modelnum,namestring,foldersuffix,versiontag,length(self.ids),done./60)
            fprintf('Output folder is: %s\n',path2ndlevel)

        end
        function Con2ndLevel(self,nrun,modelnum,namestring,versiontag,foldersuffix)
            nF = 0;
            nT = 0;
            n  = 0;
            
            path_spmmat = fullfile(self.path_second_level,sprintf('model_%02d_chrf_%01d_%s_%s%s',modelnum,versiontag,namestring,self.nrun2phase{nrun},foldersuffix),'SPM.mat');
         
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(path_spmmat);
            
            if strcmp(namestring,'CSdiff')
                n = n + 1; nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = 1;
                 matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                 
                 matlabbatch{1}.spm.stats.con.delete = 0;
                 
            elseif strcmp(namestring,'8conds')
                if ismember(nrun,[1 3])
                    nconds = 8;
                else
                    nconds = 2;
                end
                
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_F_8conds';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'main_allconds';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                if ismember(nrun,[1 3])
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 0 0 1 0 0 0 -1];
                else
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 -1];
                end
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                            
                if ismember(nrun,[1 3])
                    [VM, dVM] = self.compute_VM(-135:45:180,1,1,.001);
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = VM;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'VMtuning';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = dVM;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'dVMtuning';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [-repmat(1/7,1,3) 1 -repmat(1/7,1,3)];
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'CSP>rest';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                end
                
            elseif strcmp(namestring,'VMdVM')
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'pp';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(2);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nF = nF + 1;
                vec = eye(2); vec(logical(eye(2))) = [1 -1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'pn';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = vec;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'VM';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 0];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'dVM';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            end
            
            if nT > 0
                for tc = 1:nT
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.name =      [matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.name '_neg'];
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.weights =   -matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.weights;
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.sessrep =   'none';
                    nT = nT + 1;
                end
            end
            
           matlabbatch{1}.spm.stats.con.delete = 1;
                        
            spm_jobman('run',matlabbatch);

            ntotal = nT + nF;            
            fprintf('Done creating %d 2ndlevel contrasts (%d F, %d T) for model %d, run %s, modelname %s %s.\n',ntotal,nF,nT,modelnum,self.nrun2phase{nrun},namestring,foldersuffix)
            if nF > 0
                for nnF = 1:nF
                    disp(['(F) ' matlabbatch{1}.spm.stats.con.consess{nnF}.fcon.name])
                end
            end
            for nnT = 1:nT
                disp(['(T) ' matlabbatch{1}.spm.stats.con.consess{nF+nnT}.tcon.name])
            end
        end
        function Fit2ndlevel_FIR(self,nrun,modelnum,namestring,varargin)
            %% 2ndlevel 8conds ANOVA
            clear2ndlevel = 1;
            versiontag = 0;
            prefix = 's6_wCAT_';
            foldersuffix = '';
            bins2take = 1:14;
            
            if nargin == 5 %one varargin is given
                versiontag = varargin{1};
            elseif nargin == 6
                versiontag = varargin{1};
                foldersuffix = varargin{2};
            end
            
          
            
            
            
            
            
            
            
            
            
            
            if strcmp(namestring,'CSdiff')
                cons2collect   = 1;
            elseif strcmp(namestring,'8conds')
                %go through all conditions, i.e. 8 (B/T) or 2 (C)
                if ismember(nrun,[1 3])
                    cons2collect = 2:9;
                elseif nrun == 2
                    cons2collect = 2:3;
                end
                
            elseif strcmp(namestring,'8conds_rate')
                switch nrun
                    case 1
                        cons2collect = 1:8;%WIP
                    case 2
                        cons2collect = 2:4;
                    case 3
                        cons2collect = 2:4;
                end
            elseif strcmp(namestring,'VMdVM')
                cons2collect = 2:3;
            end
            
            dependencies = 1;
            unequalvar   = 1;
            start = tic;
            fprintf('Starting 2nd Level for model %02d, run %02d, named %s, with foldersuffix ''%s''...\n',modelnum,nrun,namestring,foldersuffix);
            
            path2ndlevel = fullfile(self.path_second_level,sprintf('model_%02d_chrf_%01d_%s_%s%s',modelnum,versiontag,namestring,self.nrun2phase{nrun},foldersuffix));
            if exist(path2ndlevel) && (clear2ndlevel==1);
                system(sprintf('rm -fr %s*',strrep(path2ndlevel,'//','/')));
            end%this is AG style.
            if ~exist(path2ndlevel)
                mkdir(path2ndlevel);
            end
              clear matlabbatch
              
                load(self.subject{1}.path_spmmat(nrun,modelnum));
            % collect all subjects' con images for every cond
            c = 0;
            for ncon = cons2collect(:)'
                c = c+1;
                clear files
                fprintf('Getting con_%04d, %s from sub ',ncon, SPM.xCon(ncon).name)
                for ns = 1:numel(self.ids)
                    files(ns,:) = strrep(self.subject{ns}.path_con(nrun,modelnum,prefix,ncon),'sub004',sprintf('sub%03d',self.ids(ns)));
                    fprintf('%d..',self.ids(ns))
                end
                fprintf('.done. (N=%02d).\n',ns)
                matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(c).scans = cellstr(files); %one cond at a time, but all subs
            end
            
            % specify rest of the model
            matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(path2ndlevel);
            matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = dependencies;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = unequalvar;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = -Inf;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {[self.path_groupmeans '/ave_wCAT_s3_ss_data.nii']};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
          
            
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {[path2ndlevel filesep 'SPM.mat']};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            spm_jobman('run',matlabbatch);
            done = toc(start);
            fprintf('\n\nDone estimating 2nd Level for HRF model %d, called %s %s, version %d, N = %02d subs in %05.2f mins. (Simple one-way anova)\n',modelnum,namestring,foldersuffix,versiontag,length(self.ids),done./60)
            fprintf('Output folder is: %s\n',path2ndlevel)
        end
    end
end
