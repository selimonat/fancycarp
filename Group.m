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
        default_run     = 3;
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
            for nrun = 1:3
                for s = 1:self.total_subjects
                    for fields = fieldnames(self.subject{s}.ratings)'
                        R = mean(self.subject{s}.ratings(self.default_run).(fields{1}),1);
                        if ~isempty(R)
                            out(nrun).(fields{1})(s,:) = R;
                        else
                            out(nrun).(fields{1})(s,:) = NaN;
                        end
                    end
                end
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
                self.feargen_plot(self.ratings.y_mean(ns,:));
                box off;
                title(sprintf('s: %d, cs+: %d',self.ids(ns),self.subject{ns}.csp),self.font_style{:}); 
                axis square;
            end
            EqualizeSubPlotYlim(gcf);
            supertitle(self.path_project,1);
        end        
        %%
        function model_ratings(self,run,funtype)
            %will fit to ratings from RUN the FUNTYPE and cache the result
            %in the midlevel folder.            
            T = [];%future tuning object
            filename = sprintf('%smidlevel/Tunings_Run_%03d_FunType_%03d_N%s.mat',self.path_project,run,funtype,sprintf('%s\b\b',sprintf('%ito',self.ids([1 end]))));
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
        
        
        
    end
end
