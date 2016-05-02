classdef BMEM < handle
    %creates a Bayesian Magnet Effect Model object
    properties (Constant)        
        options          = optimset('Display','off','maxfunevals',10000,'tolX',10^-9,'tolfun',10^-9,'MaxIter',10000,'Algorithm','interior-point');        
    end
    properties (Hidden)
        params;
        gridsize      = 20;%resolution per parameter for initial estimation.
        current_subject         = 1;
        current_data;
        current_fit;
        current_model;
        current_fit_estimates   = [];
        current_fit_mse         = [];
        current_fit_exitFlag    = [];
        current_fit_r           = [];
        current_initial         = [];
        visualization           = 1;
        LB%lower and
        UB%upper bounds for parameter estimates
        phase_names             = {'Baseline' 'Conditioning' 'Test'};
        fun_names               = {'BMEM' 'vM' 'Gau'};
    end
    properties
        x                      = linspace(-135,180,8);
        prior_function         = 'vonmises';
        param_kappa_csp        = 1;        
        param_kappa_prior      = 1;
        param_kappa_perceptual = 2;
        param_mu               = 0;
        ratings                = [];
        subjects               = [];
        csp                    = [];
        p_f_given_csp1         = [];
        f_given_csp1           = [];
        f                      = [];        
        model;
        fit_quality_NonExpVar;   %non explained variance.
        fit_quality_r          = [];     
    end
 
    methods
        function bmem = BMEM(G)
            %creates a Bayesian Magnet Effect Model object
            %G is simply a Group object
            %%
            bmem.subjects       = G.ids;
            bmem.csp            = G.csps;
            %collect the ratings, this is just for convenience and
            %performance. storing it indepedent of the group object makes
            %accessing the ratings much faster. The main point is that
            %ratings are collected as they are stored in the GROUP object.
            %So if you want to apply a different normalization scheme,
            %please do that in the Group object and re-run the BMEM object.
            for nrun = 1:3
                bmem.ratings(nrun).y = NaN(length(G.subject),8);
                for ns = 1:length(G.subject);
                    if ~isempty(G.subject{ns}.ratings(nrun).y);
                        out                        = G.subject{ns}.ratings(nrun).y_mean;%extract it
                        out                        = (out-min(out));
                        bmem.ratings(nrun).y(ns,:) = out./sum(out);
%                         bmem.ratings(nrun).y(ns,:) = out;
                    end
                end                
            end
            %% reset figures if needed.
            if bmem.visualization
                figure(9);clf;
                figure(10);clf;                
            end
        end
    end
    %% methods to run different models across all subjects
    methods 
        function fit_ss(self,nfuns)            
            %Will run models specified in NFUNS and return the fit quality.
            %Relies on self.find_params.
            % NFUN:
            % 1: BMEM model            
            % 2: vM model
            % 3: Gaussian Model                                    
            
            for nfun = nfuns(:)'
                self.current_model = nfun;
                for run = 1:3
                    for ns = 1:length(self.subjects)
                        fprintf('Fun:%i, Run:%i, Sub:%i\n',nfun,run,ns);
                        self.current_data   = self.ratings(run).y(ns,:);
                        self.current_subject = ns;
                        if ~any(isnan(self.current_data))
                            self.SetBoundary;
                            self.find_params([1 1 1]);                            
                            %
                            self.fit_quality_NonExpVar(ns,run,nfun) = self.current_fit_mse;
                            self.fit_quality_r(ns,run,nfun)         = self.current_fit_r;                            
                            %                            
                            if self.visualization
                                self.plot_subject;                               
                            end
                            %pause                            
                        else
                            self.fit_quality_NonExpVar(ns,run,nfun) = NaN;
                            self.fit_quality_r(ns,run,nfun)         = NaN;
                        end                        
                    end
                    if self.visualization
                        self.plot_group;
                    end
                end
            end
        end
        function fit_constantprior(self)
            %This function is an expansion of fit_ss to the group level
            %using the same prior function for everybody.            
            self.current_model = 1;                        
            for run = 1:3
                kappa_counter      = 0;
                for param_kappa_prior = linspace(.0001,25,16);
                    fprintf('kappa :%3.4g\n',param_kappa_prior);
                    kappa_counter          = kappa_counter +1;
                    self.param_kappa_prior = param_kappa_prior;
                    for ns = 1:length(self.subjects)
                        self.current_subject                                   = ns;
                        self.current_data                                      = self.ratings(run).y(ns,:);
                        if ~any(isnan(self.current_data))
                            %
                            self.SetBoundary;
                            self.LB(3)                                         = param_kappa_prior;
                            self.UB(3)                                         = param_kappa_prior;
                            self.find_params([ 1 1 0]);%keep prior constant
                            %
                            self.fit_quality_NonExpVar(ns,run,3+kappa_counter) = self.current_fit_mse;
                            self.fit_quality_r(ns,run,3+kappa_counter)         = self.current_fit_r;                                                        
                        else
                            self.fit_quality_NonExpVar(ns,run,4+kappa_counter-1) = NaN;
                            self.fit_quality_r(ns,run,4+kappa_counter-1)         = NaN;
                        end
                    end
                    if self.visualization
                        self.plot_subject;
                        self.plot_group;
                    end
                end
            end
        end
    end
    
    %% Methods for the probabilistic machinery.
    methods
        function out = p_given_f(self,location)
            %p(p|f) the perceptual likelihood function
            out = [];
            c = 0;
            for l = location(:)'
                c        = c + 1;
                out(:,c) = self.VonMises(self.param_kappa_perceptual,l);
            end
        end
        function out = csp1_given_face(self)
            %p(CS+ = 1 | face) the active generalization component.
            out = self.VonMises(self.param_kappa_csp,0);            
        end
        function out = get.f(self)
            %p(f) returns prior on faces            
            location               = (3.5-self.csp(self.current_subject))*45;
            location2              = (7.5-self.csp(self.current_subject))*45;
            if strcmp(self.prior_function,'double_category')                
                dummy              = self.VonMises(self.param_kappa_prior,location);            
                dummy              = dummy + self.VonMises(self.param_kappa_prior,location2);            
            elseif strcmp(self.prior_function,'male')
                dummy              = self.VonMises(self.param_kappa_prior,location);
            elseif strcmp(self.prior_function,'female')
                dummy              = self.VonMises(self.param_kappa_prior,location2);
            else
                dummy              = ones(1,8);
            end
            out                    = dummy./sum(dummy(:));
        end
        function out = csp1_f(self)
            %p(csp = 1,face)
            out = self.csp1_given_face.*self.f;
        end
        function out = get.f_given_csp1(self)
            %p(face|csp = 1)
            out = self.csp1_f;
            out = out./sum(out,2);
        end
        function out = get.p_f_given_csp1(self)
            %p(p,f | csp = 1);
            out = sum(self.p_given_f(self.x)*diag(self.f_given_csp1),2);
        end
        
        function out = VonMises(self,kappa,mu)
            %VM probability density function
            out    =  exp(kappa*cos(deg2rad(self.x-mu)));
            out    =  out./sum(out(:));
        end
    end
        
    %% Basic Parameter optimization methods
    methods
        function find_params(self,free_vector)
            % Matlab's gradient descent to find free parameters optimizing
            % the self.cost_function. Before calling the cost function here
            % we organize initial values and lower/upper bounds. FREE_VECTOR
            % distinguishes variables from constants.
            %
            % Will find the best fit for the model specified in
            % CURRENT_MODEL, and update CURRENT_FIT_R, CURRENT_FIT_MSE,
            % CURRENT_                                        
            tparam               = length(self.LB);                      
            
            % start with optimization            
            dummy  = @(params) self.cost_function_ss(params);
            % roughtly estimate initial position
            self.current_initial = self.RoughEstimator(self.x,self.current_data,dummy,self.LB,self.UB);
            %now limit the L - U to the same value if the parameter is
            %supposed to stay constant.
            self.LB  = self.LB(1:tparam).*free_vector(1:tparam) + self.current_initial(1:tparam).*(~free_vector(1:tparam));
            self.UB  = self.UB(1:tparam).*free_vector(1:tparam) + self.current_initial(1:tparam).*(~free_vector(1:tparam));                        
            
            [self.current_fit_estimates, self.current_fit_mse, self.current_fit_exitFlag] = ...
                fmincon(dummy, self.current_initial, [],[],[],[],self.LB,self.UB,[],self.options);
            
            %save the correlation
            self.current_fit_r = corr2(self.current_data(:),self.current_fit(:));
            
            %overwrite the current parameters depending on the fit results.
            %We do that so that plot methods can plot something meaningful.
             if self.current_model    == 1                                
                self.param_kappa_csp        = self.current_fit_estimates(1);
                self.param_kappa_perceptual = self.current_fit_estimates(2);
                self.param_kappa_prior      = self.current_fit_estimates(3);                
            elseif self.current_model == 2                
                self.param_kappa_csp        = self.current_fit_estimates(1);
                self.param_mu               = self.current_fit_estimates(2);                
            elseif self.current_model == 3                
                self.param_kappa_csp        = self.current_fit_estimates(1);
             end
            %
        end
        function SetBoundary(self)
            % set upper and lower bounds for different models
            if self.current_model        == 1%BMEM
                self.LB        = [0    0   0];%[kappa_csp kappa_perceptual kappa_prior]
                self.UB        = [25   25  25];
            elseif self.current_model    == 2%vM
                self.LB        = [0  -135 ];%[kappa mu]
                self.UB        = [25  180 ];
            elseif self.current_model    == 3%Immobile vM (i.e. Gaussian)
                self.LB        = 0;%[kappa_csp]
                self.UB        = 25;
            end       
        end
        function [params]=RoughEstimator(self,x,y,fun,L,U)
            %will roughly estimate free parameter values bounded between L
            %and U for FUN and x values. Error is computed according to Y.        
            %% get a grid                
            tparam   = length(L);
            for n = 1:tparam
                grid_in{n} = linspace(L(n),U(n),self.gridsize);
            end
            G      = cell(1,tparam);
            [G{:}] = ndgrid(grid_in{:});%generate a grid
            G      = cat(tparam+1,G{:});%cell to matrix
            G      = permute(G,[tparam+1 1:tparam]);%change dimension orders so that
            G      = reshape(G,[tparam,self.gridsize^tparam])';%we can make Nx3 matrix with reshape
            % compute the error for each parameter combination
            error  = zeros(1,self.gridsize^tparam);
            for npoint = 1:self.gridsize^tparam;
                error(npoint) = fun(G(npoint,:));%residual error
            end
            [m i]  = min(error);
            params = double(G(i,:));
        end
        function plot_group(self)            
            figure(10);clf
            tfun    = size(self.fit_quality_r,3);
            trun    = size(self.fit_quality_r,2);
            for nlab = 1:tfun;try;xlab{nlab} = sprintf('                    %s',self.fun_names{nlab});catch;xlab{nlab}='';end;end
            for run = 1:trun;
                subplot(1,trun,run);
                for nfun = 1:tfun
                    for ns = 1:self.current_subject
                        text(nfun+rand(1).*.6-.5,self.fit_quality_r(ns,run,nfun),mat2str(self.subjects(ns)),'fontsize',9,'color',[(nfun-1)/tfun 0 1-(nfun-1)/tfun]);
                        hold on;
                    end
                end
                box off;
                ylabel('r')
                ylim([-1 1]);
                xlim([0 tfun+1]);
                set(gca,'xgrid','on','xtick',linspace(1,tfun,tfun)-.5,'xticklabel',xlab)
                title(sprintf('Phase: %s',self.phase_names{run}));
            end
            drawnow;
        end
        
        function plot_subject(self)
            %plots active generalization, prior, and BMEM feargen profiles
            %together with the real data.
            figure(9);
            clf            
            subplot(2,3,1);
            bar(self.csp1_given_face);
            title(sprintf('Active Generalization\n(Kappa: (%3.4g) %3.4g)',self.current_initial(1),self.param_kappa_csp ));
            box off;grid on;
            set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'});
            %
            subplot(2,3,2);
            bar(self.f);
            try
                title(sprintf('Prior\n(Kappa: (%3.4g) %3.4g)',self.current_initial(2),self.param_kappa_prior));
            end
            box off;grid on;
            set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'});
            %
            subplot(2,3,3);
            bar(self.p_given_f(0));
            try
                title(sprintf('p(p|f)\n(Kappa: (%3.4g) %3.4g)',self.current_initial(3),self.param_kappa_perceptual));
            end
            box off;grid on;set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'});
            %
            subplot(2,3,4);
            bar(self.current_fit);            
            title(sprintf('Model (#%i)\nr: %3.4g\nParams:%s',self.current_model,self.current_fit_r,sprintf('%2.3g, ',self.current_fit_estimates)));
            box off;grid on;set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'});
            %
            subplot(2,3,5);
            bar(self.current_data);
            title('Ratings Data');
            box off;grid on;set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'});
            %
            EqualizeSubPlotYlim(gcf);
            supertitle(sprintf('sub:%i, CS+:%i',self.subjects(self.current_subject),self.csp(self.current_subject)),1,'fontsize',15);            
        end
        function plot_group_connected_lines(self)
            %
            figure(11);clf;
            tfun   = size(self.fit_quality_r,3);
            tphase = size(self.fit_quality_r,2);
            spsize = GetSubplotNumber(tfun);
            for nlab = 1:tfun;try;xlab{nlab} = sprintf('%s',self.phase_names{nlab}(1));end;end;
            for nfun = 1:tfun
                subplot(spsize,spsize,nfun);
                for ns = 1:size(self.fit_quality_r,1);                
                   PlotTransparentLine([1:tphase]',squeeze(self.fit_quality_r(ns,:,nfun))',.05,[(nfun-1)/tfun 0 1-(nfun-1)/tfun],'linewidth',2)
                end
                %plot also the mean
                hold on;
                plot([1:tphase]',nanmedian(squeeze(self.fit_quality_r(:,:,nfun))),'.-','color',[(nfun-1)/tfun 0 1-(nfun-1)/tfun],'markersize',40,'linewidth',3)
                hold off;   
                ylim([-1 1]);
                set(gca,'xtick',1:tfun,'xticklabel',xlab,'ytick',[-1 0 1],'xgrid','on','fontsize',25);
                ylabel('r')
                box off
                try;title(sprintf('%s',self.fun_names{nfun}));end
            end
            
        end
        function out = cost_function_ss(self,params)
            %out = cost_function_ss(self,params); Single subject cost
            %function for different models specified in self.current_model.
            %First insert the param vector (as used by fmincon) to the
            %object and then evaluate the model
            %
            % 
            
            %%
            if self.current_model == 1
                self.param_kappa_csp        = params(1);
                self.param_kappa_perceptual = params(2);
                self.param_kappa_prior      = params(3);
                self.current_fit            = self.p_f_given_csp1;                
            elseif self.current_model == 2
                self.param_kappa_csp        = params(1);
                self.param_kappa_perceptual = params(2);
                self.param_kappa_prior      = NaN;
                self.current_fit            = self.VonMises(self.param_kappa_csp,self.param_kappa_perceptual);                
            elseif self.current_model == 3
                self.param_kappa_csp        = params(1);                
                self.param_kappa_perceptual = NaN;
                self.param_kappa_prior      = NaN;
                self.current_fit            = self.VonMises(self.param_kappa_csp,0);
            end            
            %MSE
            out                             = mean(( self.current_fit(:) - self.current_data(:) ).^2);                        
        end        
    end    
end










