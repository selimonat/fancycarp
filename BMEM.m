classdef BMEM < handle
    %creates a Bayesian Magnet Effect Model object
    properties (Constant)        
        options          = optimset('Display','off','maxfunevals',10000,'tolX',10^-16,'tolfun',10^-16,'MaxIter',10000,'Algorithm','interior-point');        
    end
    properties
        x                = linspace(-135,180,8);
        phase            = 3;%Which experimental phase
        prior_function   = 'vonmises'
        param_kappa_csp        = 1;        
        param_kappa_prior      = 1;
        param_kappa_perceptual = 2;
        param_mu               = 0;
        ratings          = [];
        subjects         = [];
        csp              = [];
        p_f_given_csp1   = [];
        f_given_csp1     = [];
        f                = [];
        %p_given_f        = [];%p(p|f); perceptual likelihood
        model
        group
        current_subject = 1;
        current_data;
        current_fit
        current_model 
        fit_quality_NonExpVar   %non explained variance.
        fit_quality_r   = [];     
    end
    properties (Hidden)
        
        fit_estimates   = [];
        fit_mse         = [];
        fit_exitFlag    = [];
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
                for ns = 1:length(G.subject)
                    if ~isempty(G.subject{ns}.ratings(nrun).y)
                        out                        = G.subject{ns}.ratings(nrun).y_mean;%extract it
                        out                        = (out-min(out));
                        bmem.ratings(nrun).y(ns,:) = out./sum(out);
%                         bmem.ratings(nrun).y(ns,:) = out;
                    end
                end                
            end            
        end
    end
    %% methods to run different models across all subjects
    methods 
        function fit_ss(self,nfuns)            
            %Will run different models specified in NFUNS and return the
            %unexplained variance in EXPVARIANCE. All models 
            
            for nfun = nfuns(:)'
                self.current_model = nfun;
                for run = 3
                    for ns = 1:length(self.subjects)          
                        fprintf('Fun:%i, Run:%i, Sub:%i\n',nfun,run,ns);
                        self.current_data   = self.ratings(run).y(ns,:);
                        if ~any(isnan(self.current_data))
                            self.find_params;
                            %
                            self.fit_quality_NonExpVar(ns,run,nfun) = self.fit_mse;
                            self.fit_quality_r(ns,run,nfun)         = corr2(self.current_data(:),self.current_fit(:));
                            %
                            self.plot_subject;pause                            
                        end
                    end
                end
            end
        end
        
        function fitBMEM_ConstantPrior(self)
            %This function is an expansion of fitBMEM to the group level
            %using the same prior function for everybody.
            self.phase = 3;
            run = self.phase;
            param_mask = [1 0 0];%keep prior constant, last param is prior.
            kappa_counter = 0;
            for param_kappa_prior      = 0%linspace(.0001,25,50)                
                kappa_counter          = kappa_counter + 1;
                self.param_kappa_prior = param_kappa_prior;
                for ns = 1:length(self.subjects)
                    fprintf('kappa :%3.4g, Subject:%i\n',param_kappa_prior,ns);
                    self.current_subject                        = ns;
                    self.find_params(param_mask);
                    data                                        = self.group.ratings(run).y_mean(ns,1:8);
                    fit                                         = self.p_f_given_csp1;
                    self.fit_quality_NonExpVar(ns,run,10)       = nanmean((data(:) -fit(:)).^2);%
                    self.fit_quality_r(ns,run,10+kappa_counter) = corr2(data(:),fit(:));
                    self.plot_subject;
                end                
                plot(squeeze(nanmean(self.fit_quality_r(:,3,10:end))));
                drawnow;
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
            if strcmp(self.prior_function,'vonmises')                
                dummy                 = self.VonMises(self.param_kappa_prior,location);            
                dummy                 = dummy + self.VonMises(self.param_kappa_prior,location2);            
            else
                dummy                 = ones(1,8);
            end
            out                       = dummy./sum(dummy(:));
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
        function find_params(self,method,varargin)
            %Matlab's gradient descent to find free parameters optimizing
            %the self.cost_function. Before calling the cost function here
            %we organize initial values and lower/upper bounds.
            if nargin == 2
                free_vector = varargin{1};
            else
                free_vector = [1 1 1];
            end
            %default upper and lower bounds.
            L       = [0  0  0];
            U       = [25 25 25];
            initial = [1 1 1];            
            %            
            L       = L.*free_vector+initial.*(~free_vector);
            U       = U.*free_vector+initial.*(~free_vector);            
            %                                    
            dummy  = @(params) self.cost_function_ss(params);
            
            [self.fit_estimates, self.fit_mse, self.fit_exitFlag]  = ...
                fmincon(dummy, initial, [],[],[],[],L,U,[],self.options);
            
            %overwrite the current parameters depending on the model
             if self.current_model    == 1                                
                self.param_kappa_csp        = self.fit_estimates(1);
                self.param_kappa_perceptual = self.fit_estimates(2);
                self.param_kappa_prior      = self.fit_estimates(3);
                
            elseif self.current_model == 2
                
                self.param_kappa_csp        = self.fit_estimates(1);
                self.param_mu               = self.fit_estimates(2);
                
            elseif self.current_model == 3
                
                self.param_kappa_csp        = self.fit_estimates(1);
             end
        end
        function plot_subject(self)
            %plots active generalization, prior, and BMEM feargen profiles
            %together with the real data.
            clf
            subplot(2,2,1);
            bar(self.csp1_given_face);
            title('Active Generalization');
            %
            subplot(2,2,2);
            bar(self.f);
            title('Prior');
            %
            subplot(2,2,3);
            bar(self.current_fit);
            title('Model');
            %
            subplot(2,2,4);
            bar(self.current_data);
            title('Ratings Data');
            %
            EqualizeSubPlotYlim(gcf);
            supertitle(sprintf('sub:%i, CS+:%i',self.current_subject,self.csp(self.current_subject)),1);            
        end
        function out = cost_function_ss(self,params)
            %out = cost_function_ss(self,params); Single subject cost
            %function for different models specified in self.current_model.
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
                self.param_mu               = params(2);
                self.current_fit            = self.VonMises(self.param_kappa_csp,self.param_mu);                
            elseif self.current_model == 3
                self.param_kappa_csp        = params(1);                
                self.current_fit            = self.VonMises(self.param_kappa_csp,0);
            end            
            %MSE
            out                             = mean(( self.current_fit(:) - self.current_data(:) ).^2);            
        end        
    end    
end










