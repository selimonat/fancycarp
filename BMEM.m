classdef BMEM < handle
    %creates a Bayesian Magnet Effect Model object
    properties (Constant)        
        options          = optimset('Display','off','maxfunevals',10000,'tolX',10^-16,'tolfun',10^-16,'MaxIter',10000,'Algorithm','interior-point');        
    end
    properties (Hidden)
        params 
        current_subject  = 1;
        current_data;
        current_fit
        current_model 
        current_fit_estimates   = [];
        current_fit_mse         = [];
        current_fit_exitFlag    = [];
        current_fit_r           =[];
        current_initial         = [];
    end
    properties
        x                = linspace(-135,180,8);        
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
        model                
        fit_quality_NonExpVar   %non explained variance.
        fit_quality_r   = [];     
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
                            self.current_initial                    = .1+rand(1,3)*10;%set the initial values
                            self.find_params;                            
                            %
                            self.fit_quality_NonExpVar(ns,run,nfun) = self.current_fit_mse;
                            self.fit_quality_r(ns,run,nfun)         = self.current_fit_r;                            
                            %
                            self.current_fit_r;
%                             self.plot_subject;pause
                        else
                            self.fit_quality_NonExpVar(ns,run,nfun) = NaN;
                            self.fit_quality_r(ns,run,nfun)         = NaN;
                        end
                    end
                end
            end
        end
        function fit_constantprior(self)
            %This function is an expansion of fitBMEM to the group level
            %using the same prior function for everybody.
            
            self.current_model = 1;            
            param_mask         = [ 1 1 0];%[csp perceptual prior]
            for run = 2
                kappa_counter      = 0;
                for param_kappa_prior = linspace(.0001,25,100)                    
                    fprintf('kappa :%3.4g\n',param_kappa_prior);
                    kappa_counter          = kappa_counter +1;
                    self.param_kappa_prior = param_kappa_prior;
                    for ns = 1:length(self.subjects)
                        self.current_subject                                     = ns;
                        self.current_data                                        = self.ratings(run).y(ns,:);
                        if ~any(isnan(self.current_data))
                            %                            
                            self.current_initial                                 = [.1+rand(1,3)*10 self.param_kappa_prior];
                            self.find_params(param_mask);
                            %
                            self.fit_quality_NonExpVar(ns,run,4+kappa_counter-1) = self.current_fit_mse;
                            self.fit_quality_r(ns,run,4+kappa_counter-1)         = self.current_fit_r;
                            
%                             self.plot_subject;pause
                        else
                            self.fit_quality_NonExpVar(ns,run,4+kappa_counter-1) = NaN;
                            self.fit_quality_r(ns,run,4+kappa_counter-1)         = NaN;
                        end
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
            if strcmp(self.prior_function,'vonmises')                
                dummy              = self.VonMises(self.param_kappa_prior,location);            
                dummy              = dummy + self.VonMises(self.param_kappa_prior,location2);            
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
        function find_params(self,varargin)
            %Matlab's gradient descent to find free parameters optimizing
            %the self.cost_function. Before calling the cost function here
            %we organize initial values and lower/upper bounds.
            
            % organize initial values and bound constraints.
            if nargin == 2
                free_vector = varargin{1};
            else
                free_vector = [1 1 1];
            end
                        
            % set upper and lower bounds for different models
            if self.current_model    == 1%BMEM
                L       = [0  2  0];
                U       = [25 2 25];                
                tparam = 3;
            elseif self.current_model    == 2%vM
                L       = [0  0  ];
                U       = [25 180];                
                tparam  = 2;
            elseif self.current_model    == 3%Gaussian
                L       = 0;
                U       = 25;                
                tparam  = 1;            
            end       
            L                    = L(1:tparam).*free_vector(1:tparam)+self.current_initial(1:tparam).*(~free_vector(1:tparam));
            U                    = U(1:tparam).*free_vector(1:tparam)+self.current_initial(1:tparam).*(~free_vector(1:tparam));
            self.current_initial = self.current_initial(1:tparam);
            
            % start with optimization            
            dummy  = @(params) self.cost_function_ss(params);            
            [self.current_fit_estimates, self.current_fit_mse, self.current_fit_exitFlag]  = ...
                fmincon(dummy, self.current_initial, [],[],[],[],L,U,[],self.options);
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
        end
        function plot_subject(self)
            %plots active generalization, prior, and BMEM feargen profiles
            %together with the real data.
            clf
            subplot(2,3,1);
            bar(self.csp1_given_face);
            title(sprintf('Active Generalization\n(Kappa: (%3.4g) %3.4g)',self.current_initial(1),self.param_kappa_csp ));
            box off;grid on
            set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'})
            %
            subplot(2,3,2);
            bar(self.f);
            try
                title(sprintf('Prior\n(Kappa: %3.4g %3.4g)',self.current_initial(2),self.param_kappa_prior));
            end
            box off;grid on
            set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'})
            %
            subplot(2,3,3);
            bar(self.p_given_f(0))
            try
                title(sprintf('p(p|f)\n(Kappa: %3.4g %3.4g)',self.current_initial(3),self.param_kappa_perceptual));
            end
            box off;grid on;set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'})
            %
            subplot(2,3,4);
            bar(self.current_fit);            
            title(sprintf('Model\nr: %3.4g\n%s',self.current_fit_r,sprintf('%2.3g, ',self.current_fit_estimates)));
            box off;grid on;set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'})
            %
            subplot(2,3,5);
            bar(self.current_data);
            title('Ratings Data');
            box off;grid on;set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'})
            %
            EqualizeSubPlotYlim(gcf);
            supertitle(sprintf('sub:%i, CS+:%i',self.current_subject,self.csp(self.current_subject)),1,'fontsize',15);            
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










