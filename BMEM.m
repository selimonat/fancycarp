classdef BMEM < handle
    %creates a Bayesian Magnet Effect Model object
    properties (Constant)
        x                = linspace(-135,180,8);
        options          = optimset('Display','iter','maxfunevals',10000,'tolX',10^-16,'tolfun',10^-16,'MaxIter',10000,'Algorithm','interior-point');
        phase            = 3;%Which experimental phase
    end
    properties
        prior_function   = 'vonmises'
        param_kappa_csp        = 1;
        param_amp_csp          = 1;
        param_kappa_prior      = 1;
        param_kappa_perceptual = 2;
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
    end
    properties (Hidden)
        
        fit_estimates   = [];
        fit_likelihood  = [];
        fit_exitFlag    = [];
    end
    
    methods        
        function bmem = BMEM(G)
            %creates a Bayesian Magnet Effect Model object
            %G is simply a Group object
            %%
            bmem.subjects       = G.ids;
            bmem.csp            = G.csps;
            %collect the ratings.
            for nrun = 1:3
                bmem.ratings(nrun).y = NaN(length(G.subject),8);
                for ns = 1:length(G.subject)
                    if ~isempty(G.subject{ns}.ratings(nrun).y)
                        R                          = mean(G.subject{ns}.ratings(nrun).y,1);%extract it
                        bmem.ratings(nrun).y(ns,:) =  (R-min(R))./max(R);
                    end
                end
            end
            bmem.group = G;
        end
        function [ExpVariance] = fitTuningModels(self,nfuns)
            %Will run different models specified in NFUNS and return the
            %unexplained variance in EXPVARIANCE.
            
            for nfun = nfuns(:)'
                for run = 1:3
                    self.group.model_ratings(run,nfun);
                    fit                   = self.group.fit_results.y_fitted;
                    data                  = self.group.ratings(run).y_mean(:,1:8);
                    ExpVariance(nfun,run) = nanmean(nanmean((data -fit).^2,2));
                end
            end
        end              
        function out           = p_given_f(self,location)
            %p(p|f) the perceptual likelihood function
            out = [];
            c = 0;
            for l = location(:)'
                c        = c + 1;
                out(:,c) = Tuning.VonMises(self.x,1,self.param_kappa_perceptual,l,0);
            end
        end        
        function out           = csp1_given_face(self)
            %p(CS+ = 1 | face) the active generalization component.
            dummy        = Tuning.VonMises(self.x,self.param_amp_csp,self.param_kappa_csp,0,0);
            out          = repmat(dummy,size(self.subjects,1),1);
        end        
        function out = get.f(self)
            %p(f) returns prior on faces
            current_subject        = self.current_subject;
            location               = (3.5-self.csp(current_subject))*45;
            location2              = (7.5-self.csp(current_subject))*45;
            if strcmp(self.prior_function,'vonmises')
                dummy                 = Tuning.VonMises(self.x,1,self.param_kappa_prior,location,0);
                dummy                 = dummy + Tuning.VonMises(self.x,1,self.param_kappa_prior,location2,0);;
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
    end
    
    
    %% Basic Parameter optimization
    methods
        function find_params(self)                        
            %Matlab's gradient descent to find the free parameters of BMEM.
            %Parameters are optimized with respect to self.cost_function.
            
            dummy                 = @(params) self.cost_function(params);
            [self.fit_estimates, self.fit_likelihood, self.fit_exitFlag]  = ...            
                fmincon(dummy, [1 1 1 1], [],[],[],[],[0 0 0 0],[15 1 15 15],[],self.options);
                        
            %overwrite the current parameters.
            self.param_kappa_csp        = self.fit_estimates(1);
            self.param_amp_csp          = self.fit_estimates(2);
            self.param_kappa_perceptual = self.fit_estimates(3);
            self.param_kappa_prior       = self.fit_estimates(4);
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
            bar(self.p_f_given_csp1);
            title('Model');
            %
            subplot(2,2,4);
            bar(self.ratings(self.phase).y(self.current_subject,:));
            title('Ratings Data');
            %
            EqualizeSubPlotYlim(gcf);
            supertitle(sprintf('sub:%i, CS+:%i',self.current_subject,self.csp(self.current_subject)),1);
        end
        function out = cost_function(self,params)
            %out = cost_function(self,params);
            %
            % Measures the discrepancy between the BMEM and ratings.
            
            self.param_kappa_csp        = params(1);
            self.param_amp_csp          = params(2);
            self.param_kappa_perceptual = params(3);
            self.param_kappa_prior       = params(4);
            
            %
            data                  = self.ratings(self.phase).y(self.current_subject,:);
            current_fit           = self.p_f_given_csp1';
            out                   = mean(( current_fit - data ).^2);
        end
    end
end