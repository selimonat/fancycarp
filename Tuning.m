classdef Tuning < handle
    properties (Hidden)
        visualization =0;
        options = optimset('Display','none','maxfunevals',10000,'tolX',10^-12,'tolfun',10^-12,'MaxIter',10000,'Algorithm','interior-point');
        
    end
    %tuning object, can contain any kind of fear-tuning SCR, rating etc.
    properties
        x =[];
        y =[];
        y_mean = [];
        y_std  = []; 
        groupfit
        singlesubject
        
%         Est
%         Likelihood
%         ExitFlag
%         pval = NaN;
%         dof
%         
%         ss_residuals
%         residuals
    end
    
    methods
        function tuning = Tuning(data,varargin)
            %data is anything that has a x and y fields.
            tuning.x = data.x;
            tuning.y = data.y;            
            for x  = unique(tuning.x(:)')                
                i             = tuning.x == x;
                tuning.y_mean = [tuning.y_mean mean(tuning.y(i))];
                tuning.y_std  = [tuning.y_std  std(tuning.y(i))];
            end
        end
                
        function SingleSubjectFit(self,funtype)
            %fit FUNTYPE to each individual subject
            for ns = 1:size(self.x,1)
                fprintf('Fitting subject %03d\n',ns)
                self.singlesubject{ns} = self.Fit(self.x(ns,:),self.y(ns,:),funtype); 
            end
        end
        
        function GroupFit(self,funtype)            
            %pools different subjects and fit FUNTYPE
            self.groupfit = self.Fit(self.x(:),self.y(:),funtype);
            %
            figure;
            set(gcf,'position',[1038         479         373         327]);
            bincenters    = linspace(min(self.y(:)),max(self.y(:)),10);
            [y x]         = hist(self.y,bincenters);
            y             = y./repmat(sum(y),size(y,1),1)*100;%make it a percentage
            imagesc(unique(self.x),x,y,[0 100]);
            colormap hot;
            thincolorbar('vert');
            axis xy;
            hold on;
            errorbar(unique(self.x)',self.y_mean,self.y_std,'go');            
            plot(self.groupfit.x,self.groupfit.fit);
            plot(self.groupfit.x,self.groupfit.fit+self.groupfit.Est(end),'--');
            plot(self.groupfit.x,self.groupfit.fit-self.groupfit.Est(end),'--');
            hold off;
        end
        
        function result = Fit(self,x,y,funtype)
            %Fits FUNTYPE to a tuning defined in x and y
                        
            %% set the function to be fitted
            x        = x(:);
            y        = y(:);%make it sure to have columns
            y        = y + rand(length(y),1)*eps;
            CONSTANT = 0;%will be added to all the data points                       
            if funtype == 1
                result.fitfun = @(x,p) repmat(p(1),length(x),1);
                L           = [ min(y)  0];%mean sigma
                U           = [ max(y)  2*std(y)];
                result.dof    = 1;
                result.funname= 'null';
            elseif funtype == 2                
                result.fitfun = @(x,p) make_gaussian_fmri(x,p(1),p(2),p(3));%2 amp, std, offset
                L           = [-range(y)*2    0       mean(y)-range(y)*2          .01    ];
                U           = [range(y)*2     180      mean(y)+range(y)*2    std(y(:)+rand(length(y),1).*eps)*2 ];%
                result.dof    = 3;
                result.funname= 'gaussian';                
            elseif funtype == 3                
                result.fitfun = @(x,p) make_gaussian_fmri_zeromean(x,p(1),p(2));%2 amp fwhm
                L           = [ -range(y)*2  0     .01    ];
                U           = [  range(y)*2  180       std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 3;
                result.funname= 'gaussian_ZeroMean';
                %detect the mean, store it and subtract it
                CONSTANT    = mean(y);
                y           = y-CONSTANT;%we are not interested in the baseline, just remove it so we don't need to estimated it.                
            elseif funtype == 4                
                result.fitfun = @(x,p) make_gaussian_fmri_tau(x,p(1),p(2),p(3));%amp, tau, offset
                L           = [-range(y)*2    0    mean(y)-range(y)*2          .01    ];
                U           = [range(y)*2     180  mean(y)+range(y)*2    std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 3;
                result.funname= 'gaussian_tau';                
            elseif funtype == 5
                result.fitfun = @(x,p) VonMises_fmri(x, p(1), p(2), p(3));
                L           = [ -range(y)*2  10.^-16  mean(y)-range(y)*2   0];%amp kappa offset
                U           = [  range(y)*2  180       mean(y)+range(y)*2   std(y)];
                result.dof    = 3;
                result.funname= 'vonmises';                
            elseif funtype == 6
                result.fitfun = @(x,p) make_gabor1d_ZeroMean(x,p(1),p(2),p(3),p(4));%amp std freq baseline
                L           = [ -range(y)*2  0    1  -range(y)*2    .01    ];
                U           = [  range(y)*2  180   4   range(y)*2  std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 6;
                result.funname= 'gabor';
            elseif funtype == 7
                result.fitfun = @(x,p) p(1)*cos(x*p(2)) + p(3);%amp std freq baseline
                L           = [ -range(y)*2  1  -range(y)*2    .01    ];
                U           = [  range(y)*2  4   range(y)*2 std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 3;
                result.funname= 'cosine';
            elseif funtype == 8
                result.fitfun = @(x,p) VonMises(x,p(1),p(2),p(3),p(4));%amp,kappa,centerX,offset
                L             = [ min(y(:))-std(y)    0.1         min(x)   min(y(:))-std(y)   eps ];
                U             = [ max(y(:))+std(y)  max(abs(x))   max(x)   max(y(:))+std(y)   std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 3;
                result.funname= 'vonmisses_mobile';
            end
            %% set the objective function                        
            result.likelihoodfun  = @(params) sum(-log( normpdf( y - result.fitfun( x,params(1:end-1)) , 0,params(end)) ));
            
            %% Initial estimation of the parameters            
            %if gabor or gaussian, make a grid-estimatation
            if funtype > 1
                Init = RoughEstimator(x,y,result.funname,L,U)';
            else %null model
                Init = mean(y);
            end
           
            %% estimate initial values for sigma_noise
            %based on the likelihood of sigma given the data points assuming a Gaussian normal distribution.
            tsample      = 1000;
            sigmas       = linspace(0.001,std(y)*3,tsample);
            PsigmaGiveny = sigmas.^-length(y) .* exp( - (1./(2.*sigmas.^2) ) .* sum((y - result.fitfun(x,Init)).^2) );
            [m i]        = max(PsigmaGiveny);
            Init         = [Init sigmas(i)];%            
            %% Optimize!
            try
                [result.Est, result.Likelihood, result.ExitFlag]  = fmincon(result.likelihoodfun, Init, [],[],[],[],L,U,[],self.options);
                
                
                result.Likelihood = result.likelihoodfun(result.Est);
            catch
                result.Est         = Init;
                result.Likelihood = result.likelihoodfun(result.Est);
                result.ExitFlag   = 1;
            end
            %% get the null fit
            null.fitfun     = @(x,p) repmat(p(1),length(x),1);
            null.L          = [ min(y)  0];%mean sigma
            null.U          = [ max(y)  2*std(y)];
            null.Init       = Init;
            null.funny      = @(params) sum(-log( normpdf( y - null.fitfun( x, params(1:end-1)) , 0,params(end)) ));
            y               = y + rand(size(y))*eps;            
            [null.Est, result.null_Likelihood, null.ExitFlag]  = fmincon( null.funny, [mean(y) std(y)], [],[],[],[], null.L , null.U , [], self.options);
            %% Compute a residual
            %This is the RMS differences between the model and data points...            
            result.residuals      = y-result.fitfun(x,result.Est);
            result.ss_residuals   = sqrt(mean((y-result.fitfun(x,result.Est)).^2));                                   
            %% Get the likelihood            
            if ~isempty(result.null_Likelihood)
                result.dof   = result.dof - 1;%the DOF of the null hypothesis is 0.
                result.pval  = -log10(1-chi2cdf(-2*(result.Likelihood - result.null_Likelihood),result.dof) + eps);
            end
            result.x   = unique(x);
            result.fit = result.fitfun(result.x,result.Est);
        end
                
    end
end