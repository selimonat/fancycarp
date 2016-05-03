%% Magnet Paper: Implementation One Category

%% define a gridx
m = 0;
s = 1;
gridx_res = 1000;
gridx      = linspace(-50,50,gridx_res);

delta     = unique(diff(gridx));
delta = delta(1);
sum(exp(-((gridx-m)./s).^2))*delta


%% likelihood function: 
% p(S|T), this tells us the possible sensory signals caused by a real signal T.
mu_c     = [-8 0  20];
sigma_c  = [2 3   7 ];
T        = 0;
sigma_s  = 5;
p_s_given_t = @(x,m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
p_t_given_c = @(x,m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
p_c         = 1./length(sigma_c);
%% compute joint probabilities
for c = 1:3
    t=0;
    for T = gridx
        t = t + 1;
        p_s_t_c(:,t,c) = p_s_given_t(gridx,T,sigma_s).*p_t_given_c(T,mu_c(c),sigma_c(c)).*p_c;
    end
end
%% compute p(t,s|c): p(t,s,c)./p(c)
for c = 1:3
subplot(1,3,c)
%imagesc(p_s_t_c(:,:,c)./p_c)
tt = 0;
for T = gridx
    tt = tt+1;
    bla(:,tt) = p_s_given_t(gridx,T,sigma_s)*p_t_given_c(T,mu_c(c),sigma_c(c));
end
imagesc(bla);
end
%% compute p(c|s). 
for c = 1:3
    tt = 0;
    for T = gridx
        tt        = tt+1;
        bla(:,tt) = p_s_given_t(gridx,T,sigma_s)*p_t_given_c(T,mu_c(c),sigma_c(c));
    end
    p_c_s(:,c) = sum(bla,2);
end
%%

gridx             = linspace(-20,20,gridx_res); 
p                 = linspace(-20,20,gridx_res)';
face_id           = linspace(-8,8,8);
p_perc_given_face = @(x,m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
p_face            = @(x,m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
y                 = [];
for fi = face_id;
    y  = [ y p_perc_given_face(p,fi,4)];
end
joint = y.*repmat(p_face(1:8,4,5),1000,1);
imagesc(face_id,p,y);
plot(sum(joint,2));
%%
f   = @(m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
f   = [3 4 3 3 4 3 ];
f   = f./sum(f);
x = [10 0 0];x=x./sum(x);
pgf            = zeros(3,3);
pgf([1 5 9])   = x(1);
pgf([2 4 6 8]) = x(2);
pgf([4 6])     = sum(x(2:3))/2
pgf([3 7])     = x(3);
pf             = pgf*diag(f);
sgf            = [ .99 .99  .9 ; 0.01 0.01 .9];
sf             = sgf.*repmat(f,[2,1]);
fgs            = sf./repmat(sum(sf,2),1,3);
s              = sum(sf,2)
d              = (pgf*diag(fgs(1,:))*s(1));
c              = (pgf*diag(fgs(2,:))*s(2));
plot(sum(c,2),'color',rand(1,3))
hold on
%%
f   = @(m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
f   = ones(1,
f   = f./sum(f);
x = [10 0 0];x=x./sum(x);
pgf            = zeros(3,3);
pgf([1 5 9])   = x(1);
pgf([2 4 6 8]) = x(2);
pgf([4 6])     = sum(x(2:3))/2
pgf([3 7])     = x(3);
pf             = pgf*diag(f);
sgf            = [ .99 .99  .9 ; 0.01 0.01 .9];
sf             = sgf.*repmat(f,[2,1]);
fgs            = sf./repmat(sum(sf,2),1,3);
s              = sum(sf,2)
d              = (pgf*diag(fgs(1,:))*s(1));
c              = (pgf*diag(fgs(2,:))*s(2));
plot(sum(c,2),'color',rand(1,3))
hold on

%%
gridx_res = 1000;
gridx     = linspace(-60,60,gridx_res);
%%
precision_of_perception = 20;
prior_on_faces          = 2;
active_generalization   = 3;


x         = 1:8;
f         = @(m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
prior     = f(7,prior_on_faces) + f(3,prior_on_faces);
prior     = prior./sum(prior);
subplot(3,4,1);bar(prior);title('p(f)');
% plot(prior);
%
x         = gridx;  
f         = @(m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
m         = linspace(-30,30,8);
M         = [];
for mm = m
    x     = gridx;
    M     = [M f(mm,precision_of_perception)'];
end
perception_given_face = M;
perception_face = perception_given_face*diag(prior);
subplot(3,4,2);
imagesc(perception_given_face);axis square;title('p(p|f)');
subplot(3,4,3);
imagesc(perception_face);axis square;title('p(p,f)');
% joint     = M*diag(prior);
% imagesc(joint);
%
x                = 1:8
f                = @(m,s) exp(-((x-m)./s).^2)./(sum(exp(-((gridx)./s).^2)));
shock_given_face = f(8,active_generalization);
shock_given_face = shock_given_face./sum(shock_given_face);
%
shock_given_face = [1-shock_given_face; shock_given_face ];
subplot(3,4,4);
imagesc(shock_given_face);colorbar;colormap gray;title('p(cs+|face)');axis off

shock_face       = shock_given_face.*repmat(prior,[2 1]);
subplot(3,4,5);imagesc(shock_face);colorbar;title('p(cs+,face)');axis off
%
shock            = sum(shock_face,2);
subplot(3,4,6)
bar([0 1],shock);title('p(CS+)')
face_given_shock = shock_face./repmat(sum(shock_face,2),1,8);
subplot(3,4,7)
imagesc(face_given_shock);title('p(face|CS+)');colorbar;axis off;
%
perception_face_given_shock  = perception_given_face*diag(face_given_shock(1,:)).*shock(1);
subplot(3,4,8)
imagesc(perception_face_given_shock);colorbar;axis off;title('p(p,f | cs+ = 0)')
%
perception_face_given_shock  = perception_given_face*diag(face_given_shock(2,:)).*shock(2);
subplot(3,4,9)
imagesc(perception_face_given_shock);colorbar;axis off;title('p(p,f | cs+ = 1)');
%
subplot(3,4,10)
bar(sum(perception_face_given_shock(1:750,:),2));
hold on
[mm ii] = max(perception_given_face(:,end));
plot([ii ii],ylim,'r');
hold off













