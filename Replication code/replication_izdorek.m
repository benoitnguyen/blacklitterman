% Benoit Nguyen 10/02/2016
% Paris 1 and Banque de France

% Replication of Izdorek paper "A step-by-step implementation of
% Black-Literman model"

clear all;

%% Data

% Historical mean (table 1)
histmu = [0.0315; 0.0175 ; -0.0639 ; -0.0286 ; -0.0675 ; -0.0054 ; -0.0675 ; -0.0526];

% risk aversion coefficient (table 1)
lambda = 3.07

% market cap weights (table 2)
wmkt=[19.34; 26.13; 12.09; 12.09; 1.34; 1.34; 24.18; 3.49];

% covariance matrix (table 5)
V=[0.001005 0.001328 -0.000579 -0.000675 0.000121 0.000128 -0.000445 -0.000437
0.001328 0.007277 -0.001307 -0.000610 -0.002237 -0.000989 0.001442 -0.001535
-0.000579 -0.001307 0.059852 0.027588 0.063497 0.023036 0.032967 0.048039
-0.000675 -0.000610 0.027588 0.029609 0.026572 0.021465 0.020697 0.029854
0.000121 -0.002237 0.063497 0.026572 0.102488 0.042744 0.039943 0.065994
0.000128 -0.000989 0.023036 0.021465 0.042744 0.032056 0.019881 0.032235
-0.000445 0.001442 0.032967 0.020697 0.039943 0.019881 0.028355 0.035064
-0.000437 -0.001535 0.048039 0.029854 0.065994 0.032235 0.035064 0.079958];

%Pi is the Implied Excess Equilibrium Return Vector (N x 1 column vector);
Pi=lambda*V*wmkt; 

% tau 
tau = 0.025; 

% names
name = {'US bonds','Int bonds','US large gowth','US large value','US small growth','US small value','Int dev equity','Int EME equity'};

%% Input views, P, Q and Omega

% P is the (K x N) matrix that identifies the assets involved in the views ;
%data taken fron eqn 7 in page 15
P=  [0  0   0     0    0     0     1 0
    -1  1   0     0    0     0     0 0
     0  0   0.9   -0.9 0.1   -0.1  0 0];

% Q is the (K x 1) column vector of views
Q=[5.25; 0.25; 2.00];

% omega is the covariance matrix of error terms of the views, scaled by tau
omega=tau*diag(diag(P*V*P'));


%% Posterior vector of return

Er=inv(inv(tau*V)+P'*inv(omega)*P)*(inv(tau*V)*Pi+P'*inv(omega)*Q);

%% Portfolio construction and comparison

% wmv MV optim with historical returns
wmv = inv(lambda*V)*histmu;

% The 100% confidence solution
Er100=Pi+tau*V*P'*inv(P*tau*V*P')*(Q-P*Pi);
w100=inv(lambda*V)*Er100;

% The weight vector from BL model
wbl=inv(lambda*V)*Er;

% market capitalization weights
wmkt=inv(lambda*V)*Pi;

% create weights bars (this one alone: very high weighting on the first
% asset)
bar(wmv);
set(gca,'XTickLabel',name);
rotateXLabels( gca(), 45 );
set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPosition', [0 0 7 6]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 6]); %Set the paper to have width 5 and height 5.
saveas(gcf, './slides/wmv', 'eps') %Save figure

% matrix of weights
W = horzcat(wmkt, w100, wbl);
g=bar(W)
set(gca,'XTickLabel',name);
rotateXLabels( gca(), 45 );
l{1}='Market weights'; l{2}='100% confidence'; l{3}='Black-Litterman';   
legend(g,l);
set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPosition', [0 0 7 6]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 6]); %Set the paper to have width 5 and height 5.
saveas(gcf, './slides/allw', 'eps') %Save figure


%% Portfolio perf 

disp(' ');
disp('*****portfolio statistics*****');

disp('Mkt portfolio');
wmv=wmv/100;
ExcessReturn_mv=wmv'*Pi
Variance_mv =wmv'*V*wmv

disp('Mkt portfolio');
wmkt=wmkt/100;
ExcessReturn_mkt=wmkt'*Pi
Variance_mkt =wmkt'*V*wmkt

disp('BL portfolio');
wbl=wbl/100;
ExcessReturn_BL=wbl'*Er
Variance_BL=wbl'*V*wbl

disp('100% confidence');
w100=w100/100;
ExcessReturn_100=wbl'*Er100
Variance_100=w100'*V*w100

% create table

temp=horzcat(ExcessReturn_mv,ExcessReturn_mkt, ExcessReturn_BL, ExcessReturn_100)
temp2=horzcat(Variance_mv, Variance_mkt, Variance_BL, Variance_100)
F=vertcat(temp,temp2)