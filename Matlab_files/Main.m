%%
% small optimization problem of well  control based on MRST tool
% Created by Junjie Yu and modified by Mahammad Valiyev.


%% 1. Standard Setup
clear;   % remove items from workspace
clear all; % clears variables, catched memory
close all; % close all figures
fclose('all'); %close all open files

tic % start of timing
MainDir = cd; %cd returns current directory

% !!  change to your own MRST file location
cd('C:\Users\Mahammad\OneDrive\Desktop\USC_2021\Research\MRST\mrst-2020b');
NewPath = cd; %assign newpath to directory of MRST
cd(MainDir);
addpath(genpath(NewPath)); % Use genpath in conjunction with addpath to add a folder and its subfolders to the search path.
rng('default'); %put random number generator back to its default settings

%% 2. Reservoir Dimension Set
nx = 60; ny = 220; nz = 1;% structure of reservoir - 13200 gridsrock
dx =nx*20*ft; dy =ny*10*ft; dz = nz*50*ft; % real size of reservoir
ngrids = nx*ny*nz; %grid number

%% 3. Set up grid, load reservoir properties and assign to rock
Grid = cartGrid([nx ny nz],[dx dy dz]); % initialize grid structure
Grid  = computeGeometry(Grid); % this function computes derived quantities such as cell volumes, centroids, areas
save('Constant.mat','nx','ny','nz','dx','dy','dz','ngrids'); % save all variables from the workspace in mat file 
save('Grid','Grid') %first '' is for name of files, then variables

load('rockdata') %load file from file into workspace
% rockdata is struct, like dictionary in python shape:100*[13200*3 13200*1]
% 1 mD= 10^-15 D
rock = rockdata(82); %load 82th realization; rock is struct with 2 fields
% shape [13200*3 13200*1]
% we specified heterogeneous rock properties
%% 4. Control settings 
Controls.NumTsteps             = 10;                 % reporting times
Controls.NumofWellControlsteps = 10;                 % control time indices
Controls.totTime               = 10*365*day;         % reservoir life time in sec
Controls.totVol                = 1*sum(poreVolume(Grid, rock));  % in m^3; poreVolume is for 1 grid (dx*dy*dz*porosity)
NT                             = Controls.NumofWellControlsteps;
Controls.Max_ProdBhp  = 4500*psia; %MRST calculates in Pa; 1psia=6900 Pa
Controls.Min_ProdBhp = 500*psia; % in Pa
Controls.Min_Rate = 100*stb/day; % 1 stb/day=0.159 m3/day; in m3/s
Controls.Max_Rate = 700*stb/day; % in m3/s
% save('ControlPara.mat','-struct','Controls')
save('Controls','Controls')

%% 5. Wells setup
% get the linear indeces of positions and save into Prodpos and Injpos

% subind - convert subscripts to linear indices
% counts columnwise, 1st finishes 1st column then passes to 2nd 1,1 - 2,1
% 3,1 and indices are 1,2,3
% for function syntax sub2ind(this is size of cube; then 1st index is row,
% 2nd is column; so function takes input as (size of matrix, row, col)
% returns number

IndexProdpos =  [ sub2ind([nx,ny],1,1) ,  sub2ind([nx,ny],1,220)  ,sub2ind([nx,ny],60,220) , ... % 60x220 size of model
    sub2ind([nx,ny],60,1)]'; % returns a column vector [1 13141 13200 60]
IndexInjpos = [ sub2ind([nx,ny],5,148), sub2ind([nx,ny],50,89)]'; %returns a columns vector [8825 5330]'

% save number of wells
NumProd = size(IndexProdpos,1); %4
NumInj  = size(IndexInjpos,1);  %2

% get position information on Cartesian coordinates
[X_Inj, Y_Inj]=ind2sub([nx,ny],IndexInjpos); %returns row and column number 2 outputs: 1st row, 2nd column 3- 3,1 if size if 3x3
[X_Prod, Y_Prod]=ind2sub([nx,ny],IndexProdpos);

% Initial injector
radius = 0.1;
W = [];

for k=1:NumInj %2
    nameinj = ['inj',num2str(k)];
    W = verticalWell(W, Grid, rock,  X_Inj(k),Y_Inj(k), [],     ...
        'Type', 'rate' , 'Val', 100, ...
        'Radius', radius, 'Name', nameinj, 'Comp_i', [1, 0], 'Sign', 1, 'InnerProduct', 'ip_tpf'); %  'Comp_i', [1, 0], , 'Sign', sign(InjRates(k))
    % x(k), j(k) is position, [] z index, val is in m3/sec
    %compi fluid comp only for injector wells
    % sign 1 for inj; sign -1 for prod
    % innerProduct: used for consistent discretizations discussed in Chapter 6.
end

% Initial producer
for k = 1:NumProd
    nameprod = ['prod', num2str(k)];
    W = verticalWell(W, Grid, rock,  X_Prod(k),Y_Prod(k),[],     ...
        'Type','bhp', 'Val',1000*psia, ...
        'Radius', radius, 'Name', nameprod, 'Comp_i', [1, 0], 'Sign', -1, 'InnerProduct', 'ip_tpf');
end

%% 6. Permeability plot
myview = struct('vw',   [60, 30],    ...  % view angle
    'zoom', 1,        ...  % zoom
    'asp',  [1 1 0.5],  ...  % data aspect ratio
    'wh',   30,         ...  % well height above reservoir
    'cb',   'horiz'     ...  % colorbar location
    );


figure(1) %creates figure window

%gcf returns the handle of current figure
% outer position [ x0 y0 width height ]
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% plot celldata is for plotting grid
plotCellData(Grid,log10(convertTo(rock.perm(:,1),milli*darcy)),'EdgeColor','k','EdgeAlpha',0.1);
%plotCellData(Grid,(convertTo(rock.perm(:,1),milli*darcy)),'EdgeColor','k','EdgeAlpha',0.1);
plotWell(Grid,W(1:2),'height',50,'color','b','fontsize',30);
hold on

plotWell(Grid,W(3:6),'height',50,'color','r','fontsize',30);
set(gca,'dataasp',myview.asp); %gca returns the handle to the current axis
% relative length of data units along each axis [dx dy dz] [1 2 1] means y
% is smaller 2 times
axis off, view(myview.vw); zoom(myview.zoom), colormap(jet), axis tight %set tge axis limits to the range of data
% cs = [1 500 1000 1500 2000 2500 3000 3500 4000];
cb = colorbar
cb.Location = 'southoutside';
cb.FontSize = 30;
cb.Label.FontSize = 30;
%axis tight off, view(65,30)
hold off

%% 7. Porosity plot

figure(2) %creates figure window

%gcf returns the handle of current figure
% outer position [ x0 y0 width height ]
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% plot celldata is for plotting grid
plotCellData(Grid,rock.poro,'EdgeColor','k','EdgeAlpha',0.1);
plotWell(Grid,W(1:2),'height',50,'color','b','fontsize',30);
hold on

plotWell(Grid,W(3:6),'height',50,'color','r','fontsize',30);
set(gca,'dataasp',myview.asp); %gca returns the handle to the current axis
% relative length of data units along each axis [dx dy dz] [1 2 1] means y
% is smaller 2 times
axis off, view(myview.vw); zoom(myview.zoom), colormap(jet), axis tight %set tge axis limits to the range of data
% cs = [1 500 1000 1500 2000 2500 3000 3500 4000];
cb = colorbar
cb.Location = 'southoutside';
cb.FontSize = 30;
cb.Label.FontSize = 30;
%axis tight off, view(65,30)

%% 7. Saturation plot
load('so')
sa=ones(13200,1);

figure(3) %creates figure window

%gcf returns the handle of current figure
% outer position [ x0 y0 width height ]
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% plot celldata is for plotting grid
plotCellData(Grid,sa,'EdgeColor','k','EdgeAlpha',0.1);
plotWell(Grid,W(1:2),'height',50,'color','b','fontsize',30);
hold on

plotWell(Grid,W(3:6),'height',50,'color','r','fontsize',30);
set(gca,'dataasp',myview.asp); %gca returns the handle to the current axis
% relative length of data units along each axis [dx dy dz] [1 2 1] means y
% is smaller 2 times
axis off, view(myview.vw); zoom(myview.zoom), colormap(jet), axis tight %set tge axis limits to the range of data
% cs = [1 500 1000 1500 2000 2500 3000 3500 4000];
cb = colorbar
cb.Location = 'southoutside';
cb.FontSize = 30;
cb.Label.FontSize = 30;
cb.Limits=[0 1];
%cb.Color=[1 0 0];
%axis tight off, view(65,30)
%% 8. Injector and Producer rates
% based on partial of porevolumn
PV_injected=1;
Init_Rates = zeros(NumProd+NumInj,NT); % 6x10
InjRates   =  PV_injected*(1/NumInj ) * ones(NumInj,NT ); %2x10 initial injection rates

% provide some initial control and normalize it based on minmax scaler 
ProdBhp = (ones(NumProd, NT)*2500*psia-Controls.Min_ProdBhp)./(Controls.Max_ProdBhp-Controls.Min_ProdBhp); % 4x10
InjRates = (ones(NumInj, NT)*400*stb/day-Controls.Min_Rate)./(Controls.Max_Rate-Controls.Min_Rate);  % 2x10

Q=[InjRates;ProdBhp]; % Q represents number of wells * number of years matrix 6x10
q=reshape(Q',NT*(NumProd+NumInj),1); % q represents number of wells * number of years linear column vector reshape(A,[x,y]) 60x1;
x0 = q;

%% 9. Well control optimization

% set low and up bound
lb=[zeros(1,NT*NumInj),zeros(1,NT*NumProd)]; % 1x60
ub=[ones(1,NT*NumInj),ones(1,NT*NumProd)];  % 1x60

options = optimoptions('fmincon','Algorithm','interior-point','UseParallel',true,'HessianApproximation','bfgs','GradObj','on',...
'MaxIter',100,'Display','iter','FunValCheck','on','MaxFunEvals',5000,'Diagnostics','on','TolFun',1e-16, 'TolX', 1e-12,'TolCon',10^-16....
,'PlotFcn','optimplotfval','OutPutFcn',@outfun);

% create a folder that recored information during optimiztaioon 
history.x=[];
history.fval = [];
history.gradient = [];
save('history','history');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) forward_simulation_and_gradient(x,rock),x0,[],[],[],[],lb,ub,[],options);
%exitflag describes exit condition, output is structure
%lambda is struct with fields containing lagrange multipliers at solution
%

% save history data if you need
load('history')
save('history','history') % name as you need


%% Some tips:

% 1) The above example assumes analytical gradient is available which was
%    calculated by Adjoint appraoch, if you assume analytical gradient is
%    not avalble, use forward_simulation instead of
%    forward_simulation_and_gradient, it will calculate gradient by finite
%    difference dimension by dimension. Also the option 'GradObj' should be
%    set as 'off'.

% 2) When you want to check more details about simulation results instead
% of only NPV, you can set rerun forward_simulation() 
% but set forward_simulation(...,'SaveDetails',True),it will save a Results
% file that includes more detilas (See code for clear understand)




    















