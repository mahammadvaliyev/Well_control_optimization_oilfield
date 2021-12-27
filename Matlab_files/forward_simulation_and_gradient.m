function [cost,grad]=forward_simulation_and_gradient(q,rock,varargin)

%% 1. Wells setup
% get the liner indexs of positions and save into Prodpos and Injpos
load('Constant');
IndexProdpos =  [ sub2ind([nx,ny],1,1) ,  sub2ind([nx,ny],1,220)  ,sub2ind([nx,ny],60,220) , ...
    sub2ind([nx,ny],60,1)]';
IndexInjpos = [ sub2ind([nx,ny],5,148), sub2ind([nx,ny],50,89)]';
% save number of wells
NumProd = size(IndexProdpos,1);
NumInj  = size(IndexInjpos,1);
% get position information on Cartesian coordinates
[X_Inj, Y_Inj]=ind2sub([nx,ny],IndexInjpos);
[X_Prod, Y_Prod]=ind2sub([nx,ny],IndexProdpos);

%% 2. reshape Q
NT=10;
Q=reshape(q,NT,NumInj+NumProd)'; % q - 60x1; Q: 6x10 as '
load('Grid')
%% 3. Define well type and rates
load('Controls')
Controls.InjCtrlVals  = Q(1:NumInj,:)*(Controls.Max_Rate-Controls.Min_Rate)+Controls.Min_Rate;     % Reverse of minmax scaler - 2x10
Controls.ProdCtrlVals = Q(NumInj+1:NumInj+NumProd,:)*(Controls.Max_ProdBhp-Controls.Min_ProdBhp)+Controls.Min_ProdBhp;  % Set rate  4x10
Controls.InjCtrlTypes  = repmat({'rate'},NumInj,1); % create a column with character 'rate'; repmat repeat matrix
Controls.ProdCtrlTypes = repmat({'bhp'} ,NumProd,1); % create a column with character 'bhp'

%% 4. Difine fluid model
fluid  = initCoreyFluid('mu' , [1, 5] .* centi*poise, ...
    'rho', [1014, 859].*kilogram/meter^3, ...
    'n'  , [2, 2], 'sr', [0, 0], 'kwm', [1, 1]);
%% 5. Initial Wells

radius = 0.1;
W= [];

% Initial injector
for k=1:NumInj
    nameinj = ['inj',num2str(k)];
    W = verticalWell(W, Grid, rock,  X_Inj(k),Y_Inj(k), [],     ...
        'Type', Controls.InjCtrlTypes{k} , 'Val', abs(Controls.InjCtrlVals(k)), ...
        'Radius', radius, 'Name', nameinj, 'Comp_i', [1, 0], 'Sign', 1, 'InnerProduct', 'ip_tpf'); %  'Comp_i', [1, 0], , 'Sign', sign(InjRates(k))
end

% Initial producer

for k = 1:NumProd
    nameprod = ['prod', num2str(k)];
    W = verticalWell(W, Grid, rock,  X_Prod(k),Y_Prod(k),[],     ...
        'Type',Controls.ProdCtrlTypes{k}, 'Val',Controls.ProdCtrlVals(k), ...
        'Radius', radius, 'Name', nameprod, 'Comp_i', [1, 0], 'Sign', -1, 'InnerProduct', 'ip_tpf');
end
% System componeNTs -------------------------------------------------------
verbose = false; %verbose is for dipslaying info
S = computeMimeticIP(Grid, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, ...
    'InnerProduct', 'ip_tpf');
%mimetic finite volume solver for incompressible pressure
W = assembleWellSystem(Grid, W, 'Type', 'comp_hybrid');
% generate pressure linear system components for wells
%% Initial calculation
load sa
%horzcat horizontal concatenation
state = initResSol(Grid, 5000*psia, horzcat(1-sa,sa)); % set initial saturation of two phase
state.wellSol = initWellSol(W, 5000*psia);
%% Check sum of input and output
% 
%% box constraints for each well [min rate, max rate]
box = [repmat([0 inf], NumInj, 1); repmat([Controls.Min_ProdBhp Controls.Max_ProdBhp], NumProd, 1)];

%% Initial Schedule and Control
schedule = initSchedule(W, 'NumSteps', Controls.NumTsteps, 'TotalTime', ...
    Controls.totTime, 'Verbose', verbose);
controls = initControls(schedule, 'ControllableWells', (1:numel(W)), ...
    'MinMax', box, ...
    'Verbose', verbose, ...
    'NumControlsteps', Controls.NumTsteps);
uOpt = [Controls.InjCtrlVals;Controls.ProdCtrlVals]; %
[controls,schedule] = update_Controls_Scheduals(schedule,controls, uOpt);
verboseLevel=0;
simRes = runSchedule(state, Grid, S, W, rock, fluid, schedule, ...
    'VerboseLevel', verboseLevel);

%% calculate pressure and oil saturation
curr_id = 1:ngrids;
Press=[];
Oil_sat_final=[];
for i=1:NT+1
    Press=[Press simRes(i).resSol.pressure];
    Oil_sat_final=[Oil_sat_final simRes(i).resSol.s(curr_id,2)];
end
%

%%Save saturations
%Saturations_xopt=Oil_sat_final;
%save('Saturations_xopt','Saturations_xopt');
%% Get final cost results
objectiveFunction = str2func('simpleNPV'); %create handle
obj       = objectiveFunction(Grid, S, W, rock, fluid, simRes);
npv  = obj.val;
cost = -npv;

%% add gradient information
adjRes  = runAdjoint(simRes, Grid, S, W, rock, fluid, schedule, controls, objectiveFunction, 'VerboseLevel', verboseLevel);
    grad    = computeGradient(W, adjRes, schedule, controls);
    adjGrad = cell2mat(grad);
    GMat    = adjGrad;
    Gradient = -GMat;
    active_Grad   = Gradient*(Controls.totVol/Controls.totTime);
    gradJ    = active_Grad./norm(active_Grad);
    temp_grad           = active_Grad';
    gradient_vectorized = reshape(temp_grad, NT*size(active_Grad,1), 1);
    grad = gradient_vectorized;

%  %% save data
savedetail = false; %% defalut is setting to not save all details data
 for k = 1:length(varargin)
     if strcmpi(varargin{k},'SaveDetails') %string compare
         savedetail = varargin{k+1};
     end
 end
 
 if savedetail == true
        t=0;
        wres = cell([1, 4]);

        Prod = struct('t'  , []                  , ...
            'vpt', zeros([0, numel(W)]), ...
            'opr', zeros([0, numel(W)]), ...
            'wpr', zeros([0, numel(W)]), ...
            'wc' , zeros([0, numel(W)]));

        append_wres = @(x, t, vpt, opr, wpr, wc) ...
            struct('t'  , [x.t  ; t                  ], ...
            'vpt', [x.vpt; reshape(vpt, 1, [])], ...
            'opr', [x.opr; reshape(opr, 1, [])], ...
            'wpr', [x.wpr; reshape(wpr, 1, [])], ...
            'wc' , [x.wc ; reshape(wc , 1, [])]);


        [wres{:}] = prodCurves(W, simRes(1).resSol, fluid);
        Prod      = append_wres(Prod, t, wres{:});

        for i=2:NT+1
            t=simRes(i).timeInterval(2);
            [wres{:}] = prodCurves(W, simRes(i).resSol, fluid);
            Prod      = append_wres(Prod, t, wres{:});
        end

        wc_nan=find(isnan(Prod.wc));
        Prod.wc(wc_nan)=0;
        Wat_Prod=convertTo(Prod.wpr(:,NumInj+1:end),  stb/day);
        Wat_Inj=convertTo(Prod.wpr(:,1:NumInj),  stb/day);
        Oil_Prod=convertTo(Prod.opr(:,NumInj+1:end),  stb/day);
        
        Results.WWPR = Wat_Prod;
        Results.WOPR = Oil_Prod;
        Results.WWCT = Prod.wc(:,NumInj+1:end);
        Results.Saturation = Oil_sat_final;
        save('Results','Results')

        fid  = fopen('.\Result\Wat_Prod.json','wt')
        fprintf(fid,jsonencode(Wat_Prod));
        fclose(fid);
        fid  = fopen('.\Result\Wat_Inj.json','wt')
        fprintf(fid,jsonencode(Wat_Inj));
        fclose(fid);
        fid  = fopen('.\Result\Oil_Prod.json','wt')
        fprintf(fid,jsonencode(Oil_Prod));
        fclose(fid);
        fid  = fopen('.\Result\WCT.json','wt')
        fprintf(fid,jsonencode(Prod.wc));
        fclose(fid);
        fid  = fopen('.\Result\Pressure.json','wt')
        fprintf(fid,jsonencode(Press));
        fclose(fid);
        fid  = fopen('.\Result\Saturation.json','wt')
        fprintf(fid,jsonencode(Oil_sat_final));
        fclose(fid);
        fid  = fopen('.\Result\NPV.json','wt')
        NPV = struct('NPV',npv)
        fprintf(fid,jsonencode(NPV),'ConvertInfAndNaN',true);
        fclose(fid);
 end









