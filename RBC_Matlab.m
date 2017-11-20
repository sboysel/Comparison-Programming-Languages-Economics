%% Stochastic NGM
%
% Attribution: 
% - Author: Jesus Fernandez-Villaverde
% - Source: https://github.com/jesusfv/Comparison-Programming-Languages-Economics/blob/master/RBC_Matlab.m
%
% Adapted for ECON 602 - Assignment 7 by Sam Boysel

clear all

rng(123);

%%  Q2: Calibration

alpha = 0.36;    % Elasticity of output w.r.t. capital
beta  = 0.96;    % Discount factor

% Transition matrix
mTransition = [0.975, 0.025;
               0.025, 0.975];

% Productivity shocks
vProductivity = [exp(-sqrt(0.007^2/(1 - 0.95^2)));
                 exp(sqrt(0.007^2/(1 - 0.95^2)))]'; 

% Steady state
k_bar = (alpha*beta)^(1/(1-alpha));
y_bar = k_bar^alpha;
c_bar = y_bar-k_bar;

% Capital space
vGridCapital = [0.95*k_bar, 1.05*k_bar];

% Output space
mOutput = (vGridCapital'.^alpha)*vProductivity;                     

%% Q2: VFI

nGridCapital          = length(vGridCapital);
nGridProductivity     = length(vProductivity);
mOutput               = zeros(nGridCapital,nGridProductivity);
mValueFunction        = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew     = zeros(nGridCapital,nGridProductivity);
mPolicyFunction       = zeros(nGridCapital,nGridProductivity); % rows = [kl, kh], cols = [zl, zh]
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

while (maxDifference>tolerance)

    expectedValueFunction = mValueFunction*mTransition';

    for nProductivity = 1:nGridProductivity

        gridCapitalNextPeriod = 1;

        for nCapital = 1:nGridCapital

            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);

            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital

                consumption = mOutput(nCapital,nProductivity)-vGridCapital(nCapitalNextPeriod);
                valueProvisional = log(consumption)+beta*expectedValueFunction(nCapitalNextPeriod,nProductivity);

                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break;
                end

            end

            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;

        end

    end

    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;

    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
    end

end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
fprintf('\n')

%% Q4: expected value for output (LR mean)
e_y = 0.4875*0.9778*0.1806^0.36 + 0.0125*1.02267*0.1806^0.36 + 0.0125*0.9778*0.1996^0.36 + 0.4875*1.02267*0.1996^0.36;

%% Q5: Simulation
nPeriods = 10000;
y_sim = [];
k_sim = [];
c_sim = [];
y_sim(1) = 0;
k_sim(1) = vGridCapital(1); % Start at k_low

% Simulate shocks according to Markov chain
mc = dtmc(mTransition);
sim = simulate(mc, nPeriods);

% Simulate Stochastic NGM (T = 10,000)
for i = 1:nPeriods
    % Realize Shock
    shock = vProductivity(sim(i));
    % Policy Function
    if (shock>1)
        k_sim(i+1) = vGridCapital(2); % high state
    else
        k_sim(i+1) = vGridCapital(1); % low state
    end
    % Production
    y_sim(i) = shock*k_sim(i)^alpha;
    % Consumption
    c_sim(i) = y_sim(i) - k_sim(i+1);
end

mean(k_sim) % 0.1901
mean(y_sim) % 0.5502
mean(c_sim) % 0.3601
