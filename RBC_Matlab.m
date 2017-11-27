%% Basic RBC model with full depreciation
%
% Jesus Fernandez-Villaverde
% Haverford, July 3, 2013

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration

theta = 0.36;     % Elasticity of output w.r.t. capital
bbeta  = 0.99;    % Discount factor
delta = 0.025;    % depreciation
A = 1.72;
rho = 0.95;

% Productivity values
vProductivity = [0.9792; 0.9896; 1.0000; 1.0106; 1.0212]';

% Transition matrix
mTransition   = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
                 0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
                 0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
                 0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
                 0.0000, 0.0000, 0.0000, 0.0273, 0.9727];

%% 2. Steady State

hoursSteadyState = ((1 - theta) * ((1/beta) - 1 + delta) / ((A + 1 - theta) * ((1/beta) - 1 + delta) - theta * A * delta));
capitalSteadyState = (((1/beta) - 1 + delta) * (1 / theta))^(1 / (theta - 1)) * hoursSteadyState;
outputSteadyState = (((1/beta) - 1 + delta) * (1 / theta))^(theta / (theta - 1)) * hoursSteadyState;
consumptionSteadyState = outputSteadyState - delta * capitalSteadyState;

fprintf(' h_ss = %2.6f, k_ss = %2.6f, y_ss = %2.6f, c_ss = %2.6f\n', hoursSteadyState, capitalSteadyState, outputSteadyState, consumptionSteadyState);
fprintf('\n')

% We generate the grid of capital
vGridCapital = 0.5*capitalSteadyState:0.00001:1.5*capitalSteadyState;
vGridHours = 0.5*hoursSteadyState:0.00001:1.5*hoursSteadyState;

nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);

%% 3. Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity);
mValueFunction    = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 4. We pre-build output for each point in the grid

mOutput = (vGridCapital'.^theta)*vProductivity;

%% 5. Main iteration

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

while (maxDifference>tolerance)

    expectedValueFunction = mValueFunction*mTransition';

    for nProductivity = 1:nGridProductivity

        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;

        for nCapital = 1:nGridCapital

            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);

            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital

                consumption = mOutput(nCapital,nProductivity)-vGridCapital(nCapitalNextPeriod);
                valueProvisional = (1-bbeta)*log(consumption)+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);

                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
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

fprintf(' My check = %2.6f\n', mPolicyFunction(1000,3));
fprintf('\n')

toc

%% 6. Plotting results

figure(1)

subplot(3,1,1)
plot(vGridCapital,mValueFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Value Function')

subplot(3,1,2)
plot(vGridCapital,mPolicyFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy Function')

vExactPolicyFunction = theta*bbeta.*(vGridCapital.^theta);

subplot(3,1,3)
plot((100.*(vExactPolicyFunction'-mPolicyFunction(:,3))./mPolicyFunction(:,3)))
title('Comparison of Exact and Approximated Policy Function')

%set(gcf,'PaperOrientation','landscape','PaperPosition',[-0.9 -0.5 12.75 9])
%print('-dpdf','Figure1.pdf')
