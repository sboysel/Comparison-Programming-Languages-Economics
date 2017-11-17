%% Stochastic NGM
%
% Attribution: Jesus Fernandez-Villaverde
% Haverford, July 3, 2013
% Source: https://github.com/jesusfv/Comparison-Programming-Languages-Economics/blob/master/RBC_Matlab.m
%
% Adapted for ECON 602, Assignment 7 by
% Sam Boysel

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration

alpha = 0.36;    % Elasticity of output w.r.t. capital
beta  = 0.96;    % Discount factor

% Productivity values
vProductivity = [exp(-sqrt(0.007^2/(1 - 0.95^2)));
                 exp(sqrt(0.007^2/(1 - 0.95^2)))]';

% Transition matrix
mTransition = [0.975, 0.025;
               0.025, 0.975];

%% 2. Steady State

capitalSteadyState = (alpha*beta)^(1/(1-alpha));
outputSteadyState = capitalSteadyState^alpha;
consumptionSteadyState = outputSteadyState-capitalSteadyState;

fprintf(' y_ss = %2.6f, k_ss = %2.6f, c_ss = %2.6f\n', outputSteadyState, capitalSteadyState, consumptionSteadyState);
fprintf('\n')

% We generate the grid of capital
vGridCapital = 0.95*capitalSteadyState:0.00001:1.05*capitalSteadyState;

nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);

%% 3. Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity);
mValueFunction    = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 4. We pre-build output for each point in the grid

mOutput = (vGridCapital'.^alpha)*vProductivity;

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
                valueProvisional = log(consumption)+beta*expectedValueFunction(nCapitalNextPeriod,nProductivity);

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

fprintf(' My check = %2.6f\n', mPolicyFunction(1000,2));
fprintf('\n')

toc

%% 6. Plotting results

figure(1)

subplot(3,1,1)
plot(vGridCapital,mValueFunction)
xlabel('k')
ylabel('V(k, z)')
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Value Function')
legend('k_l','k_h','Location','southeast')

subplot(3,1,2)
plot(vGridCapital,mPolicyFunction)
xlabel('k')
ylabel("k'")
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy Function')
legend('k_l','k_h','Location','southeast')

% vExactPolicyFunction = alpha*beta.*(vGridCapital.^alpha);
% 
% subplot(3,1,3)
% plot((100.*(vExactPolicyFunction'-mPolicyFunction(:,2))./mPolicyFunction(:,2)))
% xlabel('Iteration')
% title('Comparison of Exact and Approximated Policy Function')

%set(gcf,'PaperOrientation','landscape','PaperPosition',[-0.9 -0.5 12.75 9])
%print('-dpdf','Figure1.pdf')
