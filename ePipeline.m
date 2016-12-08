%eColi Story
%master script of eColi machine learning pipeline
% T.T.O. 10/22/2016
clc
clear
%% load data
load ecoliModels
load RawData

%% choose model
% model=e_coli_core;
model=iJO1366;

%%  set up basic parameters
smallest=1e-12;%smallest number used
numflux= length(mRxns); % number of different fluxes measured
numdat=size(Fluxes,1); % number of data points
numRxns=length(model.rxns); % number of reactions

%% 13 C MFA and genome scale data integration using possibility MFA
% data_integration
load Data1a
%%
%Task #A Compute energy parameters
[EnergyData,EnergyLabels]= computeEnergyParameters(v,model,smallest);

%Task #B compute flux ratios
[FluxRatioData,RatioLabels]= computeFluxRatios(v);

%Task #C PCA analysis

