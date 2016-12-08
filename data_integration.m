%% 13 C MFA and genome scale data integration
% using possibility MFA

load RawData
%%  set up basic parameters
smallest=1e-12;%smallest number used
numflux= length(mRxns); % number of different fluxes measured
numdat=size(Fluxes,1); % number of data points
numRxns=length(model.rxns); % number of reactions
nErrorLevels=3;
ErrorLevels=[0.05 .15 .5]';
ErrorAllowed=zeros(numdat,1);

sp_Rxns={'EX_co2_e','EX_o2_e','THD2pp'};

tt_Rxns={'BIOMASS_Ec_iJO1366_core_53p95M','EX_glc__D_e'};
% tt_Rxn_ID=findRxnIDs(model,tt_Rxns);
% sp_Rxn_ID=findRxnIDs(model,sp_Rxns);
% Exclude Porin transport reactions and spontaneous reactions
Porin_ID=find(~cellfun(@isempty,strfind(model.subSystems,'Porin')));
Spont_ID=find(strcmp(model.grRules,'s0001'));
Excl_Rxn_ID=union(Porin_ID,Spont_ID);

% model.lb(Excl_Rxn_ID)=0;
% model.ub(Excl_Rxn_ID)=0;
% ttVec=[0.1 0.15 0.2 0.25 .5 1 1.5 3 4 5]';
% ttVec=[linspace(6,9,3) linspace(9, 10,10)];
possCut=0.9;
% 
v=zeros(numRxns,numdat);
% vmin=v;
% vmax=v;
poss=zeros(numdat,1);
% ttUsed=zeros(numdat,1);
Corrected=zeros(numdat,1);

%%
for j=1:numdat
%     if ttUsed(j)==0
        
        % formulate the problem
        modelP.S=model.S;
        modelP.rev=model.rev;
        modelP.lb=model.lb;
        modelP.ub=model.ub;

        %apply the measured flux data on the flux bounds using error limits
        tempData=Fluxes(j,:);
        % remove problematic reactions
        tempUsedFluxes=~isnan(tempData) & ~ismember(mRxns',sp_Rxns);

        %         tempUsedFluxes=~isnan(tempData);
        tempRxns=mRxns(tempUsedFluxes);
        tempIndex=findRxnIDs(model,tempRxns);

        tempFluxes=tempData(tempUsedFluxes);

        % Measured fluxes and their uncertainty in possibilistic terms:
        intFP = max(0.01,abs(tempFluxes.*0.05));
        intLP = max(0.02,abs(tempFluxes.*8));
%         for i=1:length(ttVec)


%             %tighter tolerances for 'tight' reactions
%             ttID=ismember(tempRxns,tt_Rxns);
%             intFP(ttID) = max(0.01,abs(tempFluxes(ttID).*.05));
%             intLP(ttID) = max(0.02,abs(tempFluxes(ttID).*ttVec(i)));


            %first try out the regular FBA optimization
            %         model.lb(tempIndex)=tempFluxes-intLP;
            %         model.ub(tempIndex)=tempFluxes+intLP;

            %
            %         testsol=optimizeCbModel(model,[],'one');
            %         testsol2=optimizeCbModel(model);
            %         rmID=testsol.x==0 & abs(testsol2.x)>=200;
            %
            %         modelP.lb(rmID)=0;
            %         modelP.ub(rmID)=0;

            %         model.lb(rmID)=0;
            %         model.ub(rmID)=0;

            % Initialize the possibilistic problem with the model constraints
            [PossProblem] = define_MOC(modelP);
            [PossMeasurements]=define_PossMeasurements(tempFluxes,intFP,intLP);

            % Extend the possibilistic problem adding the measurements-constraints
            [PossProblem] = define_MEC(PossProblem, PossMeasurements, tempIndex);

            % Perform a point_wise estimation (an unreliable estimate)
            [v(1:numRxns,j), poss(j),diagnostic(j)]=solve_maxPoss(PossProblem);
            [v(1:numRxns,j),Corrected(j)]=removeInfeasibleLoops(model,v(1:numRxns,j));

%             if poss(j)>=possCut
%                 ttUsed(j)= ttVec(i);
%                 break
%             end
%         end

        % Perform an interval estimation with maximum possibility
        %         [vmin(1:numRxns,j),vmax(1:numRxns,j)]=solve_maxPossIntervals(PossProblem);

        % Compute the flux distirbution and a estimation
        %         possibilities = [0.01:0.01:0.439];
        %         for f=1:numRxns
        %             [min_p(:,f),max_p(:,f)]=solve_PossInterval(PossProblem,possibilities,f,'none');
        %         end
        disp(j)
%     end
end


% C = strsplit(Genes1{27},',')

% glucose and growth rate tests
scatter(Fluxes(:,end),v(8,:))
figure
scatter(Fluxes(:,11),v(164,:))

save Data1b v poss Corrected  
% Data1b is with 8