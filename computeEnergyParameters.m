function [EnergyTable,EnergyLabels]= computeEnergyParameters(v,model,smallest)
% compute parameters that describe the energy metabolism of the cell
% create variable for energy parameters


EnergyLabels={'EnergyProduced','EnergyUsedForBiomassGrowth',...
    'EnergyUsedForMaintenance','EnergyProducedPerUnitCarbonSource',...
    'MaintenanceEnergyPerUnitBiomass','MaintenanceFracn',...
    'OxidativePhosEnergy','nadphProduced','nadhProduced','fadh2Produced'};


numE=length(EnergyLabels);
numdat=size(v,2);
EnergyData=zeros(numdat,numE);
EnergyTable=array2table(EnergyData,'VariableNames',EnergyLabels);

atpID=findMetIDs(model,'atp_c');
biomassID=findRxnIDs(model,'BIOMASS_Ec_iJO1366_core_53p95M');
oxdPhosID=findRxnIDs(model,'ATPS4rpp');
nadhID=findMetIDs(model,'nadh_c');
nadphID=findMetIDs(model,'nadph_c');
fadh2ID=findMetIDs(model,'fadh2_c');

atpVec=full(model.S(atpID,:));
nadhVec=full(model.S(nadhID,:));
nadphVec=full(model.S(nadphID,:));
fadh2Vec=full(model.S(fadh2ID,:));

% finding energy produced
atpMatrix=repmat(atpVec',1,numdat);
dotResult=atpMatrix.*v;
dotResult(dotResult<0 | abs(dotResult)<=smallest)=0; %only considering atp produced
EnergyTable.EnergyProduced= sum(dotResult,1)';

nadhMatrix=repmat(nadhVec',1,numdat);
dotResult=nadhMatrix.*v;
dotResult(dotResult<0 | abs(dotResult)<=smallest)=0; 
EnergyTable.nadhProduced= sum(dotResult,1)';

nadphMatrix=repmat(nadphVec',1,numdat);
dotResult=nadphMatrix.*v;
dotResult(dotResult<0 | abs(dotResult)<=smallest)=0; 
EnergyTable.nadphProduced= sum(dotResult,1)';

fadh2Matrix=repmat(fadh2Vec',1,numdat);
dotResult=fadh2Matrix.*v;
dotResult(dotResult<0 | abs(dotResult)<=smallest)=0; 
EnergyTable.fadh2Produced= sum(dotResult,1)';


EnergyTable.OxidativePhosEnergy=v(oxdPhosID,:)';
EnergyTable.EnergyUsedForBiomassGrowth=abs(model.S(atpID,biomassID))* v(biomassID,:)';
EnergyTable.EnergyUsedForMaintenance= EnergyTable.EnergyProduced-EnergyTable.EnergyUsedForBiomassGrowth;
EnergyTable.MaintenanceFracn= EnergyTable.EnergyUsedForMaintenance./ EnergyTable.EnergyProduced;
EnergyTable.MaintenanceEnergyPerUnitBiomass=EnergyTable.EnergyUsedForMaintenance./v(biomassID,:)'; % mols of ATP maintenance per g of cell

end