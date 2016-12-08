function modelNew = addFluxRatio(model, nameRatio,numRxns,denomRxns,numVal,denomVal)
%addFluxRatio adds a flux ratio constraint on a COBRA model 
%similar to addRatioReaction
% modelNew = addRatioReaction(model, listOfRxns, ratioCoeff)
%
%INPUTS
% model         COBRA model structure
% listOfRxns    List of  Reactions
% ratioCoeff    Array of ratio coefficients 
%
%OUTPUT
% modelNew      COBRA model structure containing the ratio
%
% Tola Oyetunde 12/05

modelNew = model;

[rows, ~] = size(model.S);

[~, LocNum] = ismember(numRxns,model.rxns);
[~, LocDenom] = ismember(denomRxns,model.rxns);
modelNew.S(rows+1,:) = 0;
modelNew.S(rows+1,LocNum) = denomVal;
modelNew.S(rows+1,LocDenom) = -numVal;
modelNew.b(rows+1) = 0;
modelNew.mets{rows+1} = nameRatio;
modelNew.metName{rows+1} = nameRatio;
if isfield(modelNew,'note')
    modelNew.note = strcat(modelNew.note,nameRatio, 'has been set to ',numVal/denomVal,'.');
else
    modelNew.note = strcat(nameRatio,' has been set to ',numVal/denomVal,'.');
end


