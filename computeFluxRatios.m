function [FluxRatioData,RatioLabels]= computeFluxRatios(v,model)
% compute important flux ratios

%% determination of important split ratios & other important variables
FluxRatioData=zeros(size(FluxMatrix2,2),10);
RatioLabels={'Fluxintoglycolysis';'FluxintotheEDpathway';...
    'Fluxintothemethylglyoxalpathway';'PEPtopyruvateflux';...
    'PyruvatetoacetylcoenzymeA';'FluxintoTCAcyclefromacetylcoenzymeA';...
    'Fluxintoglyoxlyateshunt';'OxaloacetatetoPEPflux';...
    'Acetatesecretion';'Ethanolsecretion'};
RatioLabels=RatioLabels';
for i=1:size(v,2)
        %#1 Flux into glycolysis
        mmets=find(strcmp(model.mets,'g6p[c]'));
        num=find(strcmp(model.rxns,'PGI'));
        FluxRatioData(i,1) = RatioCalc( mmets,[],num,FluxMatrix2,S,i );
       

        % #2 Flux into the ED pathway
        mmets=find(strcmp(iAF1260.mets,'6pgc[c]'));
        num=find(strcmp(iAF1260.rxns,'EDD'));
        FluxRatioData(i,2) = RatioCalc( mmets,[],num,FluxMatrix2,S,i );
        

        % #3 Flux into the methylglyoxal pathway
        mmets=find(strcmp(iAF1260.mets,'dhap[c]'));
        num=find(strcmp(iAF1260.rxns,'MGSA'));
        FluxRatioData(i,3) = RatioCalc( mmets,[],num,FluxMatrix2,S,i );

        % #4 PEP to pyruvate flux
        mmets1=find(strcmp(iAF1260.mets,'pep[c]'));
        mmets2=find(strcmp(iAF1260.mets,'pyr[c]'));
        FluxRatioData(i,4) = RatioCalc( mmets1,mmets2,[],FluxMatrix2,S,i );

        % #5 Pyruvate to acetylcoenzyme A
        mmets1=find(strcmp(iAF1260.mets,'pyr[c]'));
        mmets2=find(strcmp(iAF1260.mets,'accoa[c]'));
        FluxRatioData(i,5) = RatioCalc( mmets1,mmets2,[],FluxMatrix2,S,i );

        % #6 flux into TCA cycle from acetylcoenzyme A
        mmets1=find(strcmp(iAF1260.mets,'accoa[c]'));
        mmets2=find(strcmp(iAF1260.mets,'cit[c]'));
        FluxRatioData(i,6) = RatioCalc( mmets1,mmets2,[],FluxMatrix2,S,i );

        % #7 flux into glyoxlyate shunt
        mmets=find(strcmp(iAF1260.mets,'icit[c]'));
        num=find(strcmp(iAF1260.rxns,'ICL'));
        FluxRatioData(i,7) = RatioCalc( mmets,[],num,FluxMatrix2,S,i );


        % #8 Oxaloacetate to PEP flux
        mmets=find(strcmp(iAF1260.mets,'oaa[c]'));
        num=find(strcmp(iAF1260.rxns,'PPCK'));
        FluxRatioData(i,8) = RatioCalc( mmets,[],num,FluxMatrix2,S,i );

        % #9 Acetate secretion
        mmets=find(strcmp(iAF1260.mets,'accoa[c]'));
        num=find(strcmp(iAF1260.rxns,'PTAr'));
        FluxRatioData(i,9) = RatioCalc( mmets,[],num,FluxMatrix2,S,i );
       

        % #10 Ethanol secretion
        mmets=find(strcmp(iAF1260.mets,'accoa[c]'));
        num=find(strcmp(iAF1260.rxns,'ACALD')); % note the sign switch
        FluxRatioData(i,10) = RatioCalc( mmets,[],num,FluxMatrix2,S,i );
        
        disp(i)
        
end

FluxRatioData= array2table(FluxRatioData,'VariableNames',RatioLabels);



