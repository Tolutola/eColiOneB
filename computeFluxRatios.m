function [FluxRatioData,RatioLabels]= computeFluxRatios(v,model)
% compute important flux ratios

%% determination of important split ratios & other important variables
FluxRatioData=zeros(size(v,2),10);
RatioLabels={'Fluxintoglycolysis';'FluxintotheEDpathway';...
    'Fluxintothemethylglyoxalpathway';'PEPtopyruvateflux';...
    'PyruvatetoacetylcoenzymeA';'FluxintoTCAcyclefromacetylcoenzymeA';...
    'Fluxintoglyoxlyateshunt';'OxaloacetatetoPEPflux';...
    'Acetatesecretion';'Ethanolsecretion'};
RatioLabels=RatioLabels';
for i=1:size(v,2)
        %#1 Flux into glycolysis
        mmets=find(strcmp(model.mets,'g6p_c'));
        num=find(strcmp(model.rxns,'PGI'));
        
        FluxRatioData(i,1) = RatioCalc( mmets,[],num,v,model.S,i );
       

        % #2 Flux into the ED pathway
        mmets=find(strcmp(model.mets,'6pgc_c'));
        num=find(strcmp(model.rxns,'EDD'));
        FluxRatioData(i,2) = RatioCalc( mmets,[],num,v,model.S,i );
        

        % #3 Flux into the methylglyoxal pathway
        mmets=find(strcmp(model.mets,'dhap_c'));
        num=find(strcmp(model.rxns,'MGSA'));
        FluxRatioData(i,3) = RatioCalc( mmets,[],num,v,model.S,i );

        % #4 PEP to pyruvate flux
        mmets1=find(strcmp(model.mets,'pep_c'));
        mmets2=find(strcmp(model.mets,'pyr_c'));
        FluxRatioData(i,4) = RatioCalc( mmets1,mmets2,[],v,model.S,i );

        % #5 Pyruvate to acetylcoenzyme A
        mmets1=find(strcmp(model.mets,'pyr_c'));
        mmets2=find(strcmp(model.mets,'accoa_c'));
        FluxRatioData(i,5) = RatioCalc( mmets1,mmets2,[],v,model.S,i );

        % #6 flux into TCA cycle from acetylcoenzyme A
        mmets1=find(strcmp(model.mets,'accoa_c'));
        mmets2=find(strcmp(model.mets,'cit_c'));
        FluxRatioData(i,6) = RatioCalc( mmets1,mmets2,[],v,model.S,i );

        % #7 flux into glyoxlyate shunt
        mmets=find(strcmp(model.mets,'icit_c'));
        num=find(strcmp(model.rxns,'ICL'));
        FluxRatioData(i,7) = RatioCalc( mmets,[],num,v,model.S,i );


        % #8 Oxaloacetate to PEP flux
        mmets=find(strcmp(model.mets,'oaa_c'));
        num=find(strcmp(model.rxns,'PPCK'));
        FluxRatioData(i,8) = RatioCalc( mmets,[],num,v,model.S,i );

        % #9 Acetate secretion
        mmets=find(strcmp(model.mets,'accoa_c'));
        num=find(strcmp(model.rxns,'PTAr'));
        FluxRatioData(i,9) = RatioCalc( mmets,[],num,v,model.S,i );
       

        % #10 Ethanol secretion
        mmets=find(strcmp(model.mets,'accoa_c'));
        num=find(strcmp(model.rxns,'ACALD')); % note the sign switch
        FluxRatioData(i,10) = RatioCalc( mmets,[],num,v,model.S,i );
        
        disp(i)
        
end
FluxRatioData= array2table(FluxRatioData,'VariableNames',RatioLabels);
