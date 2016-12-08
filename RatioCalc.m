function [ Ratio ] = RatioCalc( mmets1, mmets2,num,FluxMatrix,S,i )
%RatioCalc determines important flux ratios in e coli central metabolism
        if isempty(mmets2) 
            Num=abs(S(mmets1,num)')* FluxMatrix(num,i);

        else
           Num1=find(((S(mmets1,:)<0)' & (S(mmets2,:)>0)' & FluxMatrix(:,i)>0) ...
            | ((S(mmets1,:)>0)' & (S(mmets2,:)<0)' & FluxMatrix(:,i)<0));
        
           Num=dot(abs(S(mmets1,Num1)'),FluxMatrix(Num1,i));
        
        end
        
        dem=((S(mmets1,:)>0)' & FluxMatrix(:,i)>0)  |  ...
                ((S(mmets1,:))'<0 & FluxMatrix(:,i)<0); 
        Dem= dot(abs(S(mmets1,dem)'),FluxMatrix(dem,i));

        Ratio=Num/Dem;
        if isnan(Ratio),Ratio=0; end % when glucose is not carbon source
        
end

