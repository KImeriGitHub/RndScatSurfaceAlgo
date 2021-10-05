classdef Spline < shape.C2boundary
    %SPLINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        AddPts %Number of pts between given points
    end
    
    methods
        function obj = Spline(D, sepPtr, AConst)
            % Pre: D is a (2 x N) Array
            %      sepPtr is (1 x M) in N, starts with 1, has distinct values
            %      AConst is a positive integer ( >= 1 )
            
            % sepPtr describes a (1xM) array starting with 1. Elements point to  
            % indices in D, which describe the beginning of a new simply
            % connected domain  
            
            if nargin<3
				AConst = 1; %one mid point
            end
            
            if sepPtr(end)<size(D,2)
               sepPtr=[sepPtr,size(D,2)+1]; 
            end
            
            M=length(sepPtr);
            
            sepMPtr = cumsum([1,diff(sepPtr)*(1+AConst)]);
            points = zeros(2, sepMPtr(end)-1);
            tvec = zeros(2, sepMPtr(end)-1);
            avec = zeros(2, sepMPtr(end)-1);
            
            
            for m=1:M-1
                %linear interpolation of underlying parametrization 'tau'
                tauExt = linspace(0, 2*pi, sepPtr(m+1)-sepPtr(m)+1);
                tau = tauExt(1:end-1);
                interpVal = (0:AConst)/(AConst+1);
                intertau = tau.*(1-interpVal.')+tauExt(2:end).*interpVal.';
                intertau = reshape(intertau, 1, numel(intertau));
                
                
                %Spline interpolation
                DmExt=[D(:, sepPtr(m):sepPtr(m+1)-1), D(:,sepPtr(m))];
                intertauExt=[intertau, 2*pi];
                %interD=spline(tauExt, DmExt, intertauExt);
                interD=pchip(tauExt, DmExt, intertauExt);
                
                
                %Shape properties
                diffInterTauExt=diff(intertauExt);
                tmvec = [diff(interD(1,:))./diffInterTauExt; diff(interD(2,:))./diffInterTauExt];
                tmvecExt = [tmvec, tmvec(:,1)];
                amvec = [diff(tmvecExt(1,:))./diffInterTauExt; diff(tmvecExt(2,:))./diffInterTauExt];
                
                
                %Add to final variables
                points(:,sepMPtr(m):sepMPtr(m+1)-1) = interD(:,1:end-1);
                tvec(:,sepMPtr(m):sepMPtr(m+1)-1) = tmvec;
                avec(:,sepMPtr(m):sepMPtr(m+1)-1) = amvec;
            end
            
            normal = [0 1;-1 0]*tvec;
            normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1);
            
            obj = obj@shape.C2boundary(points, tvec, avec, normal, sepMPtr, 'Spline');
            obj.AddPts = AConst;
        end
    end
end
