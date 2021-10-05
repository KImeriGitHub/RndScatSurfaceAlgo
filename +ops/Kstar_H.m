classdef Kstar_H < ops.Operators
	% Adjoint of the Neumann-Poincare operator for Helmholtz equation
	% H means  helmholtz
	methods
		function obj = Kstar_H(k, D1, type1, step1, D2, type2, step2)
            if nargin < 7 % If only one boundary is given
				D2 = D1;
				type2 = type1;
				step2 = step1;
            end 
            
			obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
			
            if isequal(D1,D2)
                obj.Kmat = ops.Kstar_H.make_kernel_matrix(k, D1.points, D1.tvec, D1.normal, D1.avec, D1.sigma);
            else
                obj.Kmat = ops.Kstar_H.make_kernel_matrix_disjoint(k, D2.points, D2.sigma, D1.points, D1.normal);
            end
            
			
		end
		
	end
	
	methods(Static)
		function KsH = make_kernel_matrix(k, D, tvec, normal, avec, sigma)
			% The discretization of the adjoint operator of K using P0 boundary
			% element. Precisely, we construct a matrix Ks such that Ks*phi
			% approximates the integral definition of K^*[phi].
			
			M = size(D,2);
			KsH = zeros(M, M);
            KsHx = zeros(M, M);
			tvec_norm_square = tvec(1,:).^2 + tvec(2,:).^2;
			
			for j = 1:M
				xdoty = (D(1,j)-D(1,:))*normal(1,j)+(D(2,j)-D(2,:))*normal(2,j);
				norm_xy = sqrt((D(1,j)-D(1,:)).^2+(D(2,j)-D(2,:)).^2);
				
                KsH(j,1:j-1) =  1i/4*k*besselh(1,1,k*norm_xy(1:(j-1))).*sigma(1:(j-1)).* xdoty(1:(j-1))./norm_xy(1:(j-1));
				KsH(j,j+1:M) = 1i/4*k*besselh(1,1,k*norm_xy((j+1):M)).*sigma((j+1):M).* xdoty((j+1):M)./norm_xy((j+1):M);
% 				KsH(j,1:j-1) = -1i/4*k*besselh(1,2,k*norm_xy(1:(j-1))).*sigma(1:(j-1)).* xdoty(1:(j-1))./norm_xy(1:(j-1));
% 				KsH(j,j+1:M) = -1i/4*k*besselh(1,2,k*norm_xy((j+1):M)).*sigma((j+1):M).* xdoty((j+1):M)./norm_xy((j+1):M);
                KsHx(j,1:j-1) = 1i/4*k*conj(besselh(1,2,conj(k)*norm_xy(1:(j-1)))).*sigma(1:(j-1)).* xdoty(1:(j-1))./norm_xy(1:(j-1));
                KsHx(j,j+1:M) = 1i/4*k*conj(besselh(1,2,conj(k)*norm_xy((j+1):M))).*sigma((j+1):M).* xdoty((j+1):M)./norm_xy((j+1):M);
                KsHx(j,j) = -1/4/pi*avec(:,j)'*normal(:,j)/tvec_norm_square(j)*sigma(j);				
                KsH(j,j) = -1/4/pi*avec(:,j)'*normal(:,j)/tvec_norm_square(j)*sigma(j);
            end
%             KsH = conj(KsH);
        end
        
        function KsH = make_kernel_matrix_disjoint(k, D, D_sigma, E, E_normal)
            %This is for the 'cross' matrices
            M = size(D,2);
			KsH = zeros(M, M);
            
            for j = 1:M
                for l = 1:M
% 				xdoty = (E(1,j)-D(1,:))*E_normal(1,j)+(E(2,j)-D(2,:))*E_normal(2,j);
% 				norm_xy = sqrt((E(1,j)-D(1,:)).^2+(E(2,j)-D(2,:)).^2);
                xdoty = (E(1,j)-D(1,l))*E_normal(1,j)+(E(2,j)-D(2,l))*E_normal(2,j);
                norm_xy = sqrt((E(1,j)-D(1,l)).^2+(E(2,j)-D(2,l)).^2);
                KsH(j,l) = 1i/4*k*besselh(1,1,k*norm_xy)*D_sigma(l)*xdoty/norm_xy;
                end
            end
        end
		
		function val = eval()
			error('Method not implemented!');
		end
		
	end
end

