classdef C2boundary
	% Abstract class for C2-smooth closed boundary.
	
    properties(SetAccess = protected)
        %% Manually set
        
        points % coordinates of boundary points, an array of dimension 2 X NPoints
        tvec % tangent vector
        avec % acceleration vector
        normal % outward normal vector
        seperationPtr % (1xM) array starting with 1. Elements point to indices in previous
        % variables, which describe the beginning of a new simply connected
        % domain
        name_str % name of the class
        
        
    end
    properties (Dependent)
        %% Automatically set
        
        theta % non tied-off parameterization between [0, 2pi) for every gement. (1 x N) array
        cpoints % complexification of the boundary points: points(1,:)+1i*points(2,:)
        diameter % (upper bound of) diameter of the shape
        tvec_norm % norm of the tangent vector
        sigma % element of curve integration. (1 x N) array
        center_of_mass % center of mass
        % \int f(x) ds(x) = \int_0^(2pi) f(x(t)) |x'(t)| dt
        % ~ \sum_{n=1..NPoints} f(points(n)) * sqrt(tvec_norm_square(n)) * 2pi/NPoints
        % = \sum_{n=1..NPoints} f(points(n)) * sigma(n)
        NPoints % number of discrete boundary points (for runtime compatibility with the
        % functions like get.theta)
        
        box % a minimal rectangular box [width, height] containing the shape (width: size in x-axis, height: size in y-axis)
        flag % a validity flag for the curve to be well defined
    end
	
	methods
		function obj = C2boundary(points, tvec, avec, normal, sepPtr, nstr)
			% Initialization of C2boundary object:
			% Inputs:
			% points: coordinates of boundary point, a 2 X NPoints array
			% tvec: tangent vectors of boundary point, a 2 X NPoints array
			% avec: acceleration vectors, a 2 X NPoints array
			% normal: outward normal vectors
            % sepPtr: (1 x M) in N, starts with 1, has distinct values
			% nstr: name of the boundary, a string, optional
			
			obj.points = points;
			obj.tvec = tvec;
			obj.avec = avec;
			obj.normal = normal;
			obj.seperationPtr = sepPtr;
            
			% Check the curve
% 			flag = shape.C2boundary.check_sampling(points);
%             if ~flag
%                 warning('Curve may contain singularities.');
%             end
			
			if nargin > 5
				obj.name_str = nstr;
			else
				obj.name_str = 'EmptyName';
			end
        end
        
		function val = get.theta(obj)
            sepPtr=obj.seperationPtr;
            val=zeros(1,sepPtr(end)-1);
            for m=1:length(sepPtr)-1
                temp = linspace(0,2*pi,sepPtr(m+1)-sepPtr(m)+1); % non tied-off
                val(sepPtr(m):sepPtr(m+1)-1) = temp(1:end-1);
            end
			% theta = linspace(0,2*pi,obj.NPoints); % tied-off
        end
        
        function val=get.center_of_mass(obj)
            sepPtr=obj.seperationPtr;
            pts=obj.points;
            nor=obj.normal;
            val=zeros(1,length(sepPtr)-1);
            for m =1:length(sepPtr)-1
                idxVal=sepPtr(m):sepPtr(m+1)-1;
                velocityLength=sqrt(sum(obj.tvec(idxVal).^2));
                thetaDis=sqrt(diff([obj.theta(idxVal),2*pi]).^2).';
                val(m) = [(0.5*(pts(1,idxVal).^2.*nor(1,idxVal).*velocityLength)*thetaDis)...
                    /((pts(1,idxVal).*nor(1,idxVal).*velocityLength)*thetaDis);
                    (0.5*(pts(2,idxVal).^2.*nor(2,idxVal).*velocityLength)*thetaDis)...
                    /((pts(2,idxVal).*nor(2,idxVal).*velocityLength)*thetaDis)];
            end
        end
        
		function val = get.box(obj)
            sepPtr=obj.seperationPtr;
            dd = obj.points;
            val=zeros(2,length(sepPtr)-1);
            for m =1:length(sepPtr)-1
                idxVal=sepPtr(m):sepPtr(m+1)-1;
                w = max(dd(1,idxVal)) - min(dd(1,idxVal));
                h = max(dd(2,idxVal)) - min(dd(2,idxVal));
                val(:,m) = [w; h];
            end
        end
		
		function val = get.NPoints(obj)
			val = diff(obj.seperationPtr);
		end
		
		function val = get.diameter(obj)
            sepPtr=obj.seperationPtr;
            val=zeros(1,length(sepPtr)-1);
            dd = obj.points;
            for m =1:length(sepPtr)-1
                idxVal=sepPtr(m):sepPtr(m+1)-1;
                val(m) = max(sqrt((dd(1,idxVal).'-dd(1,idxVal)).^2 + (dd(2,idxVal).'-dd(2,idxVal)).^2),[],'all');
            end
        end
		
		function val = get.cpoints(obj)
			val = obj.points(1,:)+1i*obj.points(2,:);
		end
		
		function val = get.tvec_norm(obj)
			val = sqrt(obj.tvec(1,:).^2 + obj.tvec(2,:).^2);
		end
		
		function val = get.sigma(obj)
            sepPtr=obj.seperationPtr;
            val=zeros(1,sepPtr(end)-1);
            for m=1:length(sepPtr)-1
                idxVal=sepPtr(m):sepPtr(m+1)-1;
                val(idxVal) = diff([obj.theta(idxVal), 2*pi]) .* obj.tvec_norm(idxVal);
            end
		end
		
		%% Overloading of usual operators
        %TO BE IMPLEMENTETD
		function obj = plus(obj, z0)
            %Pre: z0 is a (2 x 1) double array
			% Overload of the operator +
			if ~isa(z0, 'double')
				error('Type error: only double value can be used for translation.');
            end
			obj.points = [obj.points(1,:)+z0(1); obj.points(2,:)+z0(2)];
		end
% 		
% 		function obj = minus(obj, z0)
% 			% Overload of the operator -
% 			obj = plus(obj,-z0);
% 		end
% 		
% 		function obj = mtimes(obj, s)
% 			% Overload of the operator *
% 			if ~isa(s, 'double') || length(s) ~= 1 || s <= 0
% 				error('Type error: only positive double scalar can be used for scaling.');
% 			end
% 			
% 			obj.points = obj.points * s;
% 			obj.center_of_mass = obj.center_of_mass * s;
% 			obj.tvec = obj.tvec * s;
% 			obj.avec = obj.avec * s;
%         end
		
		
		function plot(obj, Linespec)
            if nargin<2
                sepPtr=obj.seperationPtr;
                hold on;
                for m=1:length(sepPtr)-1
                    idxVal=sepPtr(m):sepPtr(m+1)-1;
                    plot(obj.points(1,idxVal), obj.points(2,idxVal));
                end
                hold off;
            else
                
                
                sepPtr=obj.seperationPtr;
                hold on;
                for m=1:length(sepPtr)-1
                    idxVal=sepPtr(m):sepPtr(m+1)-1;
                    plot(obj.points(1,idxVal), obj.points(2,idxVal), Linespec);
                end
                hold off;
            end
        end
        
        function obj=merge(obj, dom)
            %Pre:  dom is of class C2boundary, or either one is empty
            %          Those domains should be distinct from each other
            %Post: DomRes is of class C2boundary
            %
            %Desc: Extends obj using the domain dom.
            % So especially DomRes.points == [Dom1.points, Dom2.points]
            if isempty(dom)
                return
            end
            obj.points = [obj.points, dom.points];
            obj.tvec = [obj.tvec, dom.tvec];
            obj.avec = [obj.avec, dom.avec];
            obj.normal = [obj.normal, dom.normal];
            obj.seperationPtr = [obj.seperationPtr, obj.seperationPtr(end)-1+dom.seperationPtr(2:end)];
            obj.name_str='mergedDomain';
            
        end
    end
	
	methods(Static)
		[tvec,avec,normal] = boundary_vec(D, t)
		[D, tvec, avec, normal] = boundary_vec_interpl(points0, theta0, theta)
		
		[D, tvec, avec, normal] = rescale(D0, theta0, NPoints, nsize, dspl)
		[D, tvec, avec, normal] = rescale_diff(D0, theta0, NPoints, nsize)
		
		function val = check_sampling(points)
			% A C^1 parameterized simple curve f(t) (t is the parameter) must satisfy (by Taylor expansion)
			% <f(t_n+1)-f(t_n), f(t_n)-f(t_n-1)> > 0, for sufficiently
			% small sampling step dt = t_n - t_n-1.
			%
			% This function checks this condition.
			
			NPoints = size(points, 2);
			
			val = 1;
			
			for p=1:NPoints
				x = points(:,p);
				
				if p==1
					y=points(:,NPoints);
					z=points(:,p+1);
				elseif p==NPoints
					y=points(:,p-1);
					z=points(:,1);
				else
					y=points(:,p-1);
					z=points(:,p+1);
				end
				
				toto = (z-x)'*(x-y)/(norm(y-x)*norm(z-x));
				if toto <= 0
					val = 0;
				end
			end
        end
	end
end

