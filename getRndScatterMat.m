function W = getRndScatterMat(N, n, NDisks, k, NFourier)
%PRE: N, n, NDisks are integers
%   k is a real positive number
%   NFourier is of the form 2^(integer)
%POST: W is a complex N x n matrix
%DESC: We apply a scattering in a domain as shown below and measure
%intensities at randomly chosen points in the domain and build up a matrix W.
%
%    *     S        *           |   *       R           *    
%    S     *           *        |        R          *         
%     *          S        *     |               R              
% S                *            |          *         R       
%             *         S       |    *                      R  
%--------------------------------------------------------------
%         *   S              *  |    R  *                      
%    *               S          |              R                    
%    S   *         *            |      R              *        
%             S         *       |*                       R    
%
% We have n Receivers and N Sources.

% Ideally N <(2-(NDisks/2)*r_max^2*pi)/(lambda/4)^2;  (lambda = 2pi/k)
% Ideally n <(2-(NDisks/2)*r_max^2*pi)/(lambda/4)^2;  (lambda = 2pi/k)
% ideally r_end_max*4 < 2*pi/k  (rmax<lambda/4)

%% Internal Vital Parameters
r_init = 0.001;
r_max = 0.02; 
r_min = 0.01;
r_amount = 15; % Levels of radii between r_min and r_max and between r_init and r_min.

Disk2DiskThresh=2*(sqrt(1/4/NDisks)-r_max); % Numerically stable if dis>r/4 (experience) 
Pt2PtThresh = 0.02; % ideally dis>lambda/4 (but often not possible. We need N < 2/Pt2PtThresh^2)

%% Build Up Geometry
% Note: While Loops might not end
%Disks
nd_loop = 1;
r_center_Coords_Data = zeros(2,NDisks);
while nd_loop <= NDisks
    newCoord = (rand(2,1)-0.5).*2; % Vector in [-1.0,1.0]^2
    if ~any(sum((r_center_Coords_Data(:,1:nd_loop-1)-newCoord).^2)<(r_max*sqrt(2)*2+Disk2DiskThresh)^2) % if newCoord are close to an existing disk or the origin
        r_center_Coords_Data(:,nd_loop) = newCoord;
        nd_loop=nd_loop+1;
    end
end

r_end_vals = linspace(r_min,r_max,r_amount);
r_end_Data = r_end_vals(randi(r_amount,1,size(r_center_Coords_Data,2)));

% Sources
nr_loop = 1;
z = zeros(2,N);
while nr_loop<=N
    zrnd2 = (rand(1,1)-0.5)*2;
    zrnd1 = rand(1,1)*0.9-1;
    if ~any(sum((r_center_Coords_Data-[zrnd1;zrnd2]).^2)<(r_max*sqrt(2)+Disk2DiskThresh/2)^2) ...
      && ~any(sum((z(:,1:nr_loop-1)-[zrnd1;zrnd2]).^2)<Pt2PtThresh^2)
        z(:,nr_loop) = [zrnd1;zrnd2];
        nr_loop = nr_loop + 1;
    end
end

%Receivers
nr_loop = 1;
x = zeros(2,n);
while nr_loop<=n
    xrnd2 = (rand(1,1)-0.5)*2;
    xrnd1 = (rand(1,1)*0.9+0.1);
    if  ~any(sum((r_center_Coords_Data-[xrnd1;xrnd2]).^2)<(r_max*sqrt(2)+Disk2DiskThresh/2)^2) ...
      && ~any(sum((x(:,1:nr_loop-1)-[xrnd1;xrnd2]).^2)<Pt2PtThresh^2)
        x(:,nr_loop) = [xrnd1;xrnd2];
        nr_loop = nr_loop + 1;
    end
end


%% Computing Disk to Disk Neumann function
disp('Getting BigN');
tic;
BigN = getCircNeumann(r_center_Coords_Data, r_init, r_end_Data, r_amount, NFourier, k);
BigN_Time_Data = toc;
disp(['Elapsed time = ', num2str(BigN_Time_Data)]);

[W, ~, ~, ~, ~, ~, ~, ~, ~,~] = getEvalCircNeum(z, x, BigN, r_center_Coords_Data, r_end_Data, NFourier,k);
%UnbiasedW = (W-sum(W,'all')/(size(W,1)*size(W,2)));
%NormalizedW = UnbiasedW/sqrt(var(UnbiasedW,0,'all'));
end

%figure(1); histogram(real(W(:)), 6, 'Normalization', 'probability')
%figure(2); histogram(imag(W(:)), 6, 'Normalization', 'probability')