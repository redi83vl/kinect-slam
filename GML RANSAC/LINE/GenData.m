function [Data, Model, Stat] = GenData(varargin)

g_NumOfPoints = 500;
g_ErrPointPart = 0.4;  % This points are impulse errors
%g_AbsPulseErr = 30;    % Variation of Bad points
g_NormDistrVar = 1;    % Standard deviation of Good points

if nargin == 1 && isnumeric(varargin{1})
    g_ErrPointPart = varargin{1};
end


        % Generaiting Random Data
        theta = (rand(1) + 1) * pi/6;
        
        R = ( rand([1 g_NumOfPoints]) - 0.5) * 100;
        %R= sort(R);
        Dist = randn([1 g_NumOfPoints]) * g_NormDistrVar;
        %Pulse = zeros([1 g_NumOfPoints]);
        %Pulse(1:floor(g_ErrPointPart * g_NumOfPoints)) = (rand([1 floor(g_ErrPointPart * g_NumOfPoints)]) -0.2 )* g_AbsPulseErr;
        
        Data = [cos(theta); sin(theta)] * R + [-sin(theta); cos(theta)] * Dist;
        Data(:, 1:floor(g_ErrPointPart * g_NumOfPoints)) = 2 * [max(abs(Data(1,:))), 0; 0, max(abs(Data(2,:)))] *...
                                                        (rand([2 floor(g_ErrPointPart * g_NumOfPoints)]) - 0.5);

        Model = [sin(theta) -cos(theta) 0];        
        RealLineMod = [cos(theta); sin(theta)] * sort(R);
        
        % Reconstructing
        RecLine = TLS(Data);
        RecLine = - RecLine([1 3]) ./ RecLine(2);
        
        RecLineMod = [min(Data(1, :)), max(Data(1, :)); 0 , 0];
        RecLineMod(2, :) = RecLineMod(1, :) * RecLine(1) + RecLine(2);
        
%         % Drawing
%         clf;
%         hold on;
%         plot(RealLineMod(1, :), RealLineMod(2, :), '--r');
%         plot(RecLineMod(1, :), RecLineMod(2, :), ':g');
%         
%         plot(Data(1, :), Data(2, :), '.b', 'MarkerSize', 16);        
%         hold off;
%         
%         zoom on;
%         figure(gcf);

        % Statistics
        DIST = Model * [Data; ones([1 g_NumOfPoints])];
        
        Stat.gen_inliers = ones( [1 g_NumOfPoints] );
        Stat.gen_inliers( 1 : floor( g_ErrPointPart * g_NumOfPoints ) ) = 0;
        Stat.gen_mean = mean( DIST( find( Stat.gen_inliers ) ) );
        Stat.gen_abs_mean = mean( abs( DIST( find( Stat.gen_inliers ) ) ) );
        
        Stat.res_inliers = double(abs(DIST) < 2 * g_NormDistrVar);
        Stat.res_mean = mean( DIST( find( Stat.res_inliers ) ) );
        Stat.res_abs_mean = mean( abs( DIST( find( Stat.res_inliers ) ) ) );
