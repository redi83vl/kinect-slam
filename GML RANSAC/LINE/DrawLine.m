function h = DrawLine( Line , Data, varargin)

% Determine bounding box coordinates
MinX = min(Data(1, :));
MaxX = max(Data(1, :));
MinY = min(Data(2, :));
MaxY = max(Data(2, :));

CrossMask = sign( Line * ...
    [MinX, MinX, MaxX, MaxX;
     MinY, MaxY, MaxY, MinY;
      1,    1,    1,    1  ]);

% PreInit
i = 1;
LineCoord = [0, 0; 0, 0];


CrossMask = CrossMask * CrossMask(1);
if CrossMask(2) < 0
    LineCoord(1, i) = MinX;
    LineCoord(2, i) = -(Line(3) + Line(1) * MinX) / Line(2);
    CrossMask = -CrossMask;
    i = i + 1;
end

if CrossMask(3) < 0
    LineCoord(1, i) = -(Line(3) + Line(2) * MaxY) / Line(1);
    LineCoord(2, i) = MaxY;
    CrossMask = -CrossMask;
    i = i + 1;
end

if CrossMask(4) < 0
    LineCoord(1, i) = MaxX;
    LineCoord(2, i) = -(Line(3) + Line(1) * MaxX) / Line(2);
    CrossMask = -CrossMask;
    i = i + 1;
end

if CrossMask(1) < 0
    LineCoord(1, i) = -(Line(3) + Line(2) * MinY) / Line(1);
    LineCoord(2, i) = MinY;
    CrossMask = -CrossMask;
    i = i + 1;
end

if i ~= 3
    h = 0;
else
    h = line(LineCoord(1, :), LineCoord(2, :), varargin{:});
end

return;