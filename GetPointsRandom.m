function [X, Y, D] = GetPointsRandom(nWant, XWidth, YWidth, R)
X       = zeros(nWant, 1);
Y       = zeros(nWant, 1);
dist_2  = R ^ 2;  % Squared once instead of SQRT each time
iLoop   = 1;      % Security break to yoid infinite loop
nValid  = 0;
while nValid < nWant && iLoop < 1e6
  newX = XWidth * (rand - 0.5);
  newY = YWidth * (rand - 0.5);
  if all(((X(1:nValid) - newX).^2 + (Y(1:nValid) - newY).^2) > dist_2)
    % Success: The new point does not touch existing points:
    nValid    = nValid + 1;  % Append this point
    X(nValid) = newX;
    Y(nValid) = newY;
  end
  iLoop = iLoop + 1;
end
% An error is safer than a warning:
if nValid < nWant
  error('Cannot find wanted number of points in %d iterations.', iLoop)
end
% [Edited start] "figure" inserted, 'Parent' used, 'XGrid' inserted:
% FigH  = figure;
% AxesH = axes('XTick', -100:20:100, 'YTick', -100:20:100, ...
%      'XLim', [-100, 100], 'YLim', [-100, 100], ...
%      'XGrid', 'on', 'YGrid', 'on', 'NextPlot', 'add', ...
%      'Parent', FigH);   
% plot(X, Y, 'b*', 'Parent', AxesH);
% [EDITED end]
if nargout > 2
  % D = pdist([X, Y]);   % Faster with statistics toolbox
  D = sqrt(bsxfun(@minus, X, X.') .^ 2 + bsxfun(@minus, Y, Y.') .^ 2);
end
end