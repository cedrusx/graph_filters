function RMSE = rmse(A, B)
% Root mean square error
err = A - B;
RMSE = sqrt(mean(err(~isnan(err)).^2));
