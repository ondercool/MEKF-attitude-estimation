function rmse_val = errArr2RMSE(e_arr)
    rmse_val = sqrt(mean(e_arr.^2, 2));
end