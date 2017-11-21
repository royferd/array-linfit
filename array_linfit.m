% Performs a least squares weighted fit of a line to data. Reports the 
% offset parameter and slope parameter with uncertainties. 
% author: Roy Ready
% created 2/17/2017 

function f = array_linfit(xdata, ydata, ydata_stdev, numpoints)

    %least squares weighted fit for y = a + bx
    num_files = length(xdata(:,1));
    num_rows = length(xdata(1,:));

    %linear fit arrays
    wxsq = zeros(num_files, num_rows, 1);
    wxy = zeros(num_files, num_rows, 1); 
    wy = zeros(num_files, num_rows, 1);
    wx = zeros(num_files, num_rows, 1);
    wx_sq = zeros(num_files, num_rows, 1);

    %linearly fitted parameter arrays
    weight = zeros(num_files,num_rows,1);
    del = zeros(num_files, 1);
    a = zeros(num_files, 1);
    b = zeros(num_files, 1);
    a_stdev = zeros(num_files, 1);
    b_stdev = zeros(num_files, 1);
    resistance = zeros(num_files, 1);
    resistance_stdev = zeros(num_files, 1);
    num_fitpoints = zeros(num_files, 1);
    fitline = zeros(num_files, 2, num_rows*20);
    
    residual = zeros(num_files,num_rows,1);
    
    for i = 1:num_files
        for j =1:numpoints(i)
            weight(i,j) = 1/(ydata_stdev(i,j)^2);
        end
    end



    for i = 1:num_files
        w_sum = sum(weight(i,1:numpoints(i)));

        for j = 1:numpoints(i)
            wxsq(i,j) = weight(i,j)*(xdata(i,j))^2;
            wxy(i,j) = weight(i,j)*xdata(i,j)*ydata(i,j);
            wy(i,j) = weight(i,j)*ydata(i,j);
            wx(i,j) = weight(i,j)*xdata(i,j);
        end

        wxsq_sum = sum(wxsq(i,1:numpoints(i)));
        wxy_sum = sum(wxy(i,1:numpoints(i)));
        wy_sum = sum(wy(i,1:numpoints(i)));
        wx_sum = sum(wx(i,1:numpoints(i)));
        wx_sum_sq = wx_sum^2;

        del(i) = w_sum*wxsq_sum - wx_sum_sq;
        a(i) = (wxsq_sum*wy_sum - wx_sum*wxy_sum)/del(i);
        b(i) = (w_sum*wxy_sum - wx_sum*wy_sum)/del(i);
        a_stdev(i) = sqrt(wxsq_sum/del(i));
        b_stdev(i) = sqrt(w_sum/del(i));
        resistance(i) = 1/abs(b(i));
        resistance_stdev(i) = 1/b(i)^2*b_stdev(i);

        num_fitpoints(i) = numpoints(i)*20;
        step = max(abs(xdata(i,:)))/num_fitpoints(i); % take range of x values

        for j = 1:num_fitpoints(i)
            fitline(i, j, 1) = j*step;
            fitline(i, j, 2) = a(i) +b(i)*fitline(i, j, 1);
        end

        %calculate residual
        for j = 1:numpoints(i)
            residual(i,j) = ydata(i,j) - a(i) - b(i) * xdata(i,j);
        end
    end
    
%    f = [a, a_stdev, b, b_stdev, fitline, residual];
    f = [a a_stdev b b_stdev];
end