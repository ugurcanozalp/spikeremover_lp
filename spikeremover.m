function [ x ] = spikeremover( z, anomaly_ratio, window_length)
% This function cleans data using linear programming
% z (Nx1 vector): original data
% anomaly_ratio [0~1]: expected ratio of anomaly number to all data 
% (typical value: 0.03)
% window_length [3~N]: maximum expected anomaly duration
% x : cleaned data
% Remark: low anomaly_ratio and window_length accelerates the function. So
% dont pass them unnecessarily high
% Required functions: linprog2
% example : [ x ] = spikeremover( z, 0.04,  50);
% ------------------------------------------------------------------------
% This function fits smoothed signal into spiked signal such that
% derivative of fitted signal never exceeds maximum derivative specified by
% anomaly_ratio. It is done by linear programming.
% Author    : Ugurcan Ozalp
% Mail      : uurcann94@gmail.com / Please report any problem you face.

% Change following line to use other optimizers.
OPTIMIZER = 2; % 1: linprog function (MATLAB Built-in Linear Programming)
               % 2: linprog2 function
dz = diff(z); % Derivative of signal at each point.
n = length(z);
% Reduce window_length if necessary to accelerate algorithm.
if anomaly_ratio*n < window_length 
    window_length = floor(anomaly_ratio*n);
end
% std_dz = std(diff(z));
% From number of point, data is suspected to have spike ? -> index_th
index_th = floor((1-anomaly_ratio)*n); 
% Sort derivative values and find indices.
[sorted_abs_dz,dz_ind] = sort(abs(dz));
% Maximum allowed derivative value.
dz_lim = sorted_abs_dz(index_th);
% Which points are going to be changed. Only these points and surrounding 
% points (with span of window length) will be optmized.
act_const = dz_ind(index_th:end);
x_p_ind_oz = false(n,1); % Which points will be optimized? At first none 
% of them, in for loop they will be filled.
windowd2 = ceil(window_length/2); % Half of window length.
for i = 1:length(act_const) % ind_tmp = act_const'
    ind_tmp = act_const(i);
    st = max(1,ind_tmp-windowd2);
    ed = min(ind_tmp+windowd2,n-1);
    x_p_ind_oz(st:ed) = true; % Which points will be optimized?
end
x_p_ind = find(x_p_ind_oz); % Which points will be optimized, as index.
n_const = length(x_p_ind); % Number of constraints
%% Linear Programming Formatting.
% Convert information into linear programming format.
% sigma = x-z (Error to be optimized) -> x = sigma+z;
% sigmasize :n , deltasize: n-1
% Optimization parameters
K = speye(n-1,n); K = K(x_p_ind,:);
L = K - [sparse(n_const,1),K(:,1:end-1)];
LL = L(:,x_p_ind);
I = speye(n_const);
c = [ones(2*n_const,1);zeros(2*n_const,1)];
A = [-LL;LL]*[I -I];
b = -dz_lim - [-L;L]*z ;
A_st = [A , -speye(2*n_const)]; % adding surplus variables
if OPTIMIZER==2
    [alpha,fval,~,~] = linprog2(A_st,b,c);
elseif OPTIMIZER==1
    [alpha,fval] = linprog(c,[],[],A_st,b,zeros(size(c)));
else
    error('WRONG OPTIMIZER INDEX!!');
end
disp(['Cost value: ',num2str(fval)]);
sigmapm = alpha(1:2*n_const);
sigma = sigmapm(1:n_const) - sigmapm(n_const+1:end);
% Difference between x and z
addto = zeros(size(z)); addto(x_p_ind)=full(sigma); 
x = addto + z;
end