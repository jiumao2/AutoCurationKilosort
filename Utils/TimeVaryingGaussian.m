classdef TimeVaryingGaussian < handle
    properties
        mu % Mean vector
        sigma % Variance vector
        t_bins % Time vector
    end
    
    methods
        function obj = TimeVaryingGaussian()
            obj.mu = [];
            obj.sigma = [];
            obj.t_bins = [];
        end
        
        function fit(obj, data, times_sec, binwidth_sec)
            if nargin < 4
                binwidth_sec = 180;
            end

            idx_outliers = isoutlier(data);
            data = data(~idx_outliers);
            times_sec = times_sec(~idx_outliers);

            t_min = min(times_sec);
            t_max = max(times_sec);
            t_edges = t_min:binwidth_sec:t_max;
            obj.t_bins = 0.5*(t_edges(1:end-1)+t_edges(2:end));

            % Fit the time-varying Gaussian model to the data
            % Estimate the mean 
            obj.mu = zeros(1,length(obj.t_bins));
            data_new = data;
            for k = 1:length(obj.t_bins)
                idx_this = times_sec>t_edges(k) & times_sec<t_edges(k+1);
                data_this = data(idx_this);
                if isempty(data_this)
                    continue
                end
                [f,xi] = ksdensity(data_this);
                [~, idx_max] = max(f);
                obj.mu(k) = xi(idx_max);

                data_new(idx_this) = data_this - obj.mu(k);
            end

            % Estimate the variance
            [f,xi] = ksdensity(data_new);
            x_min_up = max(xi);
            x_min_down = min(xi);

            [~, idx0] = min(abs(xi));
            for k = idx0:length(xi)-1
                if f(k+1)-f(k) > 0
                    x_min_up = xi(k);
                    break
                end
            end

            for k = idx0:-1:2
                if f(k-1)-f(k) > 0
                    x_min_down = xi(k);
                    break
                end
            end

            obj.sigma = std(data_new(data_new>x_min_down & data_new<x_min_up));
            
%             obj.sigma = std(data_new);
        end
        
        function y = pdf(obj, x, t)
            % Evaluate the PDF at time t
            mu_t = interp1(obj.t_bins, obj.mu, t, 'linear', 'extrap');
            y = normpdf(x-mu_t, 0, obj.sigma);
        end
    end
end
