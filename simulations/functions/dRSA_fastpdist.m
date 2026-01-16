function D = dRSA_fastpdist(X, metric)
    % X: m-by-n data matrix
    % metric: 'euclidean' or 'correlation'
    % Output: condensed distance vector in pdist order
    
    [m, ~] = size(X);
    numPairs = m*(m-1)/2;
    D = zeros(1, numPairs);
    
    switch lower(metric)
%         case 'euclidean'
%             % Pairwise loop for Euclidean
%             k = 1;
%             for i = 1:m-1
%                 for j = i+1:m
%                     xi = X(i,:);
%                     xj = X(j,:);
%                     valid = ~isnan(xi) & ~isnan(xj);
%                     xi_valid = double(xi(valid));
%                     xj_valid = double(xj(valid));
% 
%                     if isempty(xi_valid) || isempty(xj_valid)
%                         d = NaN;
%                     else
%                         d = sqrt(sum((xi_valid - xj_valid).^2));
%                     end
%                     D(k) = d;
%                     k = k + 1;
%                 end
%             end
            
        case 'correlation'
            % Vectorized correlation using corr + mean
            Xc = double(X) - mean(double(X),2);  %center each row
            C = corr(Xc');            % correlation between rows
            D = squareform(1 - C);    % condensed vector in pdist order
            
        otherwise
            error('Unsupported metric. Use ''euclidean'' or ''correlation''.');
    end
end
