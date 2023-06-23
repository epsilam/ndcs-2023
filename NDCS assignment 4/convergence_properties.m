function [q, mu] = convergence_properties(errors)
    num_sequences = size(errors,2);

    q = nan(num_sequences,1); % convergence orders
    mu = nan(num_sequences,1); % convergence rates

    for index = 1:num_sequences
        error_sequence = errors(:,index);
        error_sequence = error_sequence(~isnan(error_sequence));
        half_samples = round(length(error_sequence)/2);
        error_sequence = error_sequence(half_samples:end);

        errors_0 = error_sequence(1:end-2); % No shift
        errors_1 = error_sequence(2:end-1); % Shift by 1
        errors_2 = error_sequence(3:end);   % Shift by 2
    
        % Convergence order
        q(index) = round(mean(log(errors_2 ./ errors_1)./log(errors_1./errors_0)));
    
        % Convergence rate
        error_q_ratios = errors_1 ./ errors_0.^q(index);
        
        mu(index) = mean(error_q_ratios);
    end
end