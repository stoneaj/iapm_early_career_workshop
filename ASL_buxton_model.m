%% Function for buxton model
function dM = buxton_pcasl(x, input)
    % Full buxton PCASL model:
    % Buxton et al. (2008). Magnetic Resonance in Medicine, 40:383-396.
    % x(1) = cbf
    % x(2) = aat
    % ti = inversion time ... tagging duration + postlabelling delay

    % PARAMETERS
    % x(1) = CBF
    % x(2) = AAT

    tau1 = x(2); % AAT
    tau2 = x(2) + input.tau;

    for PLD_ID = 1:length(input.PLD)

    	%k = (1/t1b) - (1/t1d);
    	if ( input.PLD(PLD_ID) < tau1)
            dM(PLD_ID) = 0.0;
            %fprintf('1st loop \n')
        elseif ((input.PLD(PLD_ID) >= tau1) && (input.PLD(PLD_ID) <= tau2 ))
            qta = 1-exp(-(input.PLD(PLD_ID)-tau1)/input.t1d);
    		dM_temp = 2 * x(1) * input.m0b * input.t1d * input.alpha * exp(-tau1/input.t1b) * qta;
    		dM(PLD_ID) = dM_temp/6000;
            %fprintf('2nd loop \n')
        elseif (input.PLD(PLD_ID) > tau2)
    		qtb = 1-exp(-input.tau/input.t1d);
    		dM_temp = 2 * x(1) * input.m0b * input.t1d * input.alpha * exp(-tau1/input.t1b) * exp(-(input.PLD(PLD_ID)-input.tau-tau1)/input.t1d) * qtb;
    		dM(PLD_ID) = dM_temp/6000;
            %fprintf('3rd loop \n')
        end
    end

end
