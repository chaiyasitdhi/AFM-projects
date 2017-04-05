function [sE sR2] = stackCPD(X_1,Y_1, f_1, guess_1,stack)

%fit the entire f-d curve
if stack(1) == 0
    total_range = X_1(end);
    step_size = total_range;
else 
    step_size = stack(2); %nm
end

final_x = X_1(end);
group_x = floor( final_x/step_size );
member_x = floor( length(X_1)/group_x )-1;

range_x = 1;
itr = 1;
    while range_x + member_x <= length(X_1)
        XX = X_1(range_x:(range_x + member_x));
        YY = Y_1(range_x:(range_x + member_x));
         try
            Q = nlinfit(XX,YY,f_1,guess_1);   % perform the non-linear fit
         catch ME
            if strcmp(ME.identifier, 'stats:nlinfit:NoUsableObservations')
                display('error');
                sE = NaN;
                sR2 = NaN;
            end
         end

        E_fit_2 = Q(1)*10^6;
        E = E_fit_2;    % E is the final reported elastic modulus [=] kPa

        fit_val_1 = feval(f_1,E_fit_2/10^6,XX);
        % % Calculate an R^2 value for the fit
        SSres = sum((YY - fit_val_1).^2);
        SStot = (length(YY)-1)*var(YY);
        R2 = 1 - SSres/SStot;

        sE(itr) = E;
        sR2(itr) = R2;
        
        range_x = range_x + member_x;
        itr = itr + 1;
    end
 
end
