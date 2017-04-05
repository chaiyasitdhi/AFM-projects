function contact_point  = contact_estimate(zdist, force)

%% ### Contact-point Estimation ###
% Initially, the cantilever vibrates randomly before reaching the contact 
% point. According to the empirical rule, the 99.999~% of data originated
% from the unknown Gaussian distribution stay inside the 5*standard-
% deviation. Therefore, we set the 5*standard deviation as a threshold 
% for contact point selection, since the unusual large deviation from 
% the aveage deflection indicates that the cantilever finally reaches 
% the surface of the sample.

% ### Paramaters ###
% zdist = z-scanner distance data (X)
% force = force data (Y)
%%

% Transpose data to Nx1 dimensions
if size(force,1) > size(force,2)
else
   force = force';
   zdist = zdist';
end

diff_force = diff(force);
% Find the first 20% of the data
%need to change in some sample 
int_def = force(1:round(length(force)*0.15));

% Find the upper boundary
    upper_bound = 5.*std(int_def)+median(int_def);


index_app = find( force(:,1) >= upper_bound(1) );
contact_point = zdist(index_app);

% plot(1:numel(force),upper_bound.*ones(1,numel(force)));
% hold on; plot(1:numel(force),force, 'ro');

end