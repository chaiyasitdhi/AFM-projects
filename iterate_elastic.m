function [best_elastic]  = iterate_elastic(depth,force,lower_bound,upper_bound,increment)

global v alpha

midpoint = (lower_bound + upper_bound)/2;
direction = 1; 
local_min = midpoint;
original_trial = [local_min v alpha];
y_original_model = hertz_evaluate(depth,original_trial);
y_measured = force;

% Chi-square error criterion
weight = length(depth)/sqrt(force'*force);
delta_y = ( y_measured - y_original_model ) ;
previous_error = (delta_y'*delta_y)/weight;
previous_value = local_min;

while 1 == 1
trial_value = local_min + increment*(direction);
p_trial = [trial_value v alpha];

% Model
y_model = hertz_evaluate(depth,p_trial);

% Chi-square error criterion
weight = length(depth)/sqrt(force'*force);

delta_y = ( y_measured - y_model ) ;
trial_error = (delta_y'*delta_y)/weight;

    if trial_value - previous_value <= 0
        best_elastic = trial_value;
        break
    elseif trial_error > previous_error
         direction = -1;
         
    elseif trial_error < previous_error
         direction = 1;
         local_min = trial_value;
         previous_error = trial_error;
         previous_value = trial_value;
    end
end

% iterate_contact_point
while 1 == 1
    
end
end
