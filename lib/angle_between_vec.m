function angle=angle_between_vec(u,v)
% the angle between two vectors
% input
% u,v - [N x 3] vector
% output
% angle - [N x 1] angle between the vectors in radians

% this calculates a more numericaly stable answer
% https://stackoverflow.com/questions/16966325/python-precision-failure-in-computing-acute-angle/16996138#16996138

% TODO
% benchmark other approaches
% http://www.cs.berkeley.edu/~wkahan/MathH110/Cross.pdf


angle=atan2(vecnorm(cross(u,v,2),2,2),dot(u,v,2));
end