function est=NonLinearLeastSquares(gps_data,ref_data_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% est=NonLinearLeastSquares(gps_data,s2r)
%
% Function that calculates the single point position solution from GPS 
% pseudo range measurements using a nonlinear least squares approach. 
%
% Input
% 
% gps_data      1*M array of struct with the fields:
%               Satellite - Name of satellite
%               Satellite_Position_NED - Position of the satellite
%               PseudoRange - Measured pseudo ranges
%
% s2r		variance of range measurement error (use ref_data.s2r)
%        
% Output:
%
% est           Struct with the fields:
%               x_h - Matrix where each column holds the estimated position 
%                     and clock offset (meters) for each time instant.
%               P - Matrix where the columns holds the diagonal elements of
%               the state covariance matrix. 
% 
% Author: Isaac Skog, skog@kth.se
% Copyright (c) 2014 KTH, ISC License (open source)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N=length(gps_data(1).PseudoRange);
M=length(gps_data);
est.x_h=zeros(4,N);
est.P=zeros(4,N);



for n=1:N
    
    x=zeros(4,1);
    dx=inf(4,1);
    res=zeros(M,1);
    H=zeros(M,4);
    itr_ctr=0;
    
    %% Nonlinear least squares
    while norm(dx)>0.01 && itr_ctr<10;
        
        for m=1:M
            if ~isnan(gps_data(m).PseudoRange(n))
                dR_h=gps_data(m).Satellite_Position_NED(:,n)-x(1:3);
                res(m)=gps_data(m).PseudoRange(n)-(norm(dR_h)+x(4));
                H(m,1:3)=-dR_h'./norm(dR_h);
                H(m,4)=1;
            else
                res(m)=0;
                H(m,:)=zeros(1,4);
            end
        end
        
        
        % Calculate the correction to the state vector
        dx=(H'*H)\(H'*res);
        
        % Update the state vector
        x=x+dx;
        
        % Update the iteration counter
        itr_ctr=itr_ctr+1;
    end

    % Store the estimate
    est.x_h(:,n)=x;
    est.P(:,n)=ref_data_struct.s2r*diag(inv(H'*H));
end

end


