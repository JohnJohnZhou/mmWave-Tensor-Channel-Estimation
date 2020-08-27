function [error_theta_1,error_theta_2,error_theta_3,error_alpha,error_H]=parameter_result(est_H,est_alpha,est_theta_1,est_theta_2,est_theta_3,H,alpha,theta1,theta_2,theta_3)

error_theta_1=norm(angular_difference(sort(est_theta_1),sort(theta1)))^2;

error_theta_2=norm(angular_difference(sort(est_theta_2),sort(theta_2)))^2;

error_theta_3=norm(angular_difference(sort(est_theta_3),sort(theta_3)))^2;


[~,I]=sort(abs(est_alpha));
[~,I2]=sort(abs(alpha));
error_alpha=norm(alpha(I2)-est_alpha(I))^2;

error_H=norm(est_H-H)^2;

    function c=angular_difference(a,b)
        l=length(a);
        c=zeros(l,1);
        for jj=1:l
            if abs(a-b)>pi
                c(jj)=abs(a(jj)-b(jj))-pi;
            else
                c(jj)=abs(a(jj)-b(jj));
            end
        end
    end
end