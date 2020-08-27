%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: CRB calculation
Date: Oct./2016
Author: Zhou Zhou

%}

function CRB = CRB( theta_1,theta_2,theta_3,alpha,Q,P,S,sigma_2 )


FIM=Fisher_Information( theta_1,theta_2,theta_3,alpha,Q,P,S,1 );
FIM=FIM/sigma_2;
CRB=inv(FIM);
end

