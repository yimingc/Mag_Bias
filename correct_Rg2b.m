function [R_plus,err_flag] = correct_Rg2b(R_minus, d_angle)
I3 = eye(3);
if dot(d_angle,d_angle) < 1

	% 1st order approximation
	R_plus = R_minus * (I3 - vcross(d_angle));

	% normalize
	R_plus = (I3 - (1/2)*(R_plus*R_plus' - I3))*R_plus;
    err_flag = 1;
else
	fprintf(1, 'invalid Rb2t correction |varphi|=%g\n', norm(d_angle));
	R_plus = R_minus;
    err_flag = 0;
end