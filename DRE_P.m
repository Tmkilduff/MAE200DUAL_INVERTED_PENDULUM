function dP = DRE_P(P,E,A,C,Q1,Q2)
dP = A'*P*E + E'*P*A - E'*P*C'*Q2^(-1)*C*P*E + Q1;
dP = E'^(-1) * dP * E^(-1);
end