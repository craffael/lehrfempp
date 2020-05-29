function convergenceH

N = [39; 125; 444; 1670];

meshsizes = [3.52544; 2.22684; 1.24263; 0.669219];

eL2 = [0.555827; 0.370117; 0.20476; 0.106244];

% Normalize L2 error obtained from simulation

 eL2 = (1./sqrt(N)) .* eL2;
 
 p = polyfit(log(N), log(eL2), 1);         % result: p = [-0.9441, 1.0910]
 
 figure()
 loglog(N, eL2, '-o')
 hold on
 loglog(N, exp(p(1)*log(N) + p(2)))
 xlabel('Number of degrees of freedom')
 ylabel('L2 error with respect to solution on finest mesh');
 legend('L2 error', 'slope: -0.94');
 
 