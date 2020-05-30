function convergenceLinked

m = [100; 158; 284; 527];

tau = [0.01; 0.00631649; 0.00352475; 0.00189825];

N = [39; 125; 444; 1670];

meshsizes = [3.52544; 2.22684; 1.24263; 0.669219];

eL2 = [0.555819420839660; 0.372924; 0.207819; 0.108265];

% Normalize L2 error obtained from simulation

 eL2 = (1./sqrt(N)) .* eL2;
 
 p1 = polyfit(log(m), log(eL2), 1);         % result: p1 = [-2.1103, 7.2942]
 p2 = polyfit(log(N), log(eL2), 1);         % result: p2 = [-0.9390, 1.0734]
 p3 = polyfit(log(tau), log(eL2), 1);       % result: p3 = [-2.112, 7.3036]
 p4 = polyfit(log(meshsizes), log(eL2), 1); % result: p4 = [-2.112, -5.0836]
 
 figure()
 loglog(m, eL2, '-o')
 hold on
 loglog(m, exp(p1(1)*log(m) + p1(2)))
 xlabel('Number of time steps')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: -2.11');

 figure()
 loglog(N, eL2, '-o')
 hold on
 loglog(N, exp(p2(1)*log(N) + p2(2)))
 xlabel('Number of degrees of freedom')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: -0.94');
 
 figure()
 loglog(tau, eL2, '-o')
 hold on
 loglog(tau, exp(p3(1)*log(tau) + p3(2)))
 xlabel('Timestep sizes')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: 2.112');
 
 figure()
 loglog(meshsizes, eL2, '-o')
 hold on
 loglog(meshsizes, exp(p4(1)*log(meshsizes) + p4(2)))
 xlabel('Meshsizes')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: 2.112');

