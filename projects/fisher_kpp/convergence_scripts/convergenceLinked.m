function convergenceLinked

m = [100; 200; 400; 800];

tau = [0.01; 0.005; 0.0025; 0.00125];

N = [39; 125; 444; 1670];

meshsizes = [3.52544; 1.76272; 0.88136; 0.44068];

eL2 = [6.22189; 3.96191; 0.900926; 0.314753];

% Normalize L2 error obtained from simulation

 eL2 = sqrt(meshsizes) .* eL2;
 
 p1 = polyfit(log(m), log(eL2), 1);         % result: p1 = [-2.0052, 11.9153]
 p2 = polyfit(log(N), log(eL2), 1);         % result: p2 = [-1.1111, 6.7093]
 p3 = polyfit(log(tau), log(eL2), 1);       % result: p3 = [2.0052, 11.9153]
 p4 = polyfit(log(meshsizes), log(eL2), 1); % result: p4 = [2.0052, 0.1546]

 figure()
 loglog(m, eL2, '-o')
 hold on
 loglog(m, exp(p1(1)*log(m) + p1(2)))
 xlabel('Number of time steps')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: -2.01');

 figure()
 loglog(N, eL2, '-o')
 hold on
 loglog(N, exp(p2(1)*log(N) + p2(2)))
 xlabel('Number of degrees of freedom')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: -1.11');
 
 figure()
 loglog(tau, eL2, '-o')
 hold on
 loglog(tau, exp(p3(1)*log(tau) + p3(2)))
 xlabel('Timestep sizes')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: 2.01');
 
 figure()
 loglog(meshsizes, eL2, '-o')
 hold on
 loglog(meshsizes, exp(p4(1)*log(meshsizes) + p4(2)))
 xlabel('Meshsizes')
 ylabel('L2 error with respect to solution of smallest time step size and on finest mesh');
 legend('L2 error', 'slope: 2.01');

